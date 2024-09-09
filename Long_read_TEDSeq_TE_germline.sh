#!/bin/bash
####################################################################################
########################TED-SEQ Long-reads germline pipeline#######################
####################################################################################



####################################################################################
################################Required Tools######################################
####################################################################################
#bbmap (Version 38.18)
#cutadapt (Version 2.6)
#r-dplyr (Version 1.1.3)
#bedtools (Version 2.27.1)
#samtools (Version 1.20)
#minimap2 (Version 2.17)
#seqtk (Version 1.3)
####################################################################################
###################################Settings#########################################
####################################################################################
show_help() {
    echo "TED-SEQ Long-reads germline pipeline"
    echo "Usage: bash Long_read_TEDSeq_TE_germline.sh [options]"
    echo
    echo "Options:"
    echo "  --sample        Name of the sample (required)."
    echo "  --reads         Full path to the raw ONT fastq file (required)."
    echo "  --refDir        Path to reference genome (required)."
    echo "  --outDir        Output directory (required)."
    echo "  --te_family     Desired family analyzed (default: ATCOPIA93)."
    echo "  --barcode       Desired P7 barcode used, indicated as a number from 1 to 48 (required)."
    echo "  --cores         Number of cores to use (default: 20)."
    echo "  --min_cov       Minimum coverage for TE detection."
    echo "  --min_ratio     Minimum ratio (default: 0.8)."
    echo "  --ref_TE        Reference TE bed file (required)."
    echo "  --scriptDir     Directory containing the scripts and the flanking sequences (required)."
    echo "  --help          Show this help message."
    echo
}

# Default values
sample=""
reads=""
refDir=/mnt/data2/leduque/ONT/reference/TAIR10_allChr.fasta
outDir=""
te_family="ATCOPIA93"
barcode=""
CORES=20
min_cov=30
min_ratio=0.8
ref_TE=""
scriptDir=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample) sample=$2; shift 1;;
        --reads) reads=$2; shift 1;;
        --refDir) refDir=$2; shift 1;;
        --outDir) outDir=$2; shift 1;;
        --te_family) te_family=$2; shift 1;;
        --barcode) barcode=$2; shift 1;;
        --cores) CORES=$2; shift 1;;
        --min_cov) min_cov=$2; shift 1;;
        --min_ratio) min_ratio=$2; shift 1;;
        --ref_TE) ref_TE=$2; shift 1;;
        --scriptDir) scriptDir=$2; shift 1;;
        --help) show_help; exit 0;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1;;
    esac
    shift
done

if [[ -z "$sample" || -z "$reads" || -z "$refDir"|| -z "$barcode" || -z "$outDir"|| -z "$ref_TE"|| -z "$scriptDir" ]]; then
    echo "Error: lacking required arguments."
    exit 1
fi

if [[ -z "$outDir" ]]; then
    outDir=/mnt/data7/pvendrell/results_new_pipeline/$sample
fi

mkdir -p $outDir/$sample


# Display chosen options
echo "###############Running TED-SEQ Short-reads germline pipeline################"
echo ""
echo "############################Files###########################"
echo ""
echo "Sample: $sample"
echo "Reads: $reads"
echo "Reference Directory: $refDir"
echo "Output Directory: $outDir/$sample"
echo "TE Family: $te_family"
echo "P7 barcode: $barcode"
echo "Cores: $CORES"
echo "Minimum Coverage: $min_cov"
echo "Minimum Ratio: $min_ratio"
echo "Reference TE File: $ref_TE"
echo "Script Directory: $scriptDir"
echo ""
outTrim=$outDir/$sample/Trimmed
mkdir  $outTrim

echo "############################Starting Pipeline############################"

################################################################################################
################Trimming, demultiplexing and removal of chimeric reads #########################
################################################################################################
echo "obtaining informative reads and PCR duplicate removal"
seqtk seq -r -r $reads > $outTrim/${sample}_r.fq
cat $reads $outTrim/${sample}_r.fq > $outTrim/${sample}_fr.fq
rm  $outTrim/${sample}_r.fq
cutadapt --discard-untrimmed  -g file:$scriptDir/P7_$barcode.fa -e 0.2 -O 20 -o $outTrim/${sample}_trimmed.fq  $outTrim/${sample}_fr.fq --report=minimal -j 8
rm  $outTrim/${sample}_fr.fq
cutadapt --discard-untrimmed -a file:$scriptDir/${te_family}_P5.fa -e 0.2 -O 20 -o $outTrim/${sample}_${te_family}_trimmed.fq  $outTrim/${sample}_trimmed.fq --report=minimal -j 8
rm  $outTrim/${sample}_trimmed.fq
#removal of chimeric reads
cutadapt --discard-trimmed -g file:$scriptDir/all_adaptaters.fa -e 0.2 -O 20 -o $outTrim/${sample}_${te_family}_1st_cleaning.fq  $outTrim/${sample}_${te_family}_trimmed.fq --report=minimal -j 8
cutadapt --discard-trimmed -a file:$scriptDir/all_adaptaters.fa -e 0.2 -O 20 -o $outTrim/${sample}_${te_family}_2nd_cleaning.fq  $outTrim/${sample}_${te_family}_1st_cleaning.fq --report=minimal -j 8

rm $outTrim/${sample}_${te_family}_trimmed.fq 
rm $outTrim/${sample}_${te_family}_1st_cleaning.fq


echo "informative reads obtained and deduplicated"


##### MAPPING reads to the reference genome ##### 
minimap2 -a $refDir $outTrim/${sample}_${te_family}_2nd_cleaning.fq -o $outDir/$sample/${sample}_${te_family}.sam -t 8 -x map-ont
samtools sort $outDir/$sample/${sample}_${te_family}.sam > $outDir/$sample/${sample}_${te_family}.bam
samtools index  $outDir/$sample/${sample}_${te_family}.bam
rm $outDir/$sample/${sample}_${te_family}.sam
echo "mapping done"
echo "TE insertion detection"


#Discard reads with a mapping quality lower than 3, can be lowered if you want to recover insertions in highly repetitive sequences alhtough it may increase the number of false positive insertions
samtools view -h -q 3 --threads $CORES $outDir/$sample/${sample}_${te_family}.bam | samtools sort  --threads $CORES -  > $outDir/$sample/${sample}_clip_disc-local.bam
samtools depth  $outDir/$sample/${sample}_clip_disc-local.bam -d 5000000 >  $outDir/$sample/${sample}_depth.bed 


# get the bedfile 
bamToBed -i  $outDir/$sample/${sample}_clip_disc-local.bam -split | awk ' $1 !="ChrM" && $1 != "ChrC" {print $0}' > $outDir/$sample/${sample}_clip_disc-local_bamtobed.bed

awk '{{if ($6=="-")  print $1"\t"$2+1"\t"$6} if( $6=="+")  print $1"\t"$3"\t"$6}' $outDir/$sample/${sample}_clip_disc-local_bamtobed.bed | awk '{count[$1" "$2" "$3]++} END {for (word in count) print word, count[word]}' | sort -k1,1 -k2,2n > $outDir/$sample/${sample}_stop_site.bed
awk '{{if ($6=="-")  print $1"\t"$3"\t"$6} if( $6=="+")  print $1"\t"$2+1"\t"$6}' $outDir/$sample/${sample}_clip_disc-local_bamtobed.bed | awk '{count[$1" "$2" "$3]++} END {for (word in count) print word, count[word]}' | sort -k1,1 -k2,2n > $outDir/$sample/${sample}_start_site.bed


#R script
# calculate the ratio (nbr_start+nbr_stop)/depth for all position
Rscript $scriptDir/ratio.R $sample $outDir/$sample


# ratio 
# "Chr","start","depth","strand","nbr_start"


awk '{if ($5>=$7) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8} if ($5<=$7) { print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"$7"\t"$8}} ' $outDir/$sample/${sample}_ratio.bed | awk '{if ($4=="+") {print $1"\t"$2"\t"$3"\t""-""\t"$5"\t"$6"\t"$7} if ($4=="-") { print $1"\t"$2"\t"$3"\t""+""\t"$5"\t"$6"\t"$7}} ' > $outDir/$sample/${sample}_TEDseq_ratio_sort.bed
## run this instead if you don't have split reads ie primer far from the beginning of the transposon
awk '{if ($5>=$7) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8} if ($5<=$7) { print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"$7"\t"$8}} ' $outDir/$sample/${sample}_ratio.bed > $outDir/$sample/${sample}_TEDseq_ratio_sort.bed



rm $outDir/$sample/${sample}_ratio.bed

##### FILTERING PUTATTIVE INSERTION ##### 

# remove piles of reads
Rscript $scriptDir/filtre.R $sample $outDir/$sample

# back to bash
sort -k1,1 -k2,2n $outDir/$sample/${sample}_TEDseq_ratio_filtered.bed > $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted.bed

rm $outDir/$sample/${sample}_TEDseq_ratio_filtered.bed

#intersect with targeted_TE_sequences.bed  to remove parental copy
# Insertion with many many reads ends up in two or three insertion, here i merge them keeping only the one containing the bigest number of read 

awk  -v MINRATIO=$min_ratio ' $7>MINRATIO '  $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted.bed >  $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_sincov.bed

while read line 
 do
 echo $line  | awk '{{if ($4=="+") print $1"\t"$2-400"\t"$2+10"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} if ($4=="-")  print $1"\t"$2-10"\t"$2+400"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $outDir/$sample/temp.bed
 awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4}' $outDir/$sample/${sample}_start_site.bed |  bedtools intersect -a $outDir/$sample/temp.bed -b -  -wa -wb | awk '{ sum += $13} END { print sum}' >> $outDir/$sample/temp_cov.txt
done <  $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_sincov.bed

paste $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_sincov.bed $outDir/$sample/temp_cov.txt > $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_cov.bed
rm  $outDir/$sample/temp.bed
rm $outDir/$sample/temp_cov.txt

awk  -v MINCOV=$min_cov  '$3>=MINCOV {print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_cov.bed | bedtools intersect  -v  -a -  -b $ref_TE > $outDir/$sample/${sample}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed


awk 'BEGIN {
    chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9  }
   $1==chromosome  && $2<=interval_start+40 && $4 < depth { next}
    {print chromosome"\t"interval_start"\t"interval_end"\t"depth"\t"signe"\t"nbr_start"\t"nbr_stop"\t"ratio"\t"cov_pic  ; chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9}'  $outDir/$sample/${sample}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed  > $outDir/$sample/${sample}_temp

tail -1 $outDir/$sample/${sample}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed  >> $outDir/$sample/${sample}_temp
sort -k1,1 -k2,2nr  $outDir/$sample/${sample}_temp | tail -1 > $outDir/$sample/${sample}_temp_2nd
sort -k1,1 -k2,2nr  $outDir/$sample/${sample}_temp | awk 'BEGIN {
    chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9}
    $1==chromosome  && $2>=interval_start-40 && $4 < depth { next}
    {print chromosome"\t"interval_start"\t"interval_end"\t"depth"\t"signe"\t"nbr_start"\t"nbr_stop"\t"ratio"\t"cov_pic  ; chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9}'  >> $outDir/$sample/${sample}_temp_2nd
sort -k1,1 -k2.2n $outDir/$sample/${sample}_temp_2nd  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  |  sed '/^\s*$/d' >   $outDir/$sample/${sample}_TEDseq_germline_insertion.bed

            rm $outDir/$sample/${sample}_temp
            rm $outDir/$sample/${sample}_temp_2nd
            rm $outDir/$sample/${sample}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed 
            rm $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_cov.bed
            rm $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted_sincov.bed
            rm $outDir/$sample/${sample}_depth.bed
            rm $outDir/$sample/${sample}_start_site.bed
            rm $outDir/$sample/${sample}_stop_site.bed
            rm $outDir/$sample/${sample}_TEDseq_ratio_filtered_sorted.bed
            rm $outDir/$sample/${sample}_TEDseq_ratio_sort.bed
