# TED-Seq Germline Long-reads pipeline 

TED-Seq Germline Long-reads pipeline is the bioinformatic pipeline used to detect transposable element insertions with ultra high sensitivity using custom made libraries to enrich and capture TE insertions (TED-Seq) generated through a ONT sequencer.

## Table of contents
- [TED-Seq Germline Long-reads pipeline ](#TED-Seq Germline Long-reads pipeline )
  * [Installation](#installation)
    + [Dependencies](#dependencies)
  * [TED-Seq Long-reads main steps](#TED-Seq-long-reads-main-steps)
  * [Input data](#input-data)
  * [Usage](#usage)
	+ [Output](#output)
  * [Cite](#cite)

### Dependencies

TEDSeq Germline Long-reads pipeline requires the following softwares:

* SAMtools (v1.2 or higher) (http://samtools.sourceforge.net/)
* Seqtk (V1.3 or higher) (https://github.com/lh3/seqtk)
* Minimap2 (v2.17 or higher) (https://github.com/lh3/minimap2)
* Cutadapt (v2.6 or higher) (https://cutadapt.readthedocs.io/en/stable/)
* BBmap (v38.18 or higher) (https://sourceforge.net/projects/bbmap/)
* dplyr (Version 1.1.3) (https://dplyr.tidyverse.org/)
* bedtools (v2.20.1 or higher) (https://bedtools.readthedocs.io/en/latest/content/installation.html)


## TED-Seq long-reads main steps


1- Extraction of reads containing the desired TE family and barcode to analyze and trimming of the sequence, removal of PCR duplicates.

2- Mapping of the remaining reads to the reference genome.

3- Detection of insertion side.

4- Discarding insertions supported by few amplicons (possibly derived from somatic transposition).


## Input data

This pipeline requires as an input:

1. A raw fastq file derived from TED-seq data demultiplexed for the corresponding barcodes
2. A tab-delimited file containing TE family coordinates across the genome. The file should have the following format: 

	- Chromosome name (should be the same nomenclature as in the reference genome fasta file)
	- TE start
	- TE end
3. Indexed reference genome (bowtie2-build)
4. Nucleotide extremity of the reference copy to be analyzed

## Usage
Prepare fasta files with the donor TE flanking sequences and a bed file with the reference positions for this given TE.

Prepare index files for reference genome
Detect TE insertions
To launch the pipeline use the following script:

```
bash Long_read_TEDSeq_TE_germline.sh 
--sample        Name of the sample (required).
--reads         Full path to the raw ONT fastq file (required).
--refDir        Path to reference genome (required).
--outDir        Output directory (required).
--te_family     Desired family analyzed (default: ATCOPIA93).
--barcode       Desired P7 barcode used, indicated as a number from 1 to 48 (required).
--cores         Number of cores to use (default: 20).
--min_cov       Minimum coverage for TE detection.
--min_ratio     Minimum ratio (default: 0.8).
--ref_TE        Reference TE bed file (required).
--scriptDir     Directory containing the scripts and the flanking sequences (required).

```

### Output

TED-Seq Germline Long-reads pipeline will return as an output:

* Two bam files of the mapping of the deduplicated and trimmed long-reads to the reference genome
  $sample_clip_disc-local.bam
 
* A bed file containing the non-reference germline insertions
  $sample_TEDseq_germline_insertion_corrected.bed


The content of `$sample_TEDseq_germline_insertion_corrected.bed` contains the following fields:

```
Chromosome
Insertion position Start
Insertion position End
Coverage of the insertion
TE insertion orientation
Number of reads at the 5' extremity
Number of reads at the 3' extremity
Ratio of reads 5'/ 3'
Maximum number of reads at the interval supporting the insertion

```

## Cite

If you use this software, please cite:

Pol Vendrell-Mir, Basile Leduque, and Leandro Quadrana. “Ultra-Sensitive Detection of Transposon Insertions Across Multiple Families by Transposable Element Display Sequencing,” August 22, 2024. https://doi.org/10.1101/2024.08.21.608910.
