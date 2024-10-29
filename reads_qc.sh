#!/bin/bash
# Quality control of long and short reads before de novo genome assembly.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 29 0ctober 2024; latest update: 29 0ctober 2024

# Default parameters ###############
sample='sample'
outdir="$PWD"
threads=1
memory=512

# Help information ###############
display_usage() {
    echo "
    This script processes long and/or paired-end short reads for good quality.
    Dependencies: NanoPlot, fastp, nanoq, fastqc, seqkit. Please ensure all dependencies are
    accessible through the environmental variable PATH.

    Parameters:
      --r=*: fastq or fastq.gz file of Nanopore long reads
      --r1=*: fastq or fastq.gz file of paired-end short reads, forward
      --r2=*: fastq or fastq.gz file of paired-end short reads, reverse
      --sample=*: sample name (default: ${sample})
      --outdir=*: output directory, no forward slash (default: ${$outdir})
      --threads=*: number of threads (default: ${threads})
      --memory=*: Memory size in MB (default: ${memory})
    
    Example useage:
      Assembly_toolkit/reads_qc.sh --r=... --r1=... [other paramters]

    Output directory: 'sample_name' in the current working directory.
    "
}

if [ -z $1 ]
then
    display_usage
    exit
fi

# Read arguments ###############
for i in "$@"
do
    case $i in
        --r=*)
        long_reads="${i#*=}"  # fastq or fastq.gz file of Nanopore long reads
        ;;
        --r1=*)
        r1="${i#*=}"  # fastq or fastq.gz file of paired-end short reads, forward
        ;;
        --r2=*)
        r2="${i#*=}"  # fastq or fastq.gz file of paired-end short reads, reverse
        ;;
        --sample=*)
        sample="${i#*=}"
        ;;
        --outdir=*)
        outdir="${i#*=}"  # Output directory
        ;;
        --threads=*)
        threads="${i#*=}"  # Number of threads
        ;;
        --memory=*)
        memory="${i#*=}"  # Memory size in MB
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

if [ ! -f "$r" ] && [ ! -f "$r1" ] && [ ! -f "$r2" ]; then
    echo "[$(date)] Error: No input read file was found."

# Main ###############
if [ ! -d "$outdir" ]; then
    echo "[$(date)] Create the output directory $outdir"
    mkdir -p $outdir
fi

# Long-read QC ###############
if [ -f "$long_reads" ]; then
    echo "[$(date)] Quality control of long reads from $long_reads"

    # Raw reads
    outdir_raw="${outdir}/nanopore/raw/quality"
    mkdir -p "$outdir_raw"
    NanoPlot --fastq $long_reads --outdir "${outdir_raw}/nanoplot" --threads $threads --prefix $sample --title "Unprocessed MinION reads of $sample" --minlength 1 --drop_outliers --readtype 1D --plots kde
    fastqc --outdir $outdir_raw --noextract --format fastq --threads $threads --memory $memory $long_reads
    seqkit stats --all --threads $threads --tabular --basename --seq-type dna $long_reads > "${outdir_raw}/${sample}_seqkit_summary_nanopore_raw.tsv"

    # Processed reads
    outdir_proc="${outdir}/nanopore/processed"
    mkdir -p "${outdir_processed}/quality"
    # Process reads
    NanoPlot --fastq $long_reads --outdir "${outdir_processed}/quality/nanoplot" --threads $threads --prefix $sample --title "Processed MinION reads of $sample" --minlength 1 --drop_outliers --readtype 1D --plots kde
    fastqc --outdir "${outdir_processed}/quality" --noextract --format fastq --threads $threads --memory $memory $long_reads
    seqkit stats --all --threads $threads --tabular --basename --seq-type dna $long_reads > "${outdir_processed}/quality/${sample}_seqkit_summary_nanopore_processed.tsv"
fi
