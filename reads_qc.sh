#!/bin/bash
# Quality control of long and short reads before de novo genome assembly.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 29 0ctober 2024; latest update: 30 0ctober 2024

# Default parameters ###############
sample='sample'
outdir="$PWD"
threads=1
memory=512
long_reads_fastp_mean_qual=10
long_reads_fastp_window_size=4
long_reads_min_len=1000
long_reads_min_qual=10
long_reads_eval_raw=false
long_reads_trim=false

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
      --outdir=*: output directory, no forward slash (default: ${outdir})
      --eval_raw: which tasks to perform: raw (only evaluate the quality of raw reads); process
      --long_reads_fastp_mean_qual=*: per-window average quality in long reads for fastp (default: ${long_reads_fastp_mean_qual})
      --long_reads_fastp_window_size=*: size of fastp sliding windows for cutting the heads and tails of long reads (default: ${long_reads_fastp_window_size})
      --long_reads_min_len=*: minimum length of long reads (default: ${long_reads_min_len})
      --long_reads_min_qual=*: minimum average quality of long reads (default: ${long_reads_min_qual})
      --long_reads_eval_raw: switch on evaluation of raw-read quality
      --long_reads_trim: switch on processing long reads and evaluate the quality of processed reads (default: off)
      --threads=*: number of threads (default: ${threads})
      --memory=*: Memory size in MB (default: ${memory})
    
    Example useage:
      Assembly_toolkit/reads_qc.sh --r=... --r1=... [other paramters]

    Output directory: 'sample' in the current working directory.
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
        --long_reads_fastp_mean_qual=*)
        long_reads_fastp_mean_qual="${i#*=}"
        ;;
        --long_reads_fastp_window_size=*)
        long_reads_fastp_window_size="${i#*=}"
        ;;
        --long_reads_min_len=*)
        long_reads_min_len="${i#*=}"
        ;;
        --long_reads_min_qual=*)
        long_reads_min_qual="${i#*=}"
        ;;
        --long_reads_eval_raw)
        long_reads_eval_raw=true
        ;;
        --long_reads_trim)
        long_reads_trim=true
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

# Main ###############
if [ ! -d "$outdir" ]; then
    echo "[$(date)] Create the output directory $outdir"
    mkdir -p $outdir
fi

# Long-read QC ###############
if [ -f "$long_reads" ]; then
    # Evaluate the quality of raw reads
    if $long_reads_eval_raw; then  # Using the arithmetic evaluation for the Boolean variable (https://kodekloud.com/blog/declare-bash-boolean-variable-in-shell-script/)
        echo "[$(date)] Inspect the quality of long reads from $long_reads"
        outdir_raw="${outdir}/nanopore/raw/quality"
        mkdir -p "$outdir_raw"
        NanoPlot --fastq $long_reads --outdir "${outdir_raw}/nanoplot" --threads $threads --prefix $sample --title "Unprocessed MinION reads of $sample" --minlength 1 --drop_outliers --readtype 1D --plots kde
        fastqc --outdir $outdir_raw --noextract --format fastq --threads $threads --memory $memory $long_reads
        seqkit stats --all --threads $threads --tabular --basename --seq-type dna $long_reads > "${outdir_raw}/${sample}_seqkit_summary_nanopore_raw.tsv"
    fi

    # Quality process of reads
    if $long_reads_trim; then
        echo "[$(date)] Process long reads from $long_reads for high quality"
        outdir_processed="${outdir}/nanopore/processed"
        mkdir -p "${outdir_processed}/quality"
        quality_suffix="hQ${long_reads_fastp_mean_qual}tQ${long_reads_fastp_mean_qual}w${long_reads_fastp_window_size}L1000"
        fastp_output="${outdir_processed}/${sample}_${quality_suffix}.fastq"
        fastp -i $long_reads -o $fastp_output --html "${outdir_processed}/quality/${sample}_${quality_suffix}_fastp.html" --json /dev/null --thread $threads --cut_front --cut_tail --cut_window_size $long_reads_fastp_window_size --cut_mean_quality $long_reads_fastp_mean_qual --length_required $long_reads_min_len
        quality_suffix="${quality_suffix}_aQ${long_reads_min_qual}"
        nanoq_output="${outdir_processed}/${sample}_${quality_suffix}.fastq"
        nanoq --input $fastp_output --output-type u --output $nanoq_output --min-qual $long_reads_min_qual --report "${outdir_processed}/quality/${sample}_${quality_suffix}_nanoq_summary.txt"
        
        # Evaluate the quality of processed reads
        echo "[$(date)] Inspect the quality of processed long reads from $nanoq_output"
        NanoPlot --fastq $nanoq_output --outdir "${outdir_processed}/quality/nanoplot" --threads $threads --prefix $sample --title "Processed MinION reads of $sample (${quality_suffix})" --minlength 1 --drop_outliers --readtype 1D --plots kde
        fastqc --outdir "${outdir_processed}/quality" --noextract --format fastq --threads $threads --memory $memory $nanoq_output
        seqkit stats --all --threads $threads --tabular --basename --seq-type dna $nanoq_output > "${outdir_processed}/quality/${sample}_seqkit_summary_nanopore_${quality_suffix}.tsv"

        # Compress the output file
        echo "[$(date)] Compress $nanoq_output and remove the intermediate file $fastp_output"
        gzip "$nanoq_output"
        rm "$fastp_output"
        echo "[$(date)] Output read file: ${nanoq_output}.gz"
    fi
fi

# Short-read QC ###############
#if [ -f "$r1" ] && [ -f "$r2" ]; then
#    echo "[$(date)] This functionality will be developed."
#else
#    echo "[$(date)] Error: $r1 and/or $r2 were not accessible."
#fi