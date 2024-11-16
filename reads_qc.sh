#!/bin/bash
# Quality control of long and short reads before de novo genome assembly.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 29 0ctober 2024; latest update: 16 November 2024

# Default parameters ###############
sample='sample'
outdir="$PWD"

long_reads_eval_raw=false
long_reads_trim=false
long_reads_fastp_mean_qual=10
long_reads_fastp_window_size=4
long_reads_min_len=1000
long_reads_min_qual=10

short_reads_eval_raw=false
short_reads_trim=false
short_reads_fastp_mean_qual=30
short_reads_fastp_window_size=10
short_reads_min_len=50
short_reads_min_qual=25

rm_fastqc_zip=false
threads=1
memory=512

# Help information ###############
display_usage() {
    echo "
    This script processes long and/or paired-end short reads for good quality.
    Dependencies: NanoPlot, fastplong, fastp, nanoq, fastqc, seqkit. Please ensure all dependencies are
    accessible through the environmental variable PATH.

    Parameters:
      --r=*: fastq or fastq.gz file of Nanopore long reads
      --r1=*: fastq or fastq.gz file of paired-end short reads, forward
      --r2=*: fastq or fastq.gz file of paired-end short reads, reverse

      --outdir=*: output directory, no forward slash (default: ${outdir})
      --sample=*: sample name (default: ${sample})

      --long_reads_eval_raw: switch on the quality evaluation of raw long reads (default: off)
      --long_reads_trim: switch on processing long reads and evaluate the quality of processed reads (default: off)
      --long_reads_fastp_mean_qual=*: per-window average quality in long reads for fastp (default: ${long_reads_fastp_mean_qual})
      --long_reads_fastp_window_size=*: size of fastp sliding windows for cutting the heads and tails of long reads (default: ${long_reads_fastp_window_size})
      --long_reads_min_len=*: minimum length of long reads (default: ${long_reads_min_len})
      --long_reads_min_qual=*: minimum average quality of long reads (default: ${long_reads_min_qual})

      --short_reads_eval_raw: switch on the quality evaluation of raw short reads (default: off)
      --short_reads_trim: switch on processing short reads and evaluate the quality of processed reads (default: off)
      --short_reads_fastp_mean_qual=*: per-window average quality in short reads for fastp (default: ${short_reads_fastp_mean_qual})
      --short_reads_fastp_window_size=*: size of fastp sliding windows for cutting the heads and tails of short reads (default: ${short_reads_fastp_window_size})
      --short_reads_min_len=*: minimum length of short reads (default: ${short_reads_min_len})
      --short_reads_min_qual=*: minimum average quality of short reads (default: ${short_reads_min_qual})

      --threads=*: number of threads (default: ${threads})
      --memory=*: memory size in MB (default: ${memory})
      --rm_fastqc_zip: remove Zip files from FastQC's outputs (default: off)
    
    Example useage:
      Assembly_toolkit/reads_qc.sh --r=... --r1=... [other paramters]

    Output directory: 'sample' in the current working directory.
    "
}

if [ -z $1 ]; then
    display_usage
    exit
fi

# Other functions ###############
make_dir() {
    if [ ! -d "$1" ]; then
        echo "[$(date)] Create directory $1"
        mkdir -p "$1"
    fi
}

rm_file() {
    if [ -f "$1" ]; then
        rm "$1"
    else
        echo "Warning: $1 was not found."
    fi
}

# Read arguments ###############
for i in "$@"; do
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
        --short_reads_eval_raw)
        short_reads_eval_raw=true
        ;;
        --short_reads_trim)
        short_reads_trim=true
        ;;
        --short_reads_fastp_mean_qual=*)
        short_reads_fastp_mean_qual="${i#*=}"
        ;;
        --short_reads_fastp_window_size=*)
        short_reads_fastp_window_size="${i#*=}"
        ;;
        --short_reads_min_len=*)
        short_reads_min_len="${i#*=}"
        ;;
        --short_reads_min_qual=*)
        short_reads_min_qual="${i#*=}"
        ;;
        --threads=*)
        threads="${i#*=}"  # Number of threads
        ;;
        --memory=*)
        memory="${i#*=}"  # Memory size in MB
        ;;
        --rm_fastqc_zip)
        rm_fastqc_zip=true
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

# Main ###############
make_dir "$outdir"

# Long-read QC ###############
if [ -f "$long_reads" ]; then
    # Evaluate the quality of raw reads
    if $long_reads_eval_raw; then  # Using the arithmetic evaluation for the Boolean variable (https://kodekloud.com/blog/declare-bash-boolean-variable-in-shell-script/)
        echo "[$(date)] Inspect the quality of raw long reads in $long_reads"
        outdir_raw="${outdir}/nanopore/1_raw/quality"
        make_dir "$outdir_raw"
        NanoPlot --fastq $long_reads --outdir "${outdir_raw}/nanoplot" --threads $threads --prefix $sample --title "Unprocessed MinION reads of $sample" --minlength 1 --drop_outliers --readtype 1D --plots kde
        fastqc --outdir $outdir_raw --noextract --format fastq --threads $threads --memory $memory $long_reads
        seqkit stats --all --threads $threads --tabular --basename --seq-type dna $long_reads > "${outdir_raw}/${sample}_seqkit_stats_nanopore_raw.tsv"
        if $rm_fastqc_zip; then
            rm_file "${outdir_raw}/${sample}_fastqc.zip"
        fi
    else
        echo "[$(date)] Skip quality assessment of long reads."
    fi

    # Quality process of reads and quality report
    if $long_reads_trim; then
        echo "[$(date)] Process long reads from $long_reads for high quality"
        outdir_processed="${outdir}/nanopore/2_processed"
        make_dir "${outdir_processed}/quality"

        # Read trimming and filtering
        output_basename="${sample}_hQ${long_reads_fastp_mean_qual}tQ${long_reads_fastp_mean_qual}w${long_reads_fastp_window_size}L${long_reads_min_len}"
        fastp_output="${outdir_processed}/${output_basename}.fastq"
        fastplong -i $long_reads -o $fastp_output --html "${outdir_processed}/quality/${output_basename}_fastplong.html" --json /dev/null \
        --thread $threads --cut_front --cut_tail --cut_window_size $long_reads_fastp_window_size --cut_mean_quality $long_reads_fastp_mean_qual --length_required $long_reads_min_len
        output_basename+="_aQ$long_reads_min_qual"
        nanoq_output="${outdir_processed}/${output_basename}.fastq"
        nanoq --input $fastp_output --output-type u --output $nanoq_output --min-qual $long_reads_min_qual --report "${outdir_processed}/quality/${output_basename}_nanoq_summary.txt"
        
        # Evaluate the quality of processed reads
        echo "[$(date)] Inspect the quality of processed long reads from $nanoq_output"
        NanoPlot --fastq $nanoq_output --outdir "${outdir_processed}/quality/nanoplot" --threads $threads --prefix $output_basename --title "Processed MinION reads of $sample (${quality_suffix})" --minlength 1 --drop_outliers --readtype 1D --plots kde
        fastqc --outdir "${outdir_processed}/quality" --noextract --format fastq --threads $threads --memory $memory $nanoq_output
        seqkit stats --all --threads $threads --tabular --basename --seq-type dna $nanoq_output > "${outdir_processed}/quality/${output_basename}_seqkit_stats_nanopore.tsv"
        if $rm_fastqc_zip; then
            rm_file "${outdir_processed}/quality/${output_basename}_fastqc.zip"
        fi

        # Compress the output file
        echo "[$(date)] Compress $nanoq_output and remove the intermediate file $fastp_output"
        gzip $nanoq_output
        rm $fastp_output
        echo "[$(date)] Output read file: ${nanoq_output}.gz"
    else
        echo "[$(date)] Skip quality process of long reads."
    fi
fi

# Short-read QC ###############
if [ -f "$r1" ] && [ -f "$r2" ]; then
    # Assessment of raw-read quality
    if $short_reads_eval_raw; then
        echo "[$(date)] Evaluate the quality of raw short reads in $r1 and $r2"
        outdir_raw="${outdir}/illumina/1_raw/quality"
        make_dir "$outdir_raw"
        fastqc --outdir "$outdir_raw" --noextract --nogroup --format fastq --threads "$threads" "$r1" "$r2"
        seqkit stats --all --threads $threads --tabular --basename --seq-type dna "$r1" "$r2" > "${outdir_raw}/${sample}_seqkit_stats_illumina_raw.tsv"
        if $rm_fastqc_zip; then
            rm_file "${outdir_raw}/${sample}_1_fastqc.zip"
            rm_file "${outdir_raw}/${sample}_2_fastqc.zip"
        fi
    else
        echo "[$(date)] Skip quality assessment of raw short reads."
    fi

    if $short_reads_trim; then
        echo "[$(date)] Process paired-end short reads in $r1 and $r2 for high quality"
        outdir_processed="${outdir}/illumina/2_processed"
        make_dir "${outdir_processed}/quality"

        # Read trimming and filtering
        output_basename="${sample}_hQ${short_reads_fastp_mean_qual}tQ${short_reads_fastp_mean_qual}w${short_reads_fastp_window_size}aQ${short_reads_min_qual}L$short_reads_min_len"
        fastp_output1="${outdir_processed}/${output_basename}_1.fastq"
        fastp_output2="${outdir_processed}/${output_basename}_2.fastq"
        fastp --in1 $r1 --in2 $r2 --out1 $fastp_output1 --out2 $fastp_output2 --html "${outdir_processed}/quality/${output_basename}_fastp.html" \
        --json /dev/null --thread $threads --cut_front --cut_tail --cut_window_size $short_reads_fastp_window_size --cut_mean_quality $short_reads_fastp_mean_qual \
        --average_qual $short_reads_min_qual --length_required $short_reads_min_len
        fastqc --outdir "${outdir_processed}/quality" --noextract --format fastq --threads $threads --memory $memory $fastp_output1 $fastp_output2
        seqkit stats --all --threads $threads --tabular --basename --seq-type dna $fastp_output1 $fastp_output2 > "${outdir_processed}/quality/${output_basename}_seqkit_stats_illumina.tsv"
        if $rm_fastqc_zip; then
            rm_file "${outdir_processed}/quality/${output_basename}_1_fastqc.zip"
            rm_file "${outdir_processed}/quality/${output_basename}_2_fastqc.zip"
        fi

        # Compress the output file
        echo "[$(date)] Compress $fastp_output1 and $fastp_output2"
        gzip $fastp_output1
        gzip $fastp_output2
        echo "[$(date)] Output read file: ${fastp_output1}.gz and ${fastp_output2}.gz"
    fi
fi
