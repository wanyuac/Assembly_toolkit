#!/bin/bash
# run_polypolish.sh [isolate name] [input FASTA file] [directory of short reads] [output directory] [prefix of output files]
# The prefix of output filenames can be "1_${sample_name}" or "3_${sample_name}", etc.
# Prerequisites: python v3, BWA
# Reference: https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 1 May 2022; latest update: 12 Nov 2023

SCRIPT_VERSION=1.0.0

# Function definitions ###############
# Run './run_polypolish.sh' to display the information below.
display_parameters() {
    echo "
    run_polypolish.sh v$SCRIPT_VERSION

    This script polishes an input genome assembly with Polypolish and paired-end short reads.

    Prerequisites: please ensure bwa aligner is accessible in PATH.

    Parameters (six in total):
      -a=*: path and filename of the input assembly in FASTA format (mandatory)
      -r=*: path to the directory of paired-end short reads (without the end forward
            slash), where read files' names follow format [sample name]_[1,2].fastq.gz (mandatory)
      -i=*: isolate name (default: isolate)
      -n=*: prefix of output filenames (default: isolate)
      -o=*: path to output directory (default: polypolish_output)
      -t=*: number of threads (default: 2)
    
    Example command:
      /usr/local/bin/Assembly_toolkit/run_polypolish.sh -a=assembly.fasta -r=reads/illumina \\
        -i=isolate_1 -o=polca -t=8
    
    Outputs:
      1. a polished assembly [o]/[n]_polypolish.fna in FASTA format
      2. a log of messages from Polypolish's scripts [o]/[n]_polypolish.log
    "
}

MYDIR="$(dirname "$(readlink -f "$0")")"  # Directory of this script
source $MYDIR/modules.sh  # Function print_failure_message

# Print parameter information ###############
if [ -z "$1" ]
then
    display_parameters
    exit
fi

# Set default values ###############
i=isolate
n="$i"
outdir=polypolish_output
t=2

# Read parameters ###############
for arg in "$@"
do
    case $arg in
        -a=*)
        fasta_in="${arg#*=}"
        ;;
        -r=*)
        read_dir="${arg#*=}"
        ;;
        -i=*)
        i="${arg#*=}"
        ;;
        -n=*)
        n="${arg#*=}"
        ;;
        -o=*)
        outdir="${arg#*=}"
        ;;
        -t=*)
        t="${arg#*=}"
        ;;
        *)
        ;;
    esac
done

# Run Polypolish for the input genome assembly ###############
if [ -z "$(which bwa)" ] 
then
    echo "Error: bwa was not accessible"
    print_failure_message "$i"
    exit
fi

if [ -d "$read_dir" ]
then
    r1="${read_dir}/${i}_1.fastq.gz"
    r2="${read_dir}/${i}_2.fastq.gz"
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]
    then
        echo "Error: read file $r1 and/or $r2 were not found."
        print_failure_message "$i"
        exit
    fi
else
    echo "Error: read directory $read_dir was not found."
    print_failure_message "$i"
fi

if [ -f "$fasta_in" ]
then
    output_prefix="${outdir}/${n}_polypolish"
    echo "run_polypolish.sh v$SCRIPT_VERSION" > "${output_prefix}.txt"

    # Set up output directories and filenames
    if [ ! -d "$outdir" ]
    then
        echo "Create output directory $outdir"
        mkdir "$outdir"
    fi
    echo "Start to polish $fasta_in of isolate $i with reads from $read_dir and save outputs in ${outdir}/"

    # Filtering read alignments to exclude those of unusually large insert sizes ###############
    echo "$(date): Creating SAM files for isolate $i" >> "${output_prefix}.txt"
    tm="sams_$i"
    mkdir $tm
    bwa index $fasta_in
    bwa mem -a -t $t $fasta_in $r1 > ${tm}/unfiltered_1.sam
    bwa mem -a -t $t $fasta_in $r2 > ${tm}/unfiltered_2.sam
    polypolish_insert_filter.py --in1 ${tm}/unfiltered_1.sam --in2 ${tm}/unfiltered_2.sam --out1 ${tm}/filtered_1.sam --out2 ${tm}/filtered_2.sam 1>> "${output_prefix}.txt" 2>&1

    # Polish the input assembly ###############
    echo "$(date): Polishing assembly $fasta_in with Polypolish" >> "${output_prefix}.txt"
    polypolish $fasta_in ${tm}/filtered_1.sam ${tm}/filtered_2.sam 1>"${output_prefix}.fna" 2>>"${output_prefix}.txt"
    echo "$(date): Finished polishing $fasta_in" >> "${output_prefix}.txt"
    rm -rf $tm

    # Remove '_polypolish' from sequence headers to preserve the original ones ###############
    echo "$(date): processing the polished assembly and log files" >>"${output_prefix}.txt"
    sed -i 's/_polypolish//g' ${output_prefix}.fna

    # Remove non-character content from log files ###############
    perl -p -e 's/\x1b\[[0-9;]*[mG]//g' ${output_prefix}.txt > ${output_prefix}.log
    rm ${output_prefix}.txt

    # Clean up ###############
    rm "$fasta_in".*
    echo "This Polypolish job successfully finished. Results have been saved in ${outdir}."
fi
