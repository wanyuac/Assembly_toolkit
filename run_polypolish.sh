#!/bin/bash
# This script has been adapted from run_polypolish_legacy.sh for Polypolish v0.6.0 and onwards.
# It is not compatible with Polypolish v0.5.0, for which run_polypolish_legacy.sh should be used.
# Update note of v0.6.0: https://github.com/rrwick/Polypolish/releases/tag/v0.6.0.
# Reference: https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish
# Copyright (C) 2022-2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 1 May 2022; latest update: 17 November 2024

# Set default values ###############
SCRIPT_VERSION=2.0
i='isolate'
n="$i"
outdir=polypolish_output
t=1

# Function definitions ###############
display_parameters() {
    echo "
    run_polypolish.sh v$SCRIPT_VERSION

    This script polishes an input genome assembly with Polypolish and paired-end short reads.

    Prerequisites: please ensure bwa aligner is accessible in PATH.

    Parameters (six in total):
      -a=*: path and filename of the input assembly in FASTA format (mandatory)
      -r=*: path to the directory of paired-end short reads (without the end forward
            slash), where read files' names follow format [isolate name]_[1,2].fastq.gz (mandatory)
      -i=*: isolate name (default: ${i})
      -s=*: suffix of input FASTQ filename: [isolate name]_[suffix]_[1,2].fastq.gz (default: none)
      -n=*: prefix of output filenames (default: ${n})
      -o=*: path to output directory (default: ${outdir})
      -t=*: number of threads (default: ${t})
    
    Example command:
      ~/bin/Assembly_toolkit/run_polypolish.sh -a=1_isolate1_medaka.fasta -r=reads/illumina -i=isolate1 -n="2_isolate1" -o="\$PWD" -t=8
    
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
        -s=*)
        s="${arg#*=}"
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
    if [ -z "$s" ]
    then
        r1="${read_dir}/${i}_1.fastq.gz"
        r2="${read_dir}/${i}_2.fastq.gz"
    else
        r1="${read_dir}/${i}_${s}_1.fastq.gz"
        r2="${read_dir}/${i}_${s}_2.fastq.gz"
    fi
    if [ ! -e "$r1" ] || [ ! -e "$r2" ]  # $r1 and $r2 can be files or symbolic links
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
        echo "Create output directory $outdir" >> "${output_prefix}.txt"
        mkdir "$outdir"
    fi
    echo "[$(date)] Start to polish $fasta_in of isolate $i with reads from $read_dir and save outputs in ${outdir}/" >> "${output_prefix}.txt"
    echo "Read files: $r1 $r2" >> "${output_prefix}.txt"

    # Filtering read alignments to exclude those of unusually large insert sizes ###############
    echo "[$(date)] Creating SAM files for isolate $i" >> "${output_prefix}.txt"
    tm="sams_$i"
    mkdir $tm
    bwa index $fasta_in
    bwa mem -a -t $t $fasta_in $r1 > ${tm}/unfiltered_1.sam
    bwa mem -a -t $t $fasta_in $r2 > ${tm}/unfiltered_2.sam
    polypolish filter --in1 ${tm}/unfiltered_1.sam --in2 ${tm}/unfiltered_2.sam --out1 ${tm}/filtered_1.sam --out2 ${tm}/filtered_2.sam 1>> "${output_prefix}.txt" 2>&1

    # Polish the input assembly ###############
    echo "[$(date)] Polishing assembly $fasta_in with Polypolish" >> "${output_prefix}.txt"
    polypolish polish $fasta_in ${tm}/filtered_1.sam ${tm}/filtered_2.sam 1>"${output_prefix}.fna" 2>>"${output_prefix}.txt"
    echo "[$(date)] Finished polishing $fasta_in" >> "${output_prefix}.txt"
    rm -rf $tm

    # Remove non-character content from log files ###############
    perl -p -e 's/\x1b\[[0-9;]*[mG]//g' ${output_prefix}.txt > ${output_prefix}.log
    rm ${output_prefix}.txt

    # Clean up ###############
    rm "$fasta_in".*
    echo "This Polypolish job successfully finished. Results have been saved in ${outdir}." >> ${output_prefix}.log
fi
