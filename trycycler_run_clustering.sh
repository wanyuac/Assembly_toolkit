#!/bin/bash
# This script is a wrapper for command 'trycycler cluster'.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 13 November 2024; last update: 13 November 2024

# Help information ###############
display_usage() {
    echo "
    This script is a wrapper for command \'trycycler cluster\'.

    Parameters:
      -i=*: Input directory - the output directory containing input FASTA files (*.fasta)
      -r=*: Input file of long reads
      -o=*: Output directory (default: 2_clusters)
      -s=*: Name suffix for the output directory and log file (default: none)
      -t=*: Number of threads (default: 1)
    
    Example useage:
      ~/bin/Assembly_toolkit/trycycler_exclude_bad_contigs.sh -i=1_assemblies_min25x -r=reads/sample1.fastq.gz -s=min25x_r1 -t=8

    Dependency: trycycler
    "
}

# Main ###############
if [ -z "$1" ]
then
    display_usage
    exit
fi

# Set default values
outdir=2_clusters
t=1

for arg in "$@"
do
    case $arg in
        -i=*)
        input_dir="${arg#*=}"
        ;;
        -r=*)
        long_reads="${arg#*=}"
        ;;
        -o=*)
        outdir="${arg#*=}"
        ;;
        -s=*)
        suffix="${arg#*=}"
        ;;
        -t=*)
        threads="${arg#*=}"
        ;;
        *)
        ;;
    esac
done

if [ $threads -lt 1 ]
then
    echo "Parameter error: the number of threads is negative. Reset it to 1."
    threads=1
fi

if [ ! -z "$suffix" ]
then
    output_basename="${outdir}_${suffix}"
else
    output_basename="$outdir"
fi

if [ -f "$long_reads" ]
then
    trycycler cluster --assemblies "$input_dir"/*.fasta --reads "$long_reads" --threads "$threads" --out_dir "$output_basename" > "${output_basename}.tmp" 2>&1
    perl -pe 's/\x1b\[[0-9;]*[mG]//g' "${output_basename}.tmp" > "$output_basename/${output_basename}.log" && rm "${output_basename}.tmp"
else
    echo "Error: file $long_reads was not found."
fi
