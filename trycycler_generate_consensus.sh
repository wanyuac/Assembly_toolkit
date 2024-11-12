#!/bin/bash
# Run the last three steps of Trycycler: msa, partition, and consensus.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 12 November 2024; last update: 12 November 2024

# Help information ###############
display_usage() {
    echo "
    This script carries out the last three steps of Trycycler: msa, partition, and consensus.

    Parameters:
      -r=*: The input long reads
      -d=*: The input directory or the output directory of 'trcycyler cluster', which serves as the input directory for 'trycycler reconcile' (default: 2_clusters)
      -t=*: Number of threads (default: 1)
    
    Example useage:
      ~/bin/Assembly_toolkit/trycycler_generate_consensus.sh -r=sample1.fastq.gz -d=2_clusters -t=16 >trycycler_generate_consensus.log 2>&1

    Dependencies: trycycler, perl.
    "
}

if [ -z "$1" ]
then
    display_usage
    exit
fi

# Helper function ###############
clean_log() {
    perl -pe 's/\x1b\[[0-9;]*[mG]//g' "$1"
}

# Main procedure ###############
# Set default values ===============
dir_in='2_clusters'
t=1

# Read parameters ===============
for arg in "$@"
do
    case $arg in
        -r=*)
        long_reads="${arg#*=}"
        ;;
        -d=*)
        dir_in="${arg#*=}"
        ;;
        -t=*)
        t="${arg#*=}"
        ;;
        *)
        ;;
    esac
done

# Sanity check ===============
if [ -f "$long_reads" ]
then
    echo "[$(date)] Input read file: $long_reads"
else
    echo "[$(date)] Error: read file $long_reads was not found."
    exit 1
fi

if [ $t -lt 1 ]
then
    echo "[$(date)] Error: the number of threads is negative. Reset to 1."
    t=1
fi
echo "[$(date)] Number of threads: $t"

# Multi-sequence aligment ===============
cluster_dirs=( $(ls -1 -d "$dir_in"/cluster_*/) )
echo "[$(date)] Found ${#cluster_dirs[@]} clusters in $dir_in"
echo "[$(date)] Conduct multi-sequence alignment"

for d in "${cluster_dirs[@]}"
do
    trycycler msa --threads "$t" --cluster_dir "${d%/}"
done

# Partitioning reads
echo "[$(date)] Partition long reads"

trycycler partition --threads "$t" --reads "$long_reads" --cluster_dirs "$dir_in"/cluster_00?

# Generating a consensus ===============
echo "[$(date)] Produce consensus sequences"

for d in "${cluster_dirs[@]}"
do
    trycycler consensus --threads "$t" --cluster_dir "${d%/}" >> trycycler_consensus.txt 2>&1
done

# Finish up ===============
clean_log trycycler_consensus.txt > trycycler_consensus.log && rm trycycler_consensus.txt
