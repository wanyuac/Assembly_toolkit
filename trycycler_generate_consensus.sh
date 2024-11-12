#!/bin/bash
# Run the last three steps of Trycycler: msa, partition, and consensus.
#
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 12 November 2024; last update: 12 November 2024

# Help information ###############
display_usage() {
    echo "
    This script renames bad contigs to exclude them from inputs of command 'trycycler reconcile'.

    Parameters:
      -r=*: The input long reads
      -d=*: The input directory or the output directory of 'trcycyler cluster', which serves as the input directory for 'trycycler reconcile' (default: 2_clusters)
      -s=*: The common suffix to add to log files (default: none)
      -t=*: Number of threads (default: 1)
    
    Example useage:
      ~/bin/Assembly_toolkit/trycycler_generate_consensus.sh -i=2_clusters -n=3 -s=round1 -t=16 -p

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
# Set default values
dir_in='2_clusters'
t=1
draw_plot=false

# Read parameters
for arg in "$@"
do
    case $arg in
        -r=*)
        reads="${arg#*=}"
        ;;
        -d=*)
        dir_in="${arg#*=}"
        ;;
        -s=*)
        s="${arg#*=}"
        ;;
        -t=*)
        t="${arg#*=}"
        if [ $t -lt 1 ]
        then
            echo "Error: the number of threads is negative. Reset to 1."
            t=1
        fi
        ;;
        *)
        ;;
    esac
done

# Multi-sequence aligment
for i in `seq 1 3`; do trycycler msa --threads 16 --cluster_dir "2_clusters_min${min_depth}x_round2/cluster_00${i}"; done

# Partitioning reads
trycycler partition --threads 16 --reads $long_reads --cluster_dirs "2_clusters_min${min_depth}x_round2"/cluster_00?

# Generating a consensus
for i in `seq 1 3`
do
    trycycler consensus --threads 16 --cluster_dir "2_clusters_min${min_depth}x_round2/cluster_00$i" >> trycycler_consensus.txt 2>&1
done

clean_log trycycler_consensus.txt > trycycler_consensus.log && rm trycycler_consensus.txt
