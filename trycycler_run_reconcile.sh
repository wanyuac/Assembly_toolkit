#!/bin/bash
# This is a wrapper for the 'trycycler reconcile' command.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 6 November 2024; last update: 7 November 2024

# Help information ###############
display_usage() {
    echo "
    This script renames bad contigs to exclude them from inputs of command 'trycycler reconcile'.

    Parameters:
      -r=*: The input long reads
      -d=*: The input directory or the output directory of 'trcycyler cluster', which serves as the input directory for 'trycycler reconcile' (default: 2_clusters)
      -s=*: The common suffix to add to log files (default: none)
      -t=*: Number of threads (default: 1)
      -p: Generate a pairwise dot plot for contigs in each cluster
    
    Example useage:
      /usr/local/bin/Assembly_toolkit/trycycler_run_reconcile.sh -i=2_clusters -n=3 -s=round1 -t=16 -p

    Dependencies: trycycler, perl.
    "
}

# Main ###############
if [ -z "$1" ]
then
    display_usage
    exit
fi

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
        -p)
        draw_plot=true
        ;;
        *)
        ;;
    esac
done

# Determine the number of clusters
if [ -d "$dir_in" ]
then
    n=$(ls -d -1 $dir_in/cluster_* | wc -l)
else
    echo "Error: input directory $dir_in was not found. Exit"
    exit 1
fi

# Run 'trycycler reconcile'
if [ -f "$reads" ]
then
    echo "Read file: $reads"
    echo "Input directory: $dir_in"
    echo "Number of clusters: $n"
    echo "Number of threads: $t"
    for i in `seq 1 ${n}`
    do
        c=$(printf "%03d" "$i")  # '001', '002', '003', ...
        echo "Reconcile contigs in cluster $c"
        if [ ! -z "$s" ]
        then
            log_filename="reconcile_cluster_${c}_${s}"  # Filenames do not start from "cluster_${c}" to avoid confusion with the cluster_* glob for the 'trycycler partition' command.
        else
            log_filename="reconcile_cluster_${c}"
        fi
        cluster_dir="$dir_in/cluster_${c}"
        trycycler reconcile --reads $reads --cluster_dir "$cluster_dir" --threads $t > "${log_filename}.tmp" 2>&1
        log_file="$dir_in/${log_filename}.log"
        perl -pe 's/\x1b\[[0-9;]*[mG]//g' "${log_filename}.tmp" > "$log_file"  # Remove colour code from the original log file
        rm "${log_filename}.tmp"
        echo "Cluster $c has been reconciled. Log file: ${log_file}."
        if [ "$draw_plot" == true ]
        then
            dotplot="$cluster_dir/cluster_${c}_dotplots.png"
            echo "Create a pairwise dotplot $dotplot"
            trycycler dotplot --cluster_dir "$cluster_dir"  # Create dotplots.png in the cluster directory
            mv "$cluster_dir/dotplots.png" "$dotplot"  # To make it easier to pool figures into a single directory
        fi
    done
else
    echo "Error: read file $reads was not found. Exit."
fi
