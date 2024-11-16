#!/bin/bash
# Rename bad contigs to exclude them before running 'trycycler reconcile'.
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 6 November 2024; last update: 12 November 2024

# Help information ###############
display_usage() {
    echo "
    This script renames bad contigs to exclude them from inputs of command \'trycycler reconcile\'.

    Parameters:
      -i=*: Input directory - the output directory of \'trcycyler cluster\'
      -t=*: A two-column TSV file listing cluster names (e.g., cluster_001) and filenames without .fasta (e.g., L_utg000002c)
    
    Example useage:
      ~/bin/Assembly_toolkit/trycycler_exclude_bad_contigs.sh -i=2_clusters -t=contigs_to_exclude.tsv
    
    Format of the TSV file:
    [cluster name]\t[comma-delimited contig names]
    "
}

# Main ###############
if [ -z "$1" ]
then
    display_usage
    exit
fi

for arg in "$@"
do
    case $arg in
        -i=*)
        dir_in="${arg#*=}"  # The input directory
        ;;
        -t=*)
        tsv="${arg#*=}"
        ;;
        *)
        ;;
    esac
done

while read line
do
    if [ ! -z "$line" ]
    then
        IFS=$'\t' read cluster_name contig_list <<< "$line"
        IFS=',' read -r -a contigs <<< "$contig_list"
        for contig in "${contigs[@]}"
        do
            fasta_source="$dir_in/${cluster_name}/1_contigs/${contig}.fasta"
            fasta_target="${fasta_source}.bad"
            if [ -f "$fasta_source" ]
            then
                echo "Rename $fasta_source as $fasta_target"
                mv "$fasta_source" "$fasta_target"
            fi
        done
    fi
done < "$tsv"
