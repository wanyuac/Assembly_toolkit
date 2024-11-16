#!/bin/bash
# Extract lengths and depths of contigs from FASTA-format assemblies generated using Shovill.
# Usage:
#   ./extract_shovill_contig_stats.sh *.fna
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 7 June 2024; latest update: 7 June 2024

for fasta in "$@"
do
    fasta_filename=$(basename $fasta)
    isolate="${fasta_filename%.*}"  # Remove the filename extension. In Shell's parameter expansion syntax, ${variable%pattern} removes the shortest match of the pattern from the end of the variable's value.
    output_tsv="${isolate}__contig_stats.tsv"
    if [ -f "$output_tsv" ]
    then
        rm -f "$output_tsv"  # Delete the output file if it exists
    fi
    touch "$output_tsv"
    grep '>' "$fasta" | while read -r seq_header
    do
        fields=$(echo "$seq_header" | cut -d ' ' -f 1-3)  # Extract the first three fields (contig name, length, and coverage) from the current seq_header
        IFS=' ' read -r -a fields_array <<< "$fields"  # Parse the fields string into an array
        contig_name=${fields_array[0]#>}  # Use the parameter expansion syntax (`${variable#pattern}`) to remove the leading '>' from the first field.
        contig_len=${fields_array[1]#*=}
        contig_depth=${fields_array[2]#*=}
        echo -e "${isolate}\t${contig_name}\t${contig_len}\t${contig_depth}" >> $output_tsv
    done
done
