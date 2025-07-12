#!/bin/bash
# Extract lengths and depths of contigs from FASTA-format assemblies generated using Shovill or SPAdes.
#
# Usage:
#   ./extract_contig_stats.sh --assemblies=*.fna  # For Shovill assemblies
#   ./extract_contig_stats.sh --spades --assemblies=*.fna  # For SPAdes assemblies
#
# Copyright (C) 2024-2025 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Publication: 7 June 2024; latest update: 12 July 2025

shovill_mode=true
assemblies_glob=""

# Parse arguments ###############
for arg in "$@"; do
    case $arg in
        --spades)
            shovill_mode=false
            ;;
        --assemblies=*)
            assemblies_glob="${arg#*=}"
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            exit 1
            ;;
    esac
done

if [ -z "$assemblies_glob" ]; then
    echo "Error: --assemblies parameter is required." >&2
    exit 1
fi

# Expand the file list ###############
shopt -s nullglob  # Enable "nullglob" to make unmatched globs disappear instead of staying literal

# Expands the assemblies_glob (e.g. assemblies/*.fna) into an array of actual file paths;
# if no files matched and nullglob is enabled, assemblies will be an empty array.
# Without nullglob, assemblies would be ("assemblies/*.fna").
assemblies=($assemblies_glob)

shopt -u nullglob  # Restores Bash's default behaviour to avoid affecting the global shell environment beyond this point

if [ ${#assemblies[@]} -eq 0 ]; then
    echo "Error: No matching FASTA files found for '$assemblies_glob'" >&2
    exit 1
fi

# Main ###############
if [ "$shovill_mode" = true ]; then
    # Shovill format: >[contig name] len=[length] cov=[depth] ...
    for fasta in "${assemblies[@]}"; do
        fasta_filename=$(basename "$fasta")
        isolate="${fasta_filename%.*}"  # Remove the filename extension. In Shell's parameter expansion syntax, ${variable%pattern} removes the shortest match of the pattern from the end of the variable's value.
        output_tsv="${isolate}__contig_stats.tsv"
        : > "$output_tsv"  # A concise Bash idiom that cleans (truncates) or creates the output file to ensure it exists and is empty before appending to it
        grep -F '>' "$fasta" | while read -r seq_header; do
            fields=$(echo "$seq_header" | cut -d ' ' -f 1-3)  # Extract the first three fields (contig name, length, and coverage) from the current seq_header
            IFS=' ' read -r -a fields_array <<< "$fields"  # Parse the fields string into an array
            contig_name=${fields_array[0]#>}  # Use the parameter expansion syntax (`${variable#pattern}`) to remove the leading '>' from the first field.
            contig_len=${fields_array[1]#*=}
            contig_depth=${fields_array[2]#*=}
            echo -e "${isolate}\t${contig_name}\t${contig_len}\t${contig_depth}" >> "$output_tsv"
        done
    done
else
    # SPAdes header format example: >[contig name]_length_[length]_cov_[depth]
    for fasta in "${assemblies[@]}"; do
        fasta_filename=$(basename "$fasta")
        isolate="${fasta_filename%.*}"
        output_tsv="${isolate}__contig_stats.tsv"
        : > "$output_tsv"
        grep -F '>' "$fasta" | while read -r seq_header; do
            # contig_name=$(echo "$seq_header" | sed -n 's/^>\(.*\)_length.*/\1/p')
            contig_name="${seq_header#>}"  # It is desirable to keep the complete, unparsed contig names for the ease of contig filtering later.
            contig_len=$(echo "$seq_header" | sed -n 's/.*_length_\([0-9]\+\)_cov_.*/\1/p')
            contig_depth=$(echo "$seq_header" | sed -n 's/.*_cov_\([0-9.]\+\)/\1/p')
            echo -e "${isolate}\t${contig_name}\t${contig_len}\t${contig_depth}" >> "$output_tsv"
        done
    done
fi

# Concatenate all individual TSVs into one file with a single header
echo -e "Isolate\tContig\tLength\tDepth" > all_contig_stats.tsv
cat *__contig_stats.tsv >> all_contig_stats.tsv