#!/bin/bash
# Run miniasm and minipolish for long-read genome assembly
# Usage:
#    run_miniasm_and_minipolish [isolate or genome name] [path to long reads (fast.gz)] [number of threads]
# Dependencies: minimap2, miniasm, minipolish, any2fasta
#
# Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 29 Janurary 2024; latest update: 16 August 2024

i="$1"  # Isolate name
r="$2"  # Input long reads
t="$3"  # Number of threads

# Create temporary intermediate files.
overlaps=$(mktemp)".paf"
unpolished_assembly=$(mktemp)".gfa"

# Assemble ONT reads and polish
echo "[$(date)] Assembling reads from $r using miniasm_and_minipolish ($t threads)"
minimap2 -x ava-ont -t "$t" "$r" "$r" > "$overlaps"  # Find read overlaps with minimap2
miniasm -f "$r" "$overlaps" > "$unpolished_assembly"  # Run miniasm to make an unpolished assembly
minipolish --threads "$t" "$r" "$unpolished_assembly" > "${i}.gfa"  # Polish the assembly with minipolish, outputting the result to stdout.
any2fasta "${i}.gfa" > "${i}.fasta"  # Convert the GFA file to a FASTA file
rm "$overlaps" "$unpolished_assembly"  # Clean up
