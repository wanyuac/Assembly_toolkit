#!/bin/bash
# Run Raven assembler for long-read assemblies
# run_raven.sh [input read file] [Filename extension of the input read file (fastq/fastq.gz/etc)] [k-mer size] [window size] [number of polishing rounds] [number of threads]
# Change the number of threads for your computer.
# Copyright (C) 2022-2023 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 18 Dec 2022; latest update: 3 Jan 2023

# Read parameters ###############
r="$1"  # Input read file
e="$2"  # Filename extension of the input read file, without the leading '.' sign (namely, fastq/fastq.gz/etc)
k="$3"  # K-mer size
w="$4"  # Window length
p="$5"  # Polish rounds
t="$6"  # Number of threads

# Run Raven assembler ###############
f=$(basename "$r" ".$e")
pre="${f}_raven_k${k}w${w}p${p}"
v=`raven --version`
if [ -z "$t" ]; then
    echo 'Set the number of threads to 1'
    t=1
fi
echo "Assembling long reads from $r using Raven $v (k=${k}, w=${w}; polish=${p}; threads=${t})"
raven --threads $t --polishing-rounds $p --kmer-len $k --window-len $w --disable-checkpoints --graphical-fragment-assembly ${pre}.gfa $r > ${pre}.fasta
