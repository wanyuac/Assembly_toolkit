#!/bin/bash
# Run Raven assembler for long-read assemblies
# run_raven.sh [isolate name] [input read file] [k-mer size] [window size] [number of polishing rounds] [Filename extension of the input read file (fastq/fastq.gz/etc)]
# Change the number of threads for your computer.
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 18 Dec 2022; latest update: 18 Dec 2022

# Read parameters ###############
i="$1"  # Isolate name
r="$2"  # Input read file
k="$3"  # K-mer size
w="$4"  # Window length
p="$5"  # Polish rounds
e="$6"  # Filename extension of the input read file (without '.'; For example, 'fastq')

# Run Raven assembler ###############
f=$(basename $r ".$e")
pre="${f}_raven_k${k}w${w}p${p}"
v=`raven --version`
echo "Assembling long reads $r using Raven $v (k=${k}, w=${w}; polish $p times)"
raven --threads 4 --polishing-rounds $p --kmer-len $k --window-len $w --disable-checkpoints --graphical-fragment-assembly ${pre}.gfa $r > ${pre}.fasta
