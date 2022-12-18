#!/bin/bash
# Run Flye assembler for long-read assemblies
# run_flye.sh [isolate name] [input read file] [number of polishing rounds] [Filename extension of the input read file (fastq/fastq.gz/etc)]
# Change the number of threads for your computer.
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 18 Dec 2022; latest update: 18 Dec 2022

# Read parameters ###############
i="$1"  # Isolate name
r="$2"  # Input read file (filename extensions: fastq, fastq.gz, etc)
p="$3"  # Number of polishing rounds
e="$4"  # Filename extension of the input read file (fastq/fastq.gz/etc)
g="$5"  # Genome size (for instance, 5m)

# Generate a genome assembly from long reads ###############
f=$(basename $r ".$e")
pre="${f}_flye_p${p}"
v=`flye --version`
echo "Assembling long reads $r using Flye v$v (polish: ${p})"
mkdir ./tmp
flye --nano-raw $r --threads 4 --out-dir ./tmp --iterations $p --genome-size $g --scaffold  # --subassemblies is incompatible with --nano-raw or --nano-corr
mv ./tmp/assembly.fasta ${pre}.fasta
mv ./tmp/assembly_graph.gfa ${pre}.gfa
mv ./tmp/assembly_info.txt ${pre}.txt
mv ./tmp/flye.log ${pre}.log
rm -rf ./tmp
