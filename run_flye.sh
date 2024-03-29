#!/bin/bash
# Run Flye to assemble (Oxford Nanopore Technologies) ONT reads
# run_flye.sh [input read file] [Filename extension of the input read file (fastq/fastq.gz/etc)] [probable genome size] [number of polishing rounds] [number of threads] [error level of reads] "scaffold"
# Change the number of threads for your computer.
# Copyright (C) 2022-2023 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 18 Dec 2022; latest update: 16 Aug 2023

# Read parameters ###############
r="$1"  # Input read file (filename extensions: fastq, fastq.gz, etc)
e="$2"  # Filename extension of the input read file (fastq/fastq.gz/etc)
g="$3"  # Probable genome size (for instance, 5m)
p="$4"  # Number of polishing rounds
t="$5"  # Number of threads
m="$6"  # Mode of Nanopore reads: raw (expected error rate: <15%) or hq (<5%)
s="$7"  # Value "scaffold" or any non-null values to enable scaffolding

# Generate a genome assembly from long reads ###############
f=$(basename $r ".$e")
pre="${f}_flye_p${p}"
v=`flye --version`

if [ -z "$t" ]; then
    t=1;
fi

if [ -z "$m" ]; then
    m="raw"
fi

if [ ! -z "$s" ]; then
    scaf=true
else
    scaf=false
fi

mkdir tmp

if [ "$m" == "hq" ]; then
    echo "Assembling high-quality (error rate <5%) ONT reads $r using Flye v$v (polish: ${p}) with $t threads"
    cmd="--nano-hq $r --threads $t --out-dir tmp --iterations $p --genome-size $g"  # --subassemblies is incompatible with --nano-raw or --nano-corr
else
    echo "Assembling regular (error rate <15%) ONT reads $r using Flye v$v (polish: ${p}) with $t threads"
    cmd="--nano-raw $r --threads $t --out-dir tmp --iterations $p --genome-size $g"
fi

if [ "$scaf" == true ]; then
    echo "Enable scaffolding"
    cmd="$cmd --scaffold"
fi

echo "Run command flye $cmd"
flye $cmd

if [ -d "tmp" ]; then
    mv tmp/assembly.fasta ${pre}.fasta
    mv tmp/assembly_graph.gfa ${pre}.gfa
    mv tmp/assembly_info.txt ${pre}.txt
    mv tmp/flye.log ${pre}.log
    rm -rf tmp
fi
