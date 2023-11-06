#!/bin/bash
# run_polypolish.sh [isolate name] [input FASTA file] [directory of short reads] [output directory] [prefix of output files]
# The prefix of outputs can be "1_${sample_name}" or "2_${sample_name}", etc.
# Prerequisites: python v3, BWA
# Reference: https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 1 May 2022; latest update: 6 Nov 2023

# 1. Preparation ###############
# Read parameters
i="$1"  # Isolate name
fasta_in="$2"  # Input FASTA file (path and filename)
reads_dir="$3"  # Directory of input Illumina reads
outdir="$4"  # Output directory
t="$5"  # Number of threads
p="$6"  # Prefix of output files

# Determine the number of threads
if [ -z "$t" ]
then
    echo 'Setting number of threads to 4'
    t=4
else
    echo "Number of threads: $t"
fi

# Set up output filenames
if [ ! -z "$p" ]
then
    sam_filter_prefix="${p}_SAM_filter"
    polypolish_prefix="${p}_polypolish"
else
    sam_filter_prefix="${i}_SAM_filter"
    polypolish_prefix="${i}_polypolish"
fi

# Set up input and output paths
r1="${reads_dir}/${i}_1.fastq.gz"
r2="${reads_dir}/${i}_2.fastq.gz"

if [ ! -d "$outdir" ]
then
	mkdir $outdir
    mkdir ${outdir}/log
fi

if [ ! -d "${outdir}/log" ]
then
    mkdir "${outdir}/log"
fi

echo "$(date): Start to polish the assembly of $i (${fasta_in}) with reads from $reads_dir and save outputs in ${outdir}/"

# 2. Filtering read alignments to exclude those of unusually large insert sizes ###############
echo "Creating SAM files for isolate $i"
tm="sams_$i"
mkdir $tm
bwa index $fasta_in
bwa mem -a -t $t $fasta_in $r1 > ${tm}/unfiltered_1.sam
bwa mem -a -t $t $fasta_in $r2 > ${tm}/unfiltered_2.sam
polypolish_insert_filter.py --in1 ${tm}/unfiltered_1.sam --in2 ${tm}/unfiltered_2.sam --out1 ${tm}/filtered_1.sam --out2 ${tm}/filtered_2.sam 1>${outdir}/log/${sam_filter_prefix}.txt 2>&1

# 3. Polish the input assembly ###############
echo "$(date): Polishing assembly $fasta_in with Polypolish"
polypolish $fasta_in ${tm}/filtered_1.sam ${tm}/filtered_2.sam 1>${outdir}/${polypolish_prefix}.fna 2>${outdir}/log/${polypolish_prefix}.txt
rm -rf $tm

# 4. Remove '_polypolish' from sequence headers to preserve the original ones ###############
echo "$(date): processing the polished assembly and log files"
sed -i 's/_polypolish//g' ${outdir}/${polypolish_prefix}.fna

# 5. Remove non-character content from log files ###############
perl -p -e 's/\x1b\[[0-9;]*[mG]//g' ${outdir}/log/${polypolish_prefix}.txt > ${outdir}/log/${polypolish_prefix}.log
perl -p -e 's/\x1b\[[0-9;]*[mG]//g' ${outdir}/log/${sam_filter_prefix}.txt > ${outdir}/log/${sam_filter_prefix}.log
rm ${outdir}/log/${polypolish_prefix}.txt ${outdir}/log/${sam_filter_prefix}.txt

echo "$(date): Finished polishing the assembly of isolate ${i}"
