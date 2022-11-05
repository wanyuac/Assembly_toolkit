#!/bin/bash
# Example: ./run_polypolish.sh isolate1 $asm_medaka $illumina $asm_poli
# Prerequisites: python v3, BWA
# Reference: https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 1 May 2022; latest update: 4 Nov 2022

# 1. Preparation ###############
i="$1"  # Isolate name
fasta_in="$2"  # Input FASTA file (path and filename)
reads_dir="$3"  # Directory of input Illumina reads
outdir="$4"  # Output directory
if [ ! -d "$outdir" ]; then
	mkdir $outdir
    mkdir ${outdir}/log
fi
if [ ! -d "${outdir}/log" ]; then
    mkdir "${outdir}/log"
fi
r1="${reads_dir}/${i}_1.fastq.gz"
r2="${reads_dir}/${i}_2.fastq.gz"
echo "Start to polish the assembly of $i (${fasta_in}) with reads from $reads_dir and save the output to ${outdir}/"

# 2. Filtering read alignments to exclude those of unusually large insert sizes ###############
echo "Creating SAM files for isolate $i"
tm="sams_$i"
mkdir $tm
bwa index $fasta_in
bwa mem -a -t 4 $fasta_in $r1 > ${tm}/unfiltered_1.sam
bwa mem -a -t 4 $fasta_in $r2 > ${tm}/unfiltered_2.sam
~/bin/polypolish-v0.5.0/polypolish_insert_filter.py --in1 ${tm}/unfiltered_1.sam --in2 ${tm}/unfiltered_2.sam --out1 ${tm}/filtered_1.sam --out2 ${tm}/filtered_2.sam 1>${outdir}/log/${i}_SAM_filter.log 2>&1

# 3. Polish the input assembly ###############
echo "Polishing assembly $fasta_in with Polypolish"
~/bin/polypolish-v0.5.0/polypolish $fasta_in ${tm}/filtered_1.sam ${tm}/filtered_2.sam 1>${outdir}/${i}_polypolish.fna 2>${outdir}/log/${i}_polypolish.log
rm -rf $tm
echo -e "Finished polishing the assembly of isolate ${i}\n"
