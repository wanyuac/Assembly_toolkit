#!/bin/bash
# run_polca.sh [isolate name] [input FASTA file] [Directory of short reads] [output directory] [Polca's directory]
# Prerequisites: python v3, BWA
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 1 May 2022; latest update: 24 July 2023

# Read parameters ###############
i=$1
fasta_in=$2
reads_dir=$3
outdir=$4
polca_dir=$5  # Where polca.sh is found (no last forward slash)

# Run POLCA for the input genome ###############
fasta_name=`basename $fasta_in`
r1=${reads_dir}/${i}_1.fastq.gz
r2=${reads_dir}/${i}_2.fastq.gz
if [ -f "$r1" ] && [ -f "$r2" ]; then
    tm="polca_tmp_$i"
    mkdir $tm
    cd $tm
    echo "Polishing assembly $fasta_in for isolates $i with reads from $reads_dir"
    $polca_dir/polca.sh -a $fasta_in -r "$r1 $r2" -t 4 -m 6G
    fasta_out="${fasta_name}.PolcaCorrected.fa"
    if [ -f "$fasta_out" ]; then
        cd ..
        echo "Moving $tm/$fasta_out to ${outdir}/${i}_polca.fna"
        mv $tm/$fasta_out ${outdir}/${i}_polca.fna  # Output FASTA file and its name
        rm $tm/*.*
        rmdir $tm
    else
        cd ..
        echo "Error: the assembly of $i could not be polished."
    fi
else
    echo "Skipped isolates $i for the absence of its reads."
fi
