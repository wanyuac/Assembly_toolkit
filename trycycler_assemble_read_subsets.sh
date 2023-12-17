#!/bin/bash
# Run Flye, Raven, and Minipolish assembler to assemble subsets of Nanopore reads generated using command "trycycler subsample".
# Copyright (C) 2023 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 27 July 2023; last update: 17 December 2023

# Help information ###############
display_usage() {
    echo "
    This script iteratively runs Raven, Flye, and Miniasm-minipolish to assemble subsets (<100) of Nanopore reads
    generated using command 'trycycler subset'. Output assemblies (in subdirectory '1_assemblies' in the current
    working directory) can be used for Trycycler command 'trycycler cluster'.
    
    Preprequsites: please ensure flye, raven, minimap2, miniasm, minipolish, and any2fasta are accessible in
    your environment. These programs can be installed with Conda.

    Parameters:
      -d=*: directory of input FASTQ files (*.fastq, mandatory, and do not include the final forward slash)
      -p=*: number of assembly-polishing rounds (default: 2; optional)
      -t=*: number of threads (default: 1; optional)
      -l=*: expected genome length (e.g., 5m) for Flye (mandatory)
      -h: a flag to instruct Flye to treat input reads as of high quality (error rate <5%; optional)
      -k=*: length of minimisers used to find overlaps by Raven (default: 15; optional)
      -w=*: length of sliding window from which minimisers are sampled by Raven (default: 5; optional)
      -s=*: start sample index (default: 1; optional); useful when adding more fastq files into a previous run
    
    Example useage:
      /usr/local/bin/Assembly_toolkit/trycycler_assemble_read_subsets.sh -d=\"\$HOME/reads/subsets\" -p=2 -t=8 -l=2.5m -h -k=15 -w=5
      /usr/local/bin/Assembly_toolkit/trycycler_assemble_read_subsets.sh -d=\"\$HOME/reads/subsets_extra\" -p=2 -t=8 -l=2.5m -h -k=15 -w=5 -s=13

    Output: directory '1_assemblies', which stores assemblies in GFA and FASTA formats, will be created under the current working directory.

    Note: when parameter s > 1, the directory of input reads must contain corresponding FASTQ files. For example, when s = 16, the directory
    must have sample_16.fastq, etc. Otherwise, assemblers will return an error of missing inputs.
    "
}

output_dir='1_assemblies'

# Assembler runners ###############
run_flye() {
    local i="$1"  # Index of read subset: 01, 02, ..., 12, etc
    local d="$2"  # Directory of input reads
    local g="$3"  # Probable genome size (for instance, 5m)
    local p="$4"  # Number of polishing rounds
    local t="$5"  # Number of threads
    local h="$6"  # Accuracy level of Nanopore reads: raw (expected error rate: <15%) or hq (<5%)
    local r="$d/sample_${i}.fastq"  # Input FASTQ file
    local pre="${output_dir}/assembly_${i}"

    # Generate a genome assembly from long reads ###############
    tmp_dir=$(mktemp -d -t flye-XXXXXXXXXX)  # Create a temporary directory
    
    if [ "$h" == true ]
    then
        echo "[$(date)] Assembling high-quality (error rate <5%) ONT reads $r using Flye (polish: ${p}) with $t threads (temporary output directory: $tmp_dir)"
        flye --nano-hq "$r" --threads "$t" --out-dir "$tmp_dir" --iterations "$p" --genome-size "$g"
    else
        echo "[$(date)] Assembling regular (error rate <15%) ONT reads $r using Flye v$v (polish: ${p}) with $t threads (temporary output directory: $tmp_dir)"
        flye --nano-raw "$r" --threads "$t" --out-dir "$tmp_dir" --iterations "$p" --genome-size "$g"
    fi

    mv $tmp_dir/assembly.fasta ${pre}.fasta
    mv $tmp_dir/assembly_graph.gfa ${pre}.gfa
    mv $tmp_dir/assembly_info.txt ${pre}.txt
    mv $tmp_dir/flye.log ${pre}.log
    rm -rf $tmp_dir
}

run_raven() {
    local i="$1"  # Index of read subset: 01, 02, ..., 12, etc
    local d="$2"  # Directory of input reads
    local k="$3"  # K-mer size
    local w="$4"  # Window length
    local p="$5"  # Polish rounds
    local t="$6"  # Number of threads
    local r="$d/sample_${i}.fastq"
    local pre="${output_dir}/assembly_${i}"
    echo "[$(date)] Assembling long reads from FASTQ file $r using Raven (k=${k}, w=${w}; polish=${p}; threads=${t})"
    raven --threads "$t" --polishing-rounds "$p" --kmer-len "$k" --window-len "$w" --disable-checkpoints --graphical-fragment-assembly "${pre}.gfa" "$r" > "${pre}.fasta"
}

run_miniasm_and_minipolish() {  # Code in the function is adapted from https://github.com/rrwick/Minipolish/blob/main/miniasm_and_minipolish.sh.
    # It takes two positional arguments:
    #  1) a long read file
    #  2) the number of threads to use
    local i="$1"  # Index of read subset: 01, 02, ..., 12, etc
    local d="$2"  # Directory of input reads
    local t="$3"  # Number of threads
    local r="$d/sample_${i}.fastq"  # Input FASTQ file
    local pre="${output_dir}/assembly_${i}"  # Prefix of output files

    # Create temporary intermediate files.
    overlaps=$(mktemp)".paf"
    unpolished_assembly=$(mktemp)".gfa"

    # Assemble ONT reads and polish
    echo "[$(date)] Assembling reads from $r using miniasm_and_minipolish ($t threads)"
    minimap2 -x ava-ont -t "$t" "$r" "$r" > "$overlaps"  # Find read overlaps with minimap2
    miniasm -f "$r" "$overlaps" > "$unpolished_assembly"  # Run miniasm to make an unpolished assembly
    minipolish --threads "$t" "$r" "$unpolished_assembly" > "${pre}.gfa"  # Polish the assembly with minipolish, outputting the result to stdout.
    any2fasta "${pre}.gfa" > "${pre}.fasta"  # Convert the GFA file to a FASTA file
    rm "$overlaps" "$unpolished_assembly"  # Clean up
}

# Main ###############
if [ -z $1 ]
then
    display_usage
    exit
fi

# Read arguments ===============
index_start=1  # Default value

for i in "$@"
do
    case $i in
        -d=*)
        dir_in="${i#*=}"  # Directory of input read files sample_01.fastq, sample_02.fastq, ..., sample_xy.fastq
        ;;
        -p=*)
        polish="${i#*=}"  # Number of polishing rounds
        ;;
        -t=*)
        threads="${i#*=}"  # Number of threads
        ;;
        -l=*)
        genome_len="${i#*=}"  # Probable genome length (bp) for Flye
        ;;
        -h)
        high_accuracy_reads=true  # High-accuracy ONT reads (error rate <5%) for Flye
        ;;
        -k=*)
        raven_kmer="${i#*=}"  # K-mer length for Raven (default: 15)
        ;;
        -w=*)
        raven_window_len="${i#*=}"  # Window length for Raven (default: 5)
        ;;
        -s=*)
        index_start="${i#*=}"  # Update this value with the user-specified start index of FASTQ filenames
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

# Initialisation ===============
if [ ! -d "$dir_in" ]
then
    echo "[$(date)] Error: directory of input FASTQ files was not found."
    exit
fi

if [ -z "genome_len" ]
then
    echo "[$(date)] Error: expected genome length was not specified for Raven."
    exit
fi

if [ -z "$high_accuracy_reads" ]
then
    echo "[$(date)] Flye: treats input reads as raw (error rate <15%)."
    high_accuracy_reads=false
else
    echo "[$(date)] Flye: treats input reads as high-quality (error rate <5%)"
fi

if [ -z "$polish" ]
then
    echo "[$(date)] Set the number of assembly-polishing rounds to 2."
    polish=2
fi

if [ -z "$threads" ]
then
    echo "[$(date)] Thread number was not specified, so set the number to 1."
    threads=1
fi

if [ -z "$raven_kmer" ]
then
    echo "[$(date)] Set k-mer length to 15 bp for Raven."
    raven_kmer=15
fi

if [ -z "$raven_window_len" ]
then
    echo "[$(date)] Set window size to 5 for Raven."
    raven_window_len=5
fi

if [ ! -d "$output_dir" ]
then
    echo "[$(date)] Output directory $output_dir was not found, so create it in the current working directory."
    mkdir "$output_dir"
fi

# Assemble subsets of ONT reads ===============
n=$(find $dir_in -name '*.fastq' -type f | wc -l)

if [ $n -gt 0 ]
then
    echo "[$(date)] Found $n input FASTQ files"
    j=$index_start
    m=$((j+n))  # The maximum index + 1
    echo "FASTQ filenames start from index $j"
else
    echo "[$(date)] Error: no FASTQ file (.fastq) was found in input directory $dir_in".
    exit
fi

while [ $j -lt $m ]
do
    k=$(printf "%02d" $j)  # Add leading zeros to the index when j < 10
    run_flye "$k" "$dir_in" "$genome_len" "$polish" "$threads" "$high_accuracy_reads"
    ((j++))
    if [ $j -lt $m ]
    then
        k=$(printf "%02d" "$j")
        run_raven "$k" "$dir_in" "$raven_kmer" "$raven_window_len" "$polish" "$threads"
        ((j++))
        if [ $j -lt $m ]
        then
            k=$(printf "%02d" "$j")
            run_miniasm_and_minipolish "$k" "$dir_in" "$threads"
            ((j++))
        fi
    fi
done
