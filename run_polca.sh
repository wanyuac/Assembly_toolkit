#!/bin/bash
# Copyright (C) 2022-2023 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# First version: 1 May 2022; latest update: 17 Dec 2023

SCRIPT_VERSION=1.1.0

# Function definitions ###############
# Run './run_polca.sh' to display the information below.
display_parameters() {
    echo "
    run_polca.sh v$SCRIPT_VERSION

    This script polishes an input genome assembly with POLCA and paired-end short reads.

    Prerequisites: please ensure bwa aligner is accessible in PATH.

    Parameters (seven in total):
      -a=*: path and filename of the input assembly in FASTA format (mandatory)
      -r=*: path to the directory of paired-end short reads (without the end forward
            slash), where read files' names follow format [sample name]_[1,2].fastq.gz (mandatory)
      -i=*: isolate name (default: isolate)
      -n=*: prefix of output filenames (default: isolate)
      -o=*: path to output directory (default: polca_output)
      -p=*: MaSuRCA/bin directory (without the end forward slash), where polca.sh is
            stored (mandatory)
      -t=*: number of threads (default: 2)
    
    Example command:
      /usr/local/bin/Assembly_toolkit/run_polca.sh -a=1_isolate1_polypolish.fasta -r=reads/illumina \\
        -i=isolate1 -n=3_isolate1 -o="$PWD" -p=$HOME/bin/MaSuRCA-4.0.5/bin -t=8 > 3_isolate1_polca.log
    
    Output: a polished assembly [o]/[n]_polca.fna in FASTA format
    "
}

MYDIR="$(dirname "$(readlink -f "$0")")"  # Directory of this script
source $MYDIR/modules.sh  # Function print_failure_message

# Print parameter information ###############
if [ -z "$1" ]
then
    display_parameters
    exit
fi

# Set default values ###############
i=isolate
n="$i"
outdir=polca_output
t=2

# Read parameters ###############
for arg in "$@"
do
    case $arg in
        -a=*)
        fasta_in="${arg#*=}"
        ;;
        -r=*)
        read_dir="${arg#*=}"
        ;;
        -i=*)
        i="${arg#*=}"
        ;;
        -n=*)
        n="${arg#*=}"
        ;;
        -o=*)
        outdir="${arg#*=}"
        ;;
        -p=*)
        polca_dir="${arg#*=}"
        ;;
        -t=*)
        t="${arg#*=}"
        ;;
        *)
        ;;
    esac
done

# Run POLCA for the input genome assembly ###############
if [ -z "$(which bwa)" ] 
then
    echo "Error: bwa was not accessible"
    print_failure_message "$i"
    exit
fi

if [ -e "$fasta_in" ]
then
    fasta_name=`basename $fasta_in`
    echo "Base name of the input FASTA file: $fasta_name"
else
    echo "Error: assembly $fasta_in was not found."
    exit
fi

if [ -d "$read_dir" ]
then
    r1=${read_dir}/${i}_1.fastq.gz
    r2=${read_dir}/${i}_2.fastq.gz
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]  # '-f' works for symbolic links
    then
        echo "Error: read file $r1 and/or $r2 were not found."
        print_failure_message "$i"
        exit
    fi
else
    echo "Error: read directory $read_dir was not found."
    print_failure_message "$i"
    exit
fi

if [ -d "$polca_dir" ]
then
    polca="$polca_dir/polca.sh"
    if [ -f "$polca" ]
    then
        echo "run_polca.sh v$SCRIPT_VERSION"
        if [ ! -d "$outdir" ]
        then
            echo "Create output directory $outdir"
            mkdir $outdir
        fi
        tm="polca_tmp_$i"
        mkdir $tm  # Create a temporary directory in the current working directory
        cd $tm
        echo "$(date): Start to polish $fasta_in of isolates $i with reads $r1 and $r2 and with $t threads"
        $polca -a $fasta_in -r "$r1 $r2" -t "$t" -m 8G
        polca_out="${fasta_name}.PolcaCorrected.fa"
        fasta_out="${outdir}/${n}_polca.fna"
        if [ -f "$polca_out" ]
        then
            cd ..
            echo "Saving $tm/$polca_out to $fasta_out"
            mv "$tm/$polca_out" "$fasta_out"  # Output FASTA file and its name
            rm -r $tm  # Delete the temporary directory and its content
            rm "$fasta_in".fai
            echo "Success: polished assembly of isolate $i was saved as ${fasta_out}."
        else
            echo "Error: output file $polca_out was not found. Please check files in $tm for details."
            print_failure_message "$i"
        fi
    else
        echo "Error: script $polca was not found."
        print_failure_message "$i"
    fi
else
    echo "Error: POLCA directory $polca_dir was not found."
    print_failure_message "$i"
fi
