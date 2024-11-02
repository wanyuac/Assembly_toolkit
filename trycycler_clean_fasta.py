#!/usr/bin/env python

"""
This script removes contigs from a FASTA file by names for preparing inputs of the 'trycycler cluster'
command.

Usage:
    python trycycler_clean_fasta.py --input_dir [directory of FASTA files] --excl_tab [exclusion table]

The exclusion table is a header-less TSV file of two columns: (1) name of each Input FASTA file and
(2) comma-delimited names of contigs to exclude.

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
First edition: 2 November 2024; latest update: 2 November 2024
"""

import os
import shutil
import argparse

def parse_arguments():  # https://stackoverflow.com/questions/7066826/in-python-how-to-get-subparsers-to-read-in-parent-parsers-argument
    parser = argparse.ArgumentParser(description = "Remove contigs from FASTA files by contig names")
    parser.add_argument('-i', '--input_dir', dest = 'input_dir', type = str, required = True, help = "Directory of input FASTA files")
    parser.add_argument('-e', '--excl_tab', dest = 'excl_tab', type = str, required = True, help = "File path of an exclusion table")
    parser.add_argument('-o', '--ori_dir', dest = 'ori_dir', type = str, required = False, default = 'original', \
                        help = "Name of the directory to keep original FASTA files that are listed in the exclusion table")
    parser.add_argument('-r', '--rm_dir', dest = 'rm_dir', type = str, required = False, default = 'excluded', \
                        help = "Name of the directory to store excluded contigs")
    
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Preparation
    input_dir = args.input_dir
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Error: Input directory {input_dir} was not found.")
    org_dir = os.path.join(input_dir, args.ori_dir)
    rm_dir = os.path.join(input_dir, args.rm_dir)
    for d in [org_dir, rm_dir]:
        check_dir(d)
    
    # Filter FASTA files
    excl_tab = read_exclusion_table(args.excl_tab)
    for fasta_filename, contig_names in excl_tab.items():
        filter_fasta(fasta_filename, contig_names, input_dir, org_dir, rm_dir)
    return


def check_dir(d):
    if not os.path.exists(d):
        print(f"Create directory {d}")
        os.makedirs(d)
    return


def read_exclusion_table(e):
    """ Import the exclusion table as a dictionary """
    if not os.path.exists(e):
        raise FileNotFoundError(f"Exclusion table '{e}' was not found.")
    exclusion_table = {}
    with open(e, 'r') as tab_file:
        for line in tab_file:  # The line includes the newline character.
            line = line.strip()  # This step ensures contig names do not contain any newline character. It also removes empty lines.
            fasta, contig_names = line.strip().split('\t')
            exclusion_table[fasta] = contig_names.split(',')
    return exclusion_table


def filter_fasta(fasta_filename, contig_names, input_dir, ori_dir, rm_dir):
    """ Remove certain contigs from a FASTA file """
    fasta = os.path.join(input_dir, fasta_filename)
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"Input FASTA file {fasta} was not found.")
    
    input_fasta = os.path.join(ori_dir, fasta_filename)
    try:
        shutil.move(fasta, input_fasta)
    except Exception as err:
        raise RuntimeError(f"Error moving file {fasta} to {input_fasta}: {err}")

    with open(input_fasta, 'r') as fasta_file, \
        open(fasta, 'w') as output_fasta, \
        open(os.path.join(rm_dir, fasta_filename.split('.')[0] + '_excluded.fna'), 'w') as excl_fasta:
        contig_accept = False
        for line in fasta_file:
            if line.startswith('>'):  # Read a header line
                header_line = line.strip()  # This step removes the newline character from the contig name when there is no other fields in this line (e.g., Flye's outputs).
                contig = header_line[1:].split(' ', maxsplit = 1)[0]  # Name of the current contig
                contig_accept = contig not in contig_names
                if not contig_accept:
                    print(f"Exclude contig {contig} from {input_fasta}")
            if contig_accept:
                output_fasta.write(line)
            else:
                excl_fasta.write(line)
    
    return

    
if __name__ == '__main__':
    main()
    