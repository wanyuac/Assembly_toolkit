#!/usr/bin/env python

"""
This script identifies circular sequences in an assembly graph and draws a dot plot for
each of them.

Usage:
    dotplots_for_circular_seqs.py --graph assembly.gfa --outdir dotplots --sample isolate_A --script_dir ~/bin/Assembly_toolkit

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
First edition: 11 November 2024; latest update: 11 November 2024
"""

import os
import shutil
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # Determine where this script is located in the file system

def parse_arguments():  # https://stackoverflow.com/questions/7066826/in-python-how-to-get-subparsers-to-read-in-parent-parsers-argument
    parser = argparse.ArgumentParser(description = "Remove contigs from FASTA files by contig names")
    parser.add_argument('-g', '--graph', dest = 'graph', type = str, required = True, help = "Input assembly graph in the GFA format")
    parser.add_argument('-o', '--outdir', dest = 'outdir', type = str, required = False, default = 'dotplots',\
                        help = "Path to the output directory (default: dotplots)")
    parser.add_argument('-n', '--name', dest = 'name', type = str, required = False, default = 'sample1',\
                        help = "Sample name for output files (default: sample1)")
    parser.add_argument('-s', '--script_dir', dest = 'script_dir', type = str, required = False, default = SCRIPT_DIR, \
                        help = f"Path to the Assembly_toolkit directory (default: {SCRIPT_DIR})")
        
    return parser.parse_args()


def main():

    return


if __name__ == '__main__':
    main()
