#!/usr/bin/env python

"""
This script adds modifiers to sequence headers according to an input TSV file for submission of complete
genome assemblies to GenBank.

The input TSV file must have a column 'file' or paths of input FASTA files are stored in the first column.

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 11 February 2024; the latest update: 20 February 2024.
"""

import os
import sys
import pandas as pd
from Bio import SeqIO
from argparse import ArgumentParser
from modules.utilities import write_seq  # A local module defined in modules/write_seq.py

def parse_arguments():
    parser = ArgumentParser(description = "Adding modifiers to sequence headers in FASTA files")
    parser.add_argument('-t', '--tsv', dest = 'tsv', type = str, required = True,\
                        help = "Input TSV file specifying modifiers for sequences in complete genome assemblies")
    parser.add_argument('-o', '--outdir', dest = 'outdir', type = str, required = False, default = '.', help = "Output directory (default: .)")
    parser.add_argument('-p', '--plasmid_prefix', dest = 'plasmid_prefix', type = str, required = False, default = 'p',\
                        help = "Indicator character for plasmid names (default: p)")
    return parser.parse_args()


def main():
    args = parse_arguments()
    plasmid_prefix = sanity_check(args.tsv, args.outdir, args.plasmid_prefix)
    tab, modifiers = import_tsv(args.tsv)  # Import the input TSV file and extract modifier names from column names
    seq_summary = open(os.path.join(args.outdir, 'sequence_summary.tsv'), 'w')
    seq_summary.write('\t'.join(['File', 'Contig', 'Type', 'Length_bp', 'Modifiers\n']))
    for i in range(len(tab)):  # Iterate through rows 0..21 of tab
        r = tab.iloc[i]  # Extract the i-th row from tab
        fasta_in = r['file']  # Input FASTA file
        if os.path.exists(fasta_in):
            fasta_basename = os.path.basename(fasta_in)
            fasta_out = os.path.join(args.outdir, fasta_basename)  # Output FASTA file
            if os.path.exists(fasta_out):
                print(f"Warning: Output file {fasta_out} exists and will be overwritten.", file = sys.stderr)
            with open(fasta_in, 'r') as fasta_in_handle:
                fasta_out_handle = open(fasta_out, 'w')
                for seq in SeqIO.parse(fasta_in_handle, "fasta"):  # Go through sequences in the input FASTA file
                    seq_name = seq.id
                    modifier_fields = [construct_seq_descr(r, m) for m in modifiers]
                    modifier_fields = list(filter(None, modifier_fields))  # Remove NA cells
                    if seq_name.startswith(plasmid_prefix):  # The current sequence is a plasmid
                        modifier_fields.append(f"[plasmid-name={seq_name}]")
                        seq_type = 'plasmid'
                    else:
                        modifier_fields.append(f"[location=chromosome]")  # Define the current sequence as a chromosome.
                        seq_type = 'chromosome'
                    s = str(seq.seq)
                    modifier_string = ' '.join(modifier_fields)
                    write_seq(seq_name = seq_name, seq_descr = modifier_string, seq = s, fasta_handle = fasta_out_handle)
                    seq_summary.write('\t'.join([fasta_basename, seq_name, seq_type, str(len(s)), modifier_string + '\n']))
                fasta_out_handle.close()
        else:
            print(f"Warning: input FASTA file {fasta_in} does not exist.", file = sys.stderr)
    seq_summary.close()
    return


def sanity_check(tsv, outdir, p):
    if not os.path.exists(tsv):
        print("Error: input TSV file " + tsv + " does not exist.", file = sys.stderr)
        sys.exit(1)
    if not os.path.exists(outdir):
        print("Create the output directory " + outdir, file = sys.stdout)
        os.mkdir(outdir)
    if p == '':
        p = 'p'
    return p


def import_tsv(tsv):
    tab = pd.read_csv(tsv, sep = '\t', index_col = None, na_values = ['NA'])
    colnames = tab.columns.to_list()
    if colnames[0] != 'file':
        tab = tab.rename(columns = {tab.columns[0] : 'file'})  # Rename the first column of the data frame to 'file'
        colnames[0] = 'file'
    return tab, colnames[1 : ]


def construct_seq_descr(r, m):
    """ r: a row in tab; m: modifier """
    v = r[m]
    if pd.isnull(v):
        y = ''
    else:
        y = f"[{m}={v}]"
    return y


if __name__ == '__main__':
    main()
