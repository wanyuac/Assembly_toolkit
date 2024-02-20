#!/usr/bin/env python

"""
This script extracts plasmid and chromosome sequences from each input multi-FASTA file of a complete bacterial genome assembly
and saves each sequence in a separate FASTA file with a name [isolate name]__[chromosome/plasmid]__[chromosome/plasmid name].fna.
This script extract isolate names from names of input FASTA files, assuming the names follow the format [isolate name].[filename extension, fna/fasta].

Example command:
    extract_individual_sequences.py -i *.fna -e fna -o . -c sequence_counts.tsv -l sequence_lengths.tsv -p p

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 20 February 2024; the latest update: 20 February 2024.
"""

import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser
from modules.utilities import check_dir, write_seq

def parse_arguments():
    parser = ArgumentParser(description = "Adding modifiers to sequence headers in FASTA files")
    # Inputs
    parser.add_argument('-i', '--inputs', dest = 'inputs', nargs = '+', type = str, required = True, help = "One or multiple input multi-FASTA files")
    parser.add_argument('-e', '--filename_extension', dest = 'filename_extension', type = str, required = False, default = 'fna',\
                        help = "Filename extension of input FASTA files (default: fna)")
    # Outputs
    parser.add_argument('-o', '--outdir', dest = 'outdir', type = str, required = False, default = '.', help = "Output directory (default: .)")
    parser.add_argument('-l', '--lengths', dest = 'lengths', type = str, required = False, default = 'sequence_lengths.tsv',\
                        help = "Filename of output tab-delimited lengths of sequences per isolate (default: sequence_lengths.tsv in outdir)")
    parser.add_argument('-c', '--counts', dest = 'counts', type = str, required = False, default = 'sequence_counts.tsv',\
                        help = "Filename of output tab-delimited sequence counts per isolate (default: sequence_counts.tsv in outdir)")
    # Parameter for differentiating between chromosome and plasmid sequences
    parser.add_argument('-p', '--plasmid_prefix', dest = 'plasmid_prefix', type = str, required = False, default = 'p',\
                        help = "Indicator character for plasmid names (default: p)")
    return parser.parse_args()


def main():
    args = parse_arguments()
    check_dir(args.outdir)
    print("Read " + str(len(args.inputs)) + " FASTA files as input.")
    seq_counts = open(os.path.join(args.outdir, args.counts), 'w')
    seq_lengths = open(os.path.join(args.outdir, args.lengths), 'w')
    seq_counts.write("Isolate\tChromosome\tPlasmids\n")
    seq_lengths.write("Isolate\tSequence\tType\tLength\n")
    for f in args.inputs:
        i = extract_isolate_name(f, args.filename_extension)
        num_chr = 0
        num_plasmids = 0
        if os.path.isfile(f):
            fasta_input = open(f, 'r')
            seqs = SeqIO.to_dict(SeqIO.parse(fasta_input, 'fasta'))
            seq_names = list(seqs.keys())
            seq_names.sort()
            for s in seq_names:
                seq_record = seqs[s]
                seq_len = len(str(seq_record.seq))
                seq_id = seq_record.id
                seq_type = determine_seq_type(seq_id, args.plasmid_prefix)
                if seq_type == 'plasmid':
                    num_plasmids += 1
                else:
                    num_chr += 1
                seq_lengths.write(f"{i}\t{seq_id}\t{seq_type}\t{seq_len}\n")
                with open(os.path.join(args.outdir, f"{i}__{seq_type}__{seq_id}.fna"), 'w') as fasta_output:
                    fasta_output.write(seq_record.format('fasta'))
            fasta_input.close()
            seq_counts.write(f"{i}\t{num_chr}\t{num_plasmids}\n")
        else:
            print(f"Warning: {f} does not exist or is not a file. Skipped this file.", file = sys.stderr)
    seq_counts.close()
    seq_lengths.close()
    return


def extract_isolate_name(filename, ext):
    """ Extracts isolate name from [isolate name].[filename extension, fna/fasta] """
    filename = os.path.basename(filename)
    if ext.startswith('.'):
        ext = ext[1 : ]
    n = len(ext)
    if n > 0 and not filename.startswith('.'):
        i = filename[ : len(filename) - n - 1]
    else:
        print("Error: filename extension is not specified or starts with \'.\'.", file = sys.stderr)
        exit(1)
    return i


def determine_seq_type(seq_name, plasmid_prefix):
    """ Determines if a sequence is a chromosome or plasmid based on its sequence ID """
    return 'plasmid' if seq_name.startswith(plasmid_prefix) else 'chromosome'


if __name__ == '__main__':
    main()
