#!/usr/bin/env python3

"""
Estimate the average and individual depth of contigs in a genome assembly.

Usage:
    python estimate_read_depths.py --help
    python estimate_read_depths.py {illuminaPE, nanopore} {--help or parameters}

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 26 Oct 2024; the latest update: 27 Oct 2024.
"""

import os
import sys
import shutil
import argparse
import subprocess
from datetime import datetime
from collections import namedtuple

# Define data classes ###############
Reads = namedtuple('Reads', ['r1', 'r2'])
Software = namedtuple('Software', ['minimap2', 'samtools', 'mosdepth'])

def parse_arguments():  # https://stackoverflow.com/questions/7066826/in-python-how-to-get-subparsers-to-read-in-parent-parsers-argument
    # Shared arguments
    parser_common = argparse.ArgumentParser(description = "Estimate read depths for genome assemblies.", add_help=False)
    parser_common.add_argument('--assembly', dest = 'assembly', type = str, required = True, help = "Input FASTA file of the genome assembly to be analysed")
    parser_common.add_argument('--prefix', dest = 'prefix', type = str, required = True, help = "Prefix for output files")
    parser_common.add_argument('--output_dir', dest = 'output_dir', type = str, required = False, default = '.', help = "Directory to save output files (default: current working directory)")
    parser_common.add_argument('--minimap2_dir', dest = 'minimap2_dir', type = str, required = False, default = '', help = "Directory of minimap2 (default: none)")
    parser_common.add_argument('--samtools_dir', dest = 'samtools_dir', type = str, required = False, default = '', help = "Directory of samtools (default: none)")
    parser_common.add_argument('--mosdepth_dir', dest = 'mosdepth_dir', type = str, required = False, default = '', help = "Directory of mosdepth (default: none)")
    parser_common.add_argument('--threads', dest = 'threads', type = str, required = False, default = '1', help = "Number of threads (default: 1)")
    parser_common.add_argument('--window_size', dest = 'window_size', type = str, required = False, default = '200', help = "Window size (bp) for mosdepth (default: 200)")
    
    # Subcommands
    parser_parent = argparse.ArgumentParser(description = "Parental argument parser")
    subparsers = parser_parent.add_subparsers(dest = 'subcommand', required = True, help = "Subcommand (nanopore or illuminaPE)")
    
    # Subcommand: nanopore
    parser_nanopore = subparsers.add_parser(name = 'nanopore', parents = [parser_common], help = "Estimate read depth for Nanopore reads")
    parser_nanopore.add_argument('--r', dest = 'r', help = "Nanopore read file")

    # Subcommand: illuminaPE
    parser_illumina = subparsers.add_parser(name = 'illuminaPE', parents = [parser_common], help = "Estimate read depth for Illumina paired-end reads")
    parser_illumina.add_argument('--r1', dest = 'r1', help = "Illumina read file 1 for forward reads")
    parser_illumina.add_argument('--r2', dest = 'r2', help = "Illumina read file 2 for reverse reads")

    return parser_parent.parse_args()


def main():
    args = parse_arguments()
    log = initiate_outputs(args.output_dir, '_'.join([args.prefix, args.subcommand]))
    software = check_dependencies(args, log)
    print_run_info(software, args, log)
    if args.subcommand == "nanopore":
        reads = Reads(r1 = args.r, r2 = None)
    else:
        reads = Reads(r1 = args.r1, r2 = args.r2)
    estimate_depths(reads, args, software, log)

    return


def initiate_outputs(outdir, prefix):
    """ Check output directory and create log files """
    new_dir = False
    if not os.path.exists(outdir):
        os.makedirs(outdir)  # Same as `mkdir -p [outdir]` in Shell
        new_dir = True
    log = os.path.join(outdir, prefix + '_estimate_read_depths.log')
    if new_dir:
        with open(log, 'w') as log_f:
            log_f.write(f"Warning: created output directory {outdir} since was not found.\n")
    else:
        open(log, 'w').close()  # Create an empty log file
    
    return log


def add_log(message, log):
    with open(log, 'a') as log_file:
        print(message + '\n', file = log_file)
    return


def find_software(n, p, log):
    """ Search for software named n in directory p or in system PATH """
    found = False
    if p != '':  # If a specific path is provided
        software_path = os.path.join(p, n)
        if os.path.isfile(software_path) and os.access(software_path, os.X_OK):
            found = True
        else:
            add_log(f"Error: {n} was not found in the specified directory {p}", log)
    if not found or p == '':
        software_path = shutil.which(n)
        if software_path is None:
            add_log(f"Error: {n} was not found in the system PATH", log)
            sys.exit(1)
    
    return software_path


def check_dependencies(args, log):
    minimap2_path = find_software('minimap2', args.minimap2_dir, log)
    samtools_path = find_software('samtools', args.samtools_dir, log)
    mosdepth_path = find_software('mosdepth', args.mosdepth_dir, log)
    software_paths = Software(minimap2 = minimap2_path, samtools = samtools_path, mosdepth = mosdepth_path)
    
    return software_paths


def print_run_info(software, args, log):
    lines = []
    run_type = args.subcommand
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines.append(f"Genome name: {args.prefix}")
    lines.append(f"Run started at: {current_time}")
    lines.append(f"Input type: {run_type}")
    if run_type == 'nanopore':
        lines.append(f"Input reads: {args.r}")
    else:
        lines.append(f"Input reads: {args.r1}; {args.r2}")
    lines.append(f"Output directory: {args.output_dir}")
    lines.append(f"Number of threads: {args.threads}")
    lines.append(f"Mosdepth: {software.mosdepth}")
    lines.append(f"Minimap2: {software.minimap2}")
    lines.append(f"Samtools: {software.samtools}")
    add_log('\n'.join(lines) + '\n', log)
    
    return


def execute_command(command, command_name, log_file):
    try:
        result = subprocess.run(command, shell = True, check = True, text = True, capture_output = True)
        add_log(f"{command_name} command executed successfully: {command}", log_file)
        add_log(f"{command_name} stdout output:\n{result.stdout}", log_file)
    except subprocess.CalledProcessError as err:  # Log error if command fails
        add_log(f"Error:\n{err.stderr}", log_file)
        sys.exit(1)
    
    return


def estimate_depths(reads, args, software, log):
    """ Run mosdepth to estimate read depths """
    output_bam = os.path.join(args.output_dir, f"{args.prefix}_{args.subcommand}_sorted.bam")
    if args.subcommand == 'nanopore':
        command = f"{software.minimap2} -x map-ont --MD -a -t {args.threads} {args.assembly} {reads.r1} | " + \
                  f"{software.samtools} view -S -b -u -@ {args.threads} - | " + \
                  f"{software.samtools} sort -@ {args.threads} -o {output_bam} -"
    else:  # github.com/lh3/minimap2/blob/master/cookbook.md#map-sr
        command = f"{software.minimap2} -x sr --MD -a -t {args.threads} {args.assembly} {reads.r1} {reads.r2} | " + \
                  f"{software.samtools} view -S -b -u -@ {args.threads} - | " + \
                  f"{software.samtools} sort -@ {args.threads} -o {output_bam} -"
    
    # Read mapping
    execute_command(command, 'Minimap & samtools', log)

    # Depth estimation
    if os.path.exists(output_bam):
        # Create an index file
        command = f"{software.samtools} index -@ {args.threads} {output_bam}"
        execute_command(command, 'samtools index', log)
        # Run mosdepth
        mosdepth_output_prefix = os.path.join(args.output_dir, '_'.join([args.prefix, args.subcommand]))
        command = f"{software.mosdepth} --no-per-base --fast-mode --by {args.window_size} --threads {args.threads} " + \
                  f"{mosdepth_output_prefix} {output_bam}"
        execute_command(command, 'Mosdepth', log)
    else:
        add_log(f"Error: output file {output_bam} does not exist. Exit.", log)
        sys.exit(1)

    # Filter the output depth summary to remove rows containing the keyword '_region'
    depth_summary = mosdepth_output_prefix + '.mosdepth.summary.txt'
    filtered_summary = mosdepth_output_prefix + '_mosdepth_depth_per_contig.tsv'
    command = f"awk 'NR==1 || $1 !~ /_region/' {depth_summary} > {filtered_summary}"
    execute_command(command, 'awk', log)
    
    return


if __name__ == '__main__':
    main()
