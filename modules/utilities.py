#!/usr/bin/env python

"""
Modules for package Assembly_toolkit.

Copyright (C) 2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 20 February 2024; the latest update: 11 November 2024.
"""

import os
import sys


def write_seq(seq_name, seq_descr, seq, fasta_handle):
    """ Writes a sequence in the FASTA format """
    fasta_handle.write(f">{seq_name} {seq_descr}\n")
    fasta_handle.write(seq + '\n')
    return


def check_dir(p):
    """ Creates a directory if p is not found in the system """
    dir_exists = os.path.isdir(p)
    if dir_exists:
        print(f"Directory {p} exists." )
    elif not os.path.isfile(p):
        print("Create directory " + p)
        os.makedirs(p)
    else:
        print(f"Error: {p} exists as a file, so no directory could be created with this path.", file = sys.stderr)
        sys.exit(1)
    return
