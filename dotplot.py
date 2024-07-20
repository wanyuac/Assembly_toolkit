"""
Creating a dot plot from each sequence in a multi-FASTA file.

This script is adapted from Ryan Wick's script dotplot.py, which is part of Trycycler (github.com/rrwick/Trycycler).

Useage: dotplot.py -i assembly.fna -o dotplots
Dependencies: Python v3, BioPython, PIL (pip install Pillow)

Copyright (C) 2022 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 28 September 2022; the latest update: 20 July 2024

Appendix: Copyright information of Ryan Wick's dotplot.py
'Copyright 2020 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Trycycler

This file is part of Trycycler. Trycycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Trycycler is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Trycycler.
If not, see <http://www.gnu.org/licenses/>.'
"""

import os
import sys
import collections
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PIL import Image, ImageDraw, ImageFont, ImageColor

# Some hard-coded constants for creating the dotplot. Values between 0 and 1 are relative to full image resolution.
INITIAL_TOP_LEFT_GAP = 0.1
BORDER_GAP = 0.015
OUTLINE_WIDTH = 0.0015
TEXT_GAP = 0.005
MAX_FONT_SIZE = 0.020
BACKGROUND_COLOUR = 'white'
FILL_COLOUR = 'lightgrey'
TEXT_COLOUR = 'black'
FORWARD_STRAND_DOT_COLOUR = 'mediumblue'
REVERSE_STRAND_DOT_COLOUR = 'firebrick'
REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def parse_argument():
    parser = ArgumentParser(description = "Creating dot plots from a multi-FASTA file")
    parser.add_argument('-i', '--input', dest = 'input', type = str, required = True, help = "An input FASTA file")
    parser.add_argument('-o', '--outdir', dest = 'outdir', type = str, required = False, default = '.', help = "Output directory (default: the current working directory)")
    parser.add_argument('-s', '--sample', dest = 'sample', type = str, required = False, default = None, help = "Sample name to be used as the prefix of output filenames (default: None)")
    parser.add_argument('-p', '--pixels', dest = 'pixels', type = int, required = False, default = 960, help = "Size of output figures (default: 960 pixels)")
    parser.add_argument('-k', '--kmer', dest = 'kmer', type = int, required = False, default = 32, help = "K-mer size for creating dot plots (default: 32 bp)")  # Based on the default k-mer size in https://github.com/rrwick/Trycycler/blob/main/trycycler/__main__.py
    parser.add_argument('-n', '--no_plot', dest = 'no_plot', action = 'store_true', help = "A flag for producing only a summary of contig lengths without dot plots (default: off)")
    return parser.parse_args()


def main():
    args = parse_argument()

    """ Sanity check """
    if not os.path.isdir(args.outdir):
        print(f"\nWarning: output directory ({args.outdir}) does not exist. An output directory will be created.")
        os.mkdir(args.outdir)
    if not os.path.exists(args.input):
        sys.exit(f"\nError: input file {args.input} does not exist.")

    """ Set some sizes in terms of pixels """
    image_sizes = {"top_left_gap" : int(round(INITIAL_TOP_LEFT_GAP * args.pixels)),\
                   "border_gap" : max(2, int(round(BORDER_GAP * args.pixels))),\
                   "text_gap" : max(1, int(round(TEXT_GAP * args.pixels))),\
                   "outline_width" : max(1, int(round(OUTLINE_WIDTH * args.pixels))),\
                   "max_font_size" : max(1, int(round(MAX_FONT_SIZE * args.pixels)))}

    """ Process input sequences and create dotplots """
    log('\t'.join(['Sample', 'Sequence', 'Length_bp']))  # Print the header line in the log file
    seqs = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    for i, s in seqs.items():
        log('\t'.join([args.sample, i, str(len(s.seq))]))
        if not args.no_plot:
            draw_dotplot(i, s.seq, args, image_sizes)  # One dot plot per sequence (contig)
    return


def draw_dotplot(seq_id, seq, args, image_sizes):
    """ Making a dot plot for a sequence (contig) """

    # We create an initial image to test the label sizes.
    start_position, end_position, bp_per_pixel = get_positions(seq_id, seq, args.pixels, args.kmer, image_sizes["top_left_gap"], image_sizes["border_gap"])
    image = Image.new('RGB', (args.pixels, args.pixels), BACKGROUND_COLOUR)  # Create a blank image
    min_font_size, max_text_height = draw_labels(image, seq_id, start_position, end_position, image_sizes["text_gap"], image_sizes["outline_width"], image_sizes["max_font_size"])

    # Now that we know the values for min_font_size and max_text_height, we start over, this time
    # limiting the font size to the minimum (so all text is the same size) and readjusting the
    # top-left gap (so it isn't bigger than necessary).
    top_left_gap = max_text_height + image_sizes["border_gap"]
    start_position, end_position, bp_per_pixel = get_positions(seq_id, seq, args.pixels, args.kmer, top_left_gap, image_sizes["border_gap"])
    image = Image.new('RGB', (args.pixels, args.pixels), BACKGROUND_COLOUR)  # Recreate the blank image
    draw_labels(image, seq_id, start_position, end_position, image_sizes["text_gap"], image_sizes["outline_width"], min_font_size)
    draw_dots(image, seq_id, seq, start_position, bp_per_pixel, args.kmer)
    draw_sequence_box(image, seq_id, start_position, end_position, image_sizes["outline_width"], False)  # This step covers all dots that leak into the border line.
    if args.sample != None:
        image.save(os.path.join(args.outdir, f'{args.sample}__{seq_id}__{args.kmer}-mer.png'))
    else:
        image.save(os.path.join(args.outdir, f'{seq_id}__{args.kmer}-mer.png'))
    return


def get_positions(seq_id, seq, pixels, kmer_size, top_left_gap, bottom_right_gap):
    """
    This function returns the image coordinates that start/end each sequence. Since the dot plot is
    symmetrical, there is only one start/end per sequence (used for both x and y coordinates).
    """
    seq_length = len(seq) - kmer_size + 1  # Range of dot positions in the X axis
    all_gaps = top_left_gap + bottom_right_gap
    width = pixels - all_gaps
    bp_per_pixel = seq_length / width  # A floating number
    start_position = top_left_gap
    end_position = start_position + width
    return start_position, end_position, bp_per_pixel


def draw_labels(image, seq_id, start_position, end_position, text_gap, outline_width, font_size):
    draw = ImageDraw.Draw(image)
    font, text_width, text_height, font_size = get_font(draw, seq_id, font_size, start_position, end_position)

    # Horizontal labels on the top side.
    pos = start_position - text_height - outline_width - text_gap
    draw.text((start_position, pos), seq_id, font = font, fill = TEXT_COLOUR)

    # Vertical labels on the left side.
    image_2 = Image.new('RGBA', (text_width, text_height), BACKGROUND_COLOUR)
    draw_2 = ImageDraw.Draw(image_2)
    draw_2.text((0, 0), text = seq_id, font = font, fill = TEXT_COLOUR)
    image_2 = image_2.rotate(90, expand = 1)
    sx, sy = image_2.size
    image.paste(image_2, (pos, end_position - sy, pos + sx, end_position), image_2)
    return font_size, text_height


def draw_sequence_box(image, seq_id, start_position, end_position, outline_width, fill = True):
    """
    This function draws the box for a dot plot in the full image.
    """
    draw = ImageDraw.Draw(image)
    start = start_position - outline_width
    end = end_position + outline_width
    if fill:
        draw.rectangle([(start, start), (end, end)], fill = FILL_COLOUR, outline = 'black', width = outline_width)
    else:
        draw.rectangle([(start, start), (end, end)], outline = 'black', width = outline_width)
    return


def draw_dots(image, seq_id, seq, start_position, bp_per_pixel, kmer_size):
    pixels = image.load()
    forward_colour = ImageColor.getrgb(FORWARD_STRAND_DOT_COLOUR)
    reverse_colour = ImageColor.getrgb(REVERSE_STRAND_DOT_COLOUR)
    forward_kmers, reverse_kmers = get_all_kmer_positions(kmer_size, seq)

    for j in range(len(seq) - kmer_size + 1):
        j_pixel = int(round(j / bp_per_pixel)) + start_position
        k = seq[j : (j + kmer_size)]
        if k in reverse_kmers:
            for i in reverse_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + start_position
                draw_dot(pixels, i_pixel, j_pixel, reverse_colour)

        if k in forward_kmers:
            for i in forward_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + start_position
                draw_dot(pixels, i_pixel, j_pixel, forward_colour)
    return


def get_all_kmer_positions(kmer_size, seq):
    forward_kmers, reverse_kmers = collections.defaultdict(list), collections.defaultdict(list)  # Dictionaries of lists
    rev_comp_seq = reverse_complement(seq)
    seq_len = len(seq) - kmer_size + 1
    for i in range(seq_len):
        k = seq[i : (i + kmer_size)]  # The sequence of the i-th k-mer
        forward_kmers[k].append(i)
        k = rev_comp_seq[i : i + kmer_size]
        reverse_kmers[k].append(seq_len - i - 1)
    assert len(forward_kmers) < len(seq)
    assert len(reverse_kmers) < len(seq)
    return forward_kmers, reverse_kmers


def draw_dot(pixels, i, j, colour):
    try:
        pixels[i, j] = colour
        pixels[i + 1, j] = colour
        pixels[i - 1, j] = colour
        pixels[i, j + 1] = colour
        pixels[i, j - 1] = colour
    except IndexError:
        pass


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def get_font(draw, label, font_size, start_position, end_position):
    font, is_default_font = load_font(font_size)

    # If we had to resort to the default font, then we can't size it.
    if is_default_font:
        left, top, right, bottom = font.getbbox(label)  # Updated for new Pillow package; depreciated code: "text_width, text_height = draw.textsize(label, font=font)"
        text_width, text_height = right - left, bottom - top
        return font, text_width, text_height, font_size

    # If we have a size-able font, then we adjust the size down until the text fits in the
    # available space.
    available_width = end_position - start_position
    left, top, right, bottom = font.getbbox(label)
    text_width, text_height = right - left, bottom - top  # Depreciated: text_width, text_height = draw.textsize(label, font = font)
    while text_width > available_width:
        font_size -= 1
        font, _ = load_font(font_size)
        left, top, right, bottom = font.getbbox(label)
        text_width, text_height = right - left, bottom - top  # Depreciated: text_width, text_height = draw.textsize(label, font = font)
    return font, text_width, text_height, font_size


def load_font(font_size):
    """
    This function loads a font, but it has to try a few different ones because different platforms
    have different fonts.
    """
    try:
        return ImageFont.truetype('DejaVuSans.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('OpenSans-Regular.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('Arial.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('LiberationSans-Regular.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('NimbusSans-Regular.otf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('Verdana.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('Lato-Regular.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('FreeSans.ttf', font_size), False
    except OSError:
        pass
    return ImageFont.load_default(), True


def log(message = '', end = '\n'):
    print(message, file = sys.stdout, flush = True, end = end)


if __name__ == '__main__':
    main()
