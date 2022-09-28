"""
Creating a dot plot from each sequence in a multi-FASTA file.

This script is adapted from Ryan Wick's script dotplot.py, which is part of Trycycler (github.com/rrwick/Trycycler).

Useage: dotplot.py -i assembly.fna -o dotplots
Dependencies: Python v3, BioPython, PIL (pip install Pillow)

Copyright (C) 2022 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 28 September 2022; the latest update: 28 September 2022

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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PIL import Image, ImageDraw, ImageFont, ImageColor

# Some hard-coded constants for creating the dotplot. Values between 0 and 1 are relative to full image resolution.
INITIAL_TOP_LEFT_GAP = 0.1
BORDER_GAP = 0.015
BETWEEN_SEQ_GAP = 0.0125
OUTLINE_WIDTH = 0.0015
TEXT_GAP = 0.0025
MAX_FONT_SIZE = 0.025
BACKGROUND_COLOUR = 'white'
SELF_VS_SELF_COLOUR = 'lightgrey'
SELF_VS_OTHER_COLOUR = 'whitesmoke'
TEXT_COLOUR = 'black'
FORWARD_STRAND_DOT_COLOUR = 'mediumblue'
REVERSE_STRAND_DOT_COLOUR = 'firebrick'

def parse_argument():
	parser = ArgumentParser(description = "Creating dot plots from a multi-FASTA file")
	parser.add_argument('-i', '--input', dest = 'input', type = str, required = True, help = "An input FASTA file")
	parser.add_argument('-o', '--outdir', dest = 'outdir', type = str, required = False, default = '.', help = "Output directory (default: the current working directory)")
	return parser.parse_args()


def main():
    args = parse_argument()

    """ Sanity check """
    if not os.path.isdir(args.outdir):
        sys.exit(f"\nError: output directory ({args.outdir}) does not exist")
    if not os.path.exists(args.input):
        sys.exit(f"\nError: input file {args.input} does not exist.")

    """ Process input sequences """
    seqs = load_fasta(args.i)

    return


if __name__ == '__main__':
    main()


def dotplot(args):
    
    seq_names = sorted(seqs.keys())
    image = create_dotplots(seq_names, seqs, args)

    dotplots_filename = args.cluster_dir / 'dotplots.png'
    image.save(dotplots_filename)
    log(f'Saving dotplots to: {dotplots_filename}')
    log()
    return


def create_dotplots(seq_names, seqs, args):
    seq_names = sorted(seq_names)

    # Set some sizes in terms of pixels
    top_left_gap = int(round(INITIAL_TOP_LEFT_GAP * args.res))
    border_gap = max(2, int(round(BORDER_GAP * args.res)))
    between_seq_gap = max(2, int(round(BETWEEN_SEQ_GAP * args.res)))
    text_gap = max(1, int(round(TEXT_GAP * args.res)))
    outline_width = max(1, int(round(OUTLINE_WIDTH * args.res)))
    max_font_size = max(1, int(round(MAX_FONT_SIZE * args.res)))

    # We create an initial image to test the label sizes.
    start_positions, end_positions, bp_per_pixel = \
        get_positions(args, seq_names, seqs, args.kmer, top_left_gap, border_gap, between_seq_gap)
    image = Image.new('RGB', (args.res, args.res), BACKGROUND_COLOUR)
    min_font_size, max_text_height = \
        draw_labels(image, seq_names, start_positions, end_positions, text_gap, outline_width,
                    max_font_size)

    # Now that we know the values for min_font_size and max_text_height, we start over, this time
    # limiting the font size to the minimum (so all text is the same size) and readjusting the
    # top-left gap (so it isn't bigger than necessary).
    top_left_gap = max_text_height + border_gap
    start_positions, end_positions, bp_per_pixel = \
        get_positions(args, seq_names, seqs, args.kmer, top_left_gap, border_gap, between_seq_gap)
    image = Image.new('RGB', (args.res, args.res), BACKGROUND_COLOUR)
    draw_sequence_boxes(image, seq_names, start_positions, end_positions, outline_width)
    draw_labels(image, seq_names, start_positions, end_positions, text_gap, outline_width,
                min_font_size)

    for name_a in seq_names:
        seq_a = seqs[name_a]
        for name_b in seq_names:
            seq_b = seqs[name_b]
            log(f'  {name_a} vs {name_b}')
            draw_dots(image, name_a, name_b, seq_a, seq_b, start_positions, bp_per_pixel, args.kmer)

    # The boxes are drawn once more, this time with no fill. This is to overwrite any dots which
    # leaked into the outline, which would look messy.
    draw_sequence_boxes(image, seq_names, start_positions, end_positions, outline_width, fill=False)

    log()
    return image


def get_positions(args, seq_names, seqs, kmer_size, top_left_gap, bottom_right_gap,
                  between_seq_gap):
    """
    This function returns the image coordinates that start/end each sequence. Since the dot plot is
    symmetrical, there is only one start/end per sequence (used for both x and y coordinates).
    """
    seq_lengths = {n: len(seqs[n]) - kmer_size + 1 for n in seq_names}
    all_gaps = top_left_gap + bottom_right_gap + between_seq_gap * (len(seq_names) - 1)
    pixels_for_sequence = args.res - all_gaps
    total_seq_length = sum(seq_lengths[n] for n in seq_names)
    bp_per_pixel = total_seq_length / pixels_for_sequence

    start_positions, end_positions = {}, {}
    current_pos = top_left_gap
    for name in seq_names:
        start_positions[name] = current_pos
        rect_size = int(round(seq_lengths[name] / bp_per_pixel))
        current_pos += rect_size
        end_positions[name] = current_pos
        current_pos += between_seq_gap

    return start_positions, end_positions, bp_per_pixel


def draw_sequence_boxes(image, seq_names, start_positions, end_positions, outline_width, fill=True):
    """
    This function draws the box for each of the dot plots in the full image.
    """
    draw = ImageDraw.Draw(image)
    for name_a in seq_names:
        start_a = start_positions[name_a] - outline_width
        end_a = end_positions[name_a] + outline_width
        for name_b in seq_names:
            start_b = start_positions[name_b] - outline_width
            end_b = end_positions[name_b] + outline_width
            if fill:
                if name_a == name_b:
                    fill_colour = SELF_VS_SELF_COLOUR
                else:
                    fill_colour = SELF_VS_OTHER_COLOUR
                draw.rectangle([(start_a, start_b), (end_a, end_b)],
                               fill=fill_colour, outline='black', width=outline_width)
            else:
                draw.rectangle([(start_a, start_b), (end_a, end_b)],
                               outline='black', width=outline_width)


def draw_labels(image, seq_names, start_positions, end_positions, text_gap, outline_width,
                font_size):
    draw = ImageDraw.Draw(image)

    min_pos = min(p for p in start_positions.values())
    font_sizes, text_heights = [], []
    for name in seq_names:
        font, text_width, text_height, font_size = \
            get_font(draw, name, font_size, start_positions[name], end_positions[name])
        font_sizes.append(font_size)
        text_heights.append(text_height)

        # Horizontal labels on the top side.
        pos = min_pos - text_height - outline_width - text_gap
        draw.text((start_positions[name], pos), name, font=font, fill=TEXT_COLOUR)

        # Vertical labels on the left side.
        image_2 = Image.new('RGBA', (text_width, text_height), BACKGROUND_COLOUR)
        draw_2 = ImageDraw.Draw(image_2)
        draw_2.text((0, 0), text=name, font=font, fill=TEXT_COLOUR)
        image_2 = image_2.rotate(90, expand=1)
        sx, sy = image_2.size
        image.paste(image_2, (pos, end_positions[name] - sy, pos + sx, end_positions[name]),
                    image_2)
    return min(font_sizes), max(text_heights)


def draw_dots(image, name_a, name_b, seq_a, seq_b, start_positions, bp_per_pixel, kmer_size):
    pixels = image.load()
    a_start_pos = start_positions[name_a]
    b_start_pos = start_positions[name_b]

    forward_colour = ImageColor.getrgb(FORWARD_STRAND_DOT_COLOUR)
    reverse_colour = ImageColor.getrgb(REVERSE_STRAND_DOT_COLOUR)

    a_forward_kmers, a_reverse_kmers = get_all_kmer_positions(kmer_size, seq_a)

    for j in range(len(seq_b) - kmer_size + 1):
        j_pixel = int(round(j / bp_per_pixel)) + b_start_pos
        k = seq_b[j:j+kmer_size]
        if k in a_reverse_kmers:
            for i in a_reverse_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + a_start_pos
                draw_dot(pixels, i_pixel, j_pixel, reverse_colour)

        if k in a_forward_kmers:
            for i in a_forward_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + a_start_pos
                draw_dot(pixels, i_pixel, j_pixel, forward_colour)


def draw_dot(pixels, i, j, colour):
    try:
        pixels[i, j] = colour
        pixels[i+1, j] = colour
        pixels[i-1, j] = colour
        pixels[i, j+1] = colour
        pixels[i, j-1] = colour
    except IndexError:
        pass


def get_all_kmer_positions(kmer_size, seq):
    forward_kmers, reverse_kmers = collections.defaultdict(list), collections.defaultdict(list)
    rev_comp_seq = reverse_complement(seq)
    seq_len = len(seq) - kmer_size + 1
    for i in range(seq_len):
        k = seq[i:i+kmer_size]
        forward_kmers[k].append(i)
        k = rev_comp_seq[i:i+kmer_size]
        reverse_kmers[k].append(seq_len - i - 1)
    assert len(forward_kmers) < len(seq)
    assert len(reverse_kmers) < len(seq)
    return forward_kmers, reverse_kmers


def load_fasta(fasta_filename, include_full_header=False):
    fasta_seqs = []
    with open(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    if include_full_header:
                        fasta_seqs.append((name.split()[0], name, ''.join(sequence)))
                    else:
                        fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            if include_full_header:
                fasta_seqs.append((name.split()[0], name, ''.join(sequence)))
            else:
                fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


def load_contig_sequences(cluster_dir):
    filenames = get_contigs_from_cluster_dir(cluster_dir)
    contig_seqs, fasta_names = {}, {}
    for f in filenames:
        seqs = load_fasta(f)
        assert len(seqs) == 1
        name, seq = seqs[0]
        contig_seqs[name] = seq
        fasta_names[name] = f
    return contig_seqs, fasta_names


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def get_font(draw, label, font_size, start_position, end_position):
    font, is_default_font = load_font(font_size)

    # If we had to resort to the default font, then we can't size it.
    if is_default_font:
        text_width, text_height = draw.textsize(label, font=font)
        return font, text_width, text_height, font_size

    # If we have a size-able font, then we adjust the size down until the text fits in the
    # available space.
    available_width = end_position - start_position
    text_width, text_height = draw.textsize(label, font=font)
    while text_width > available_width:
        font_size -= 1
        font, _ = load_font(font_size)
        text_width, text_height = draw.textsize(label, font=font)
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
    print(message, file=sys.stderr, flush = True, end = end)
