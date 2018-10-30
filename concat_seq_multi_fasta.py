#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
concat_seq_multi_fasta.py - Concatenate all sequences found in a fasta file into a single sequence

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: October 30, 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os.path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature


version = '0.1'


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 concat_seq_multi_fasta.py"')

    parser = argparse.ArgumentParser(prog='concat_seq_multi_fasta.py',
                                     description='Concatenate all sequences found in a fasta file into a single'
                                                 ' sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-f', '--fasta', nargs=1, type=argparse.FileType('r'), required=True,
                                 metavar='/path/to/multi/fasta/file.fasta',
                                 help='Path to the input multi fasta file')

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outfile', type=str, metavar='/path/to/output/concatenated/file.fasta',
                                         help='Path to the directory where the sequences will be stored',
                                         required=False, default='concatenated.fasta')
    parser_optional_general.add_argument('-s', '--spacer', type=int, metavar='10',
                                         help='Number of "N"s to be added between sequences',
                                         required=False, default=100)

    args = parser.parse_args()

    print('\n'
          '===>  RUNNING  concat_seq_multi_fasta.py  <===')

    args.fasta = os.path.abspath(args.fasta[0].name)

    args.outfile = os.path.abspath(args.outfile)
    if not os.path.isdir(os.path.dirname(args.outfile)):
        os.makedirs(os.path.dirname(args.outfile))

    concatenated = ''
    features = {}
    for seq in SeqIO.parse(args.fasta, 'fasta'):
        id_seq = seq.id
        start = len(concatenated)
        end = len(concatenated) + len(seq)
        features[id_seq] = (start, end)
        concatenated += str(seq.seq)
        concatenated += 'N' * args.spacer
    if args.spacer > 0:
        concatenated = concatenated[:-args.spacer]

    concatenated = SeqRecord(Seq(concatenated, generic_dna),
                             id=os.path.splitext(os.path.basename(args.outfile))[0],
                             description='')

    with open(args.outfile, 'wt', newline='\n') as writer:
        _ = SeqIO.write(concatenated, writer, "fasta")

    for id_seq, info in list(features.items()):
        my_start_pos = SeqFeature.ExactPosition(info[0])
        my_end_pos = SeqFeature.ExactPosition(info[1])

        my_feature_location = SeqFeature.FeatureLocation(my_start_pos, my_end_pos)

        my_feature_type = "misc_feature"
        qualifiers = {'label': id_seq, 'note': 'Geneious type: Concatenated sequence'}
        strand = 0

        my_feature = SeqFeature.SeqFeature(location=my_feature_location, type=my_feature_type, qualifiers=qualifiers,
                                           strand=strand)

        concatenated.features.append(my_feature)

    with open('{base}.{new_type}'.format(base=args.outfile, new_type='gb'), 'wt', newline='\n') as writer:
        _ = SeqIO.write(concatenated, writer, 'genbank')


if __name__ == "__main__":
    main()
