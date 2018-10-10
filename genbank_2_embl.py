#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
spyogenes_nga_promoter_variant.py - Determines Streptococcus pyogenes
nga gene promoter variant using fasta sequences

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: September 03, 2018

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
import argparse
import os
import sys

from Bio import SeqIO


version = '0.1'


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 genbank_2_embl.py"')

    parser = argparse.ArgumentParser(prog='genbank_2_embl.py',
                                     description='Converts a genbank file into a embl file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-g', '--genbank', nargs=1, type=argparse.FileType('r'), required=True,
                                 metavar='/path/to/genbank/file.gb',
                                 help='Path to the genbank file')
    parser_required.add_argument('-e', '--embl', nargs=1, type=str,
                                 metavar='/path/to/output/embl/file.embl',
                                 help='Path to the output embl file',
                                 required=True)

    args = parser.parse_args()

    args.genbank = os.path.abspath(args.genbank[0].name)
    args.embl = os.path.abspath(args.embl[0])

    if not os.path.isdir(os.path.dirname(args.embl)):
        os.makedirs(os.path.dirname(args.embl))

    count = SeqIO.convert(args.genbank, 'genbank', args.embl, 'embl')
    print('Converted {} records'.format(count))


if __name__ == "__main__":
    main()
