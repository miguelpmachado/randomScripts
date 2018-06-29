#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
create_emm_type_db_from_cdc.py - Creates the emm type DB from CDC to be used with emm_pipeline-streplab

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: June 29, 2018

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
import argparse
import os
import urllib.request
from Bio import SeqIO


version = '0.1'


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 create_emm_type_db_from_cdc.py"')

    parser = argparse.ArgumentParser(prog='python3 create_emm_type_db_from_cdc.py',
                                     description="Creates the emm type DB from CDC to be used with"
                                                 " emm_pipeline-streplab",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_optional = parser.add_argument_group('Facultative options')
    parser_optional.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                 help='Path to the directory where the data will be stored',
                                 required=False, default='.')
    parser_optional.add_argument('-f', '--fasta', type=str, metavar='ftp://path/to/file/emm_type.fasta',
                                 help='Path to the fasta file that contains emm type sequences that will be downloaded'
                                      ' and parsed',
                                 required=False,
                                 default='ftp://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/trimmed.tfa')

    args = parser.parse_args()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    emm_fasta = os.path.join(args.outdir, 'temp_emm.fasta')
    print('Downloading emm types fasta file')
    urllib.request.urlretrieve(args.fasta, emm_fasta)
    print('Creating emm DB files for emm_pipeline-streplab')
    with open(emm_fasta, mode='rt') as reader_fasta, \
            open(os.path.join(args.outdir, 'emm.fasta'), mode='wt') as writer_fasta, \
            open(os.path.join(args.outdir, 'emm_reference.txt'), mode='wt') as writer_srst2_tab:
        writer_srst2_tab.write('\t'.join(['ST', 'emm']) + '\n')
        for index, record in enumerate(SeqIO.parse(reader_fasta, 'fasta')):
            record.id = record.id.lower()
            writer_srst2_tab.write('\t'.join([record.id, record.id]) + '\n')
            record.id = 'emm-{id}'.format(id=record.id)
            record.description = ''
            writer_fasta.write(record.format('fasta'))
    if os.path.isfile(emm_fasta):
        os.remove(emm_fasta)


if __name__ == "__main__":
    main()
