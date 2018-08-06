#!/usr/bin/env python2

# -*- coding: utf-8 -*-

"""
get_read_count_nominalLength.py - Get number of reads and reads
nominal length (the expected size of the insert) from ENA IDs

Copyright (C) 2017 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: March 09, 2017

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

import os
import argparse
import time


version = '0.1'


def get_taxon_run_ids(ids_list_seq_from_web_taxon_file):
    list_ids = []
    with open(ids_list_seq_from_web_taxon_file, 'rtU') as reader:
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')
                    list_ids.append(line[0])
    return list_ids


def get_list_ids_from_file(list_ids_file):
    list_ids = []
    with open(list_ids_file, 'rtU') as lines:
        for line in lines:
            line = line.splitlines()[0]
            if len(line) > 0:
                list_ids.append(line)
    return list_ids


def get_read_run_info(ena_id):
    import urllib

    url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=' + ena_id + '&result=read_run'

    read_run_info = None
    try:
        url = urllib.urlopen(url)
        read_run_info = url.read().splitlines()
        if len(read_run_info) <= 1:
            read_run_info = None
    except Exception as error:
        print(error)

    return read_run_info


def get_sequencing_information(read_run_info):
    header_line = read_run_info[0].split('\t')
    info_line = read_run_info[1].split('\t')

    sequencing_information = {'run_accession': None, 'nominal_length': None, 'read_count': None, 'base_count': None}

    for i in range(0, len(header_line)):
        header = header_line[i].lower()
        if header in sequencing_information.keys():
            if len(info_line[i]) > 0:
                sequencing_information[header] = info_line[i]

    return sequencing_information


def get_read_count_nominal_length(ena_id):
    sequencing_information = {'run_accession': None, 'nominal_length': None, 'read_count': None, 'base_count': None}

    read_run_info = get_read_run_info(ena_id)
    if read_run_info is not None:
        sequencing_information = get_sequencing_information(read_run_info)

    return sequencing_information


def main():
    parser = argparse.ArgumentParser(prog='get_read_count_nominalLength.py',
                                     description='Get number of reads and reads nominal length from ENA IDs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-l', '--listIDs', type=argparse.FileType('r'), metavar='/path/to/list_IDs.txt',
                                 help='Path to list containing the IDs to be downloaded (one per line)', required=True)

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the information will be stored',
                                         required=False, default='.')
    parser_optional_general.add_argument('--taxon', action='store_true',
                                         help='Specifies that the listIDs file is from'
                                              ' ReMatCh <https://github.com/B-UMMI/ReMatCh/> seqFromWebTaxon.py')

    args = parser.parse_args()

    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    print('\n' + 'Getting list of IDs')
    if args.taxon:
        list_ids = get_taxon_run_ids(os.path.abspath(args.listIDs.name))
    else:
        list_ids = get_list_ids_from_file((os.path.abspath(args.listIDs.name)))
    print(str(len(list_ids)) + ' IDs to process')

    print('\n' + 'Getting samples information')
    with open(os.path.join(outdir, 'ENA_read_count_nominalLength.' + time.strftime("%Y%m%d-%H%M%S") + '.tab'), 'wt')\
            as writer:
        header = ['run_accession', 'nominal_length', 'read_count', 'base_count']
        writer.write('#' + '\t'.join(header) + '\n')
        counter = 0
        for sample in list_ids:
            sequencing_information = get_read_count_nominal_length(sample)
            if sequencing_information['run_accession'] is not None:
                writer.write('\t'.join([str(sequencing_information[i]) for i in header]) + '\n')
            counter += 1
            if counter % 100 == 0:
                print(str(counter) + ' IDs already processed')


if __name__ == "__main__":
    main()
