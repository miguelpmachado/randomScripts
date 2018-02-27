#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
sample_fastq.py - Randomly sampling reads from fastq files
<https://github.com/miguelpmachado/randomScripts/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: January 03, 2018

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

import random
import argparse
import gzip
import bz2
import zipfile
import sys
import os
import shlex
import subprocess
from threading import Timer

version = '0.1'


def kill_subprocess_popen(subprocess_popen, command):
    print('Command run out of time: ' + str(command))
    subprocess_popen.kill()


def run_command_popen_communicate(command, shell_true, timeout_sec_none, print_comand_true):
    run_successfully = False
    if not isinstance(command, str):
        command = ' '.join(command)

    if print_comand_true:
        print('Running: ' + command)

    if shell_true:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        command = shlex.split(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    not_killed_by_timer = True
    if timeout_sec_none is None:
        stdout, stderr = proc.communicate()
    else:
        timer = Timer(timeout_sec_none, kill_subprocess_popen, args=(proc, command,))
        timer.start()
        stdout, stderr = proc.communicate()
        timer.cancel()
        not_killed_by_timer = timer.isAlive()

    if proc.returncode == 0:
        run_successfully = True
    else:
        if not print_comand_true and not_killed_by_timer:
            print('Running: ' + str(command))
        if len(stdout) > 0:
            print('STDOUT')
            print(stdout.decode("utf-8"))
        if len(stderr) > 0:
            print('STDERR')
            print(stderr.decode("utf-8"))
    return run_successfully, stdout.decode("utf-8"), stderr.decode("utf-8")


def required_length(tuple_length_options, argument_name):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values) not in tuple_length_options:
                msg = 'Option {argument_name} requires one of the following number of' \
                      ' arguments: {tuple_length_options}'.format(argument_name=argument_name,
                                                                  tuple_length_options=tuple_length_options)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def guess_file_compression(file_path):
    magic_dict = {
        b'\x1f\x8b\x08': 'gz',
        b'\x42\x5a\x68': 'bz2',
        b'\x50\x4b\x03\x04': 'zip'
    }

    max_len = max(len(x) for x in magic_dict)

    with open(file_path, 'rb') as f:
        file_start = f.read(max_len)

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return file_type

    return None


def open_fastq_files(fastq):
    copen = {
        'gz': gzip.open,
        'bz2': bz2.open,
        'zip': zipfile.ZipFile
    }

    file_objects = []
    for fastq_files in fastq:
        ftype = guess_file_compression(fastq_files)

        if ftype is not None:
            file_objects.append(copen[ftype](fastq_files, mode='rt'))
        else:
            file_objects.append(open(fastq_files))

    return file_objects


def count_reads(fastq):
    read_fastq_files = open_fastq_files(fastq)

    total_records = []
    for fastq_file in read_fastq_files:
        total_records.append(int(sum([1 for _ in fastq_file]) / 4))
        fastq_file.close()

    if len(total_records) > 1 and len(set(total_records)) > 1:
        sys.exit('\n'.join(['Different number of reads found in the provided fastq files'] +
                           ['{fastq_file}: {n_reads} reads'.format(fastq_file=fastq_file, n_reads=total_records[x])
                            for x, fastq_file in enumerate(fastq)]))

    return total_records[0]


def count_reads_gzip(fastq):
    total_records = []
    for fastq_file in fastq:
        command = ['gunzip', '--keep', '--stdout', fastq_file, '|', 'awk', "\'{if(NR%4==1) print}\'", '|', 'wc', '-l']
        run_successfully, stdout, stderr = run_command_popen_communicate(command, True, None, True)
        if run_successfully:
            total_records.append(int(stdout.split()[0]))

    if len(total_records) == 0:
        sys.exit('It was not possible to determine the number of reads')
    if len(total_records) > 1 and len(set(total_records)) > 1:
        sys.exit('\n'.join(['Different number of reads found in the provided fastq files'] +
                           ['{fastq_file}: {n_reads} reads'.format(fastq_file=fastq_file, n_reads=total_records[x])
                            for x, fastq_file in enumerate(fastq)]))

    return total_records[0]


def fastq_sampling(fastq, out_fastq, total_records, number):
    read_fastq_files = open_fastq_files(fastq)
    write_fastq_files = [gzip.open(fastq_file, 'wt') for fastq_file in out_fastq]

    output_sequence_sets = set(random.sample(range(total_records + 1), number))

    for x, reader in enumerate(read_fastq_files):
        print('-' * 15)
        writer = write_fastq_files[x]
        record_number = 0
        for line in reader:
            line1 = line.strip()
            line2 = reader.readline().strip()
            line3 = reader.readline().strip()
            line4 = reader.readline().strip()
            if record_number in output_sequence_sets:
                writer.write(line1 + '\n')
                writer.write(line2 + '\n')
                writer.write(line3 + '\n')
                writer.write(line4 + '\n')
            record_number += 1
            if record_number % 500000 == 0:
                print('{file}: {percentage}%'.format(file=fastq[x], percentage=round((record_number /
                                                                                      float(total_records)) * 100,
                                                                                     2)))
        reader.close()
        writer.close()
        print('{file}: 100%'.format(file=fastq[x]))


def main():
    if sys.version_info.major != 3:
        sys.exit('\n' + '"sample_fastq.py" requires Python 3' + '\n' + 'Try running with "python3 sample_fastq.py"')

    parser = argparse.ArgumentParser(prog='python3 sample_fastq.py',
                                     description='Randomly sampling reads from fastq files',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-f', '--fastq', nargs='+', action=required_length((1, 2), '--fastq'),
                                 type=argparse.FileType('r'), metavar='/path/to/input/file.fq.gz',
                                 help='Path to single OR paired-end fastq files. If two files are passed, they will be'
                                      ' assumed as being the paired fastq files. Can be used compressed (gzip, bzip2,'
                                      ' zip) or decompressed files.', required=True)
    parser_required.add_argument('-o', '--outFastq', nargs='+', action=required_length((1, 2), '--outFastq'), type=str,
                                 metavar='/path/to/output/file_downsampled.fq.gz',
                                 help='Path to the downsampled fastq file gzip compressed. If paired-end fastq files'
                                      ' are provided, write the path to both downsampled fastq files in the same order'
                                      ' as the input fastq files.', required=True)

    downsampling_options = parser.add_mutually_exclusive_group()
    downsampling_options.add_argument('--fraction', type=float, metavar='0.1',
                                      help='Fraction of reads to sample [0, 1]', required=False)
    downsampling_options.add_argument('--number', type=int, metavar='1000', help='Number of reads to sample',
                                      required=False)

    args = parser.parse_args()

    if len(args.fastq) != len(args.outFastq):
        parser.error('Provide the same number of files with --fastq and --outFastq')

    if args.fraction is None and args.number is None:
        parser.error('Provide one of the following options: --fraction or --number')

    args.outFastq = [os.path.abspath(outFastq) for outFastq in args.outFastq]

    errors_found = []
    for outFastq in args.outFastq:
        if os.path.isfile(outFastq):
            errors_found.append('{} file already exist'.format(outFastq))
        else:
            if not os.path.isdir(os.path.dirname(outFastq)):
                errors_found.append('{} directory does not exist'.format(os.path.dirname(outFastq)))
    if len(errors_found) > 0:
        sys.exit('\n'.join(['Found some problems with output files'] + list(set(errors_found))))

    args.fastq = [os.path.abspath(fastq.name) for fastq in args.fastq]

    print('\n' + 'Counting reads...')
    total_records = count_reads(args.fastq)
    # total_records = count_reads_gzip(args.fastq)
    print('{} reads found'.format(total_records))

    if args.fraction is not None:
        args.number = int(total_records * args.fraction)
        print('{} reads will be sampled corresponding to {}%'.format(args.number, args.fraction * 100))

    if args.number >= total_records:
        sys.exit('The number of reads to subsample are equal or higher than the number of reads found')

    fastq_sampling(args.fastq, args.outFastq, total_records, args.number)

    print('-' * 15)
    print('DONE')


if __name__ == "__main__":
    main()
