#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
concatenate_fastq_files.py

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: February 02, 2018

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
import sys
import argparse
import time
import multiprocessing
import pickle
import functools
import traceback
import gzip
import bz2
import zipfile

# import subprocess
# from threading import Timer
# import shlex


version = '0.2'


def check_create_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


def start_logger(workdir):
    time_str = time.strftime("%Y%m%d-%H%M%S")
    sys.stdout = Logger(workdir, time_str)
    logfile = sys.stdout.get_log_file()
    return logfile, time_str


class Logger(object):
    def __init__(self, out_directory, time_str):
        self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
        self.terminal = sys.stdout
        self.log = open(self.logfile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        pass

    def get_log_file(self):
        return self.logfile


def file_matching_extension(file_name, fastq_ends):
    matching = False
    for fastq_end in fastq_ends:
        if len(rchop(os.path.basename(file_name), fastq_end)) < len(os.path.basename(file_name)):
            matching = True
            break

    return matching


def read_pair_0_1_2(file_name, fastq_ends):
    orientation_0_1_2 = 0

    if file_name.endswith(fastq_ends[0]):
        orientation_0_1_2 = 1
    elif file_name.endswith(fastq_ends[1]):
        orientation_0_1_2 = 2
    else:
        print('WARNING: no paired information was found for ' + file_name + ' file!')

    return orientation_0_1_2


def gather_files(indir, fastq_ends):
    fastq_files = []
    directories = [d for d in os.listdir(indir) if not d.startswith('.') and os.path.isdir(os.path.join(indir, d, ''))]
    fastq_in_indir = False
    if len(directories) == 0:
        directories.append(os.path.basename(indir))
        fastq_in_indir = True
    for directory in directories:
        if fastq_in_indir:
            directory_path = str(indir)
        else:
            directory_path = os.path.join(indir, directory, '')

        files = [f for f in os.listdir(directory_path) if
                 not f.startswith('.') and os.path.isfile(os.path.join(directory_path, f))]
        for file_found in files:
            if file_matching_extension(file_found, fastq_ends):
                file_path = os.path.join(directory_path, file_found)
                fastq_files.append(file_path)

    return fastq_files


def rchop(string, ending):
    """
    This function removes the ending of a string

    Parameters
    ----------
    string : str
        String from which the ending part will be removed, e.g. 'sample_1.fastq.gz'
    ending : str
        Ending part of the string to be removed, e.g. '_1.fastq.gz'

    Returns
    -------
    string : str
        The original string without the ending part, e.g. 'sample'
    """
    if string.endswith(ending):
        string = string[:-len(ending)]
    return string


def organize_sample_files(list_files, sample_name_delimiter, fastq_ends):
    sample_files = {}

    for fastq_file in list_files:
        fastq_file_name = os.path.basename(fastq_file)

        sample_fastq_file = None
        for fastq_end in fastq_ends:
            sample_fastq_file = rchop(fastq_file_name, fastq_end)
            if len(sample_fastq_file) < len(os.path.basename(fastq_file)):
                break
            else:
                sample_fastq_file = None

        if sample_fastq_file is not None:
            sample = sample_fastq_file.split(sample_name_delimiter)[0]
            file_size = os.path.getsize(fastq_file)

            if sample not in sample_files:
                sample_files[sample] = {}

            if fastq_file_name not in sample_files[sample]:
                sample_files[sample][fastq_file_name] = [{'file_size': file_size, 'file_path': fastq_file}]
            else:
                file_different_size = True
                for file_properties in sample_files[sample][fastq_file_name]:
                    if file_size == file_properties['file_size']:
                        file_different_size = False
                        break

                if file_different_size:
                    sample_files[sample][fastq_file_name].append({'file_size': file_size, 'file_path': fastq_file})

    return sample_files


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


def concatenate_compressed_files(compressed_files, combined_gz):
    """
    Concatenate compressed files into a single file

    Parameters
    ----------
    compressed_files : :obj:`list` of :obj:`str`
        List of files (path) that will be concatenated
    combined_gz : str
        Path to the combined gziped file that will be created

    Returns
    -------
    run_successfully : bool
        Reports if the concatenation run successfully
    """

    # Check if all the files have the same compression type
    files_type = {}
    for file_compressed in compressed_files:
        files_type[file_compressed] = guess_file_compression(file_compressed)
    types = set(files_type.values())

    run_successfully = False
    if len(types) == 1:
        if os.path.isfile(combined_gz):
            os.remove(combined_gz)

        with gzip.open(combined_gz, 'at') as writer:
            for file_compressed in compressed_files:
                reader = open_fastq_files([file_compressed])[0]
                for line in reader:
                    writer.write(line)
            run_successfully = True
    else:
        print('Different types of files were found: {files_type}'.format(files_type=files_type))

    return run_successfully


'''
def decompress(path_original_compressed_file, path_decompressed_file):
    reader = open_fastq_files([path_original_compressed_file])[0]

    with open(path_decompressed_file, 'at') as writer:
        for line in reader:
            writer.write(line)


def compress(path_file_to_compress):
    with gzip.open(path_file_to_compress + '.gz', 'wt') as writer:
        with open(path_file_to_compress, 'rt') as reader:
            for line in reader:
                writer.write(line)
'''


'''
def kill_subprocess_popen(subprocess_popen, command):
    print('Command run out of time: ' + str(command))
    subprocess_popen.kill()


def run_command_popen_communicate(command, shell_true, timeout_sec_none, print_comand_true):
    run_successfully = False
    if not isinstance(command, str):
        command = ' '.join(command)
    command = shlex.split(command)

    if print_comand_true:
        print('Running: ' + ' '.join(command))

    if shell_true:
        command = ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
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
    return run_successfully, stdout, stderr
'''


def save_variable_to_pickle(variable_to_store, outdir, prefix):
    pickle_file = os.path.join(outdir, str(prefix + '.pkl'))
    with open(pickle_file, 'wb') as writer:
        pickle.dump(variable_to_store, writer)


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            print('Exception in ' + func.__name__)
            print(e)
            traceback.print_exc()

    return wrapped_func


@trace_unhandled_exceptions
def save_single_sample_fastq(sample, dict_sample_files, outdir, fastq_ends):
    print('\n' + '\n' + sample)

    save_problems = True

    files_path_list = {0: [], 1: [], 2: []}

    for fastq_file_name in dict_sample_files:
        for file_properties in dict_sample_files[fastq_file_name]:
            files_path_list[read_pair_0_1_2(fastq_file_name, fastq_ends)].append(file_properties['file_path'])

    if len(files_path_list[1]) != len(files_path_list[2]) or (
            len(files_path_list[0]) > 0 and len(files_path_list[1]) > 0):
        print('ERROR: ' + sample + ' sample with pair-end files problems!')
    else:
        for paired_information in files_path_list:
            if paired_information == 0:
                new_file_name = sample + '.fastq'
            else:
                new_file_name = sample + '_' + str(paired_information) + '.fastq'

            if len(files_path_list[paired_information]) == 0:
                save_problems = False
                continue
            elif len(files_path_list[paired_information]) == 1:
                link_path = os.path.join(outdir, str(
                    new_file_name + '.' + files_path_list[paired_information][0].rsplit('.', 1)[1]))
                if os.path.islink(link_path):
                    os.remove(link_path)
                os.symlink(files_path_list[paired_information][0], link_path)
                save_problems = False
            else:
                # combined_file_path = os.path.join(outdir, new_file_name)

                # import time

                # for fastq_file in sorted(files_path_list[paired_information]):
                #     start = time.time()
                #     decompress(fastq_file, combined_file_path)
                #     end = time.time()
                #     print(fastq_file, end - start, os.path.getsize(combined_file_path))
                #
                # start = time.time()
                # compress(combined_file_path)
                # end = time.time()
                # print(combined_file_path, end - start)

                combined_file_path = os.path.join(outdir, new_file_name + '.gz')

                save_problems = not concatenate_compressed_files(sorted(files_path_list[paired_information]),
                                                                 combined_file_path)
                if save_problems:
                    print('ERROR: it was not possible to concatenate the fastq files'
                          ' for {sample} sample!'.format(sample=sample))
                    if os.path.isfile(combined_file_path):
                        os.remove(combined_file_path)

    save_variable_to_pickle(save_problems, outdir, str(sample + '_save_problems'))


def extract_variable_from_pickle(pickle_file):
    with open(pickle_file, 'rb') as reader:
        variable = pickle.load(reader)
    return variable


def save_sample_fastq(dict_sample_files, outdir, threads, fastq_ends):
    save_problems = False

    pool = multiprocessing.Pool(processes=threads)
    for sample in dict_sample_files:
        pool.apply_async(save_single_sample_fastq, args=(sample, dict_sample_files[sample], outdir, fastq_ends,))
    pool.close()
    pool.join()

    with open(os.path.join(outdir, 'sample_with_problems.txt'), 'wt') as writer:
        files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
        for file_found in files:
            if file_found.endswith('_save_problems.pkl'):
                file_path = os.path.join(outdir, file_found)

                sample_save_problems = extract_variable_from_pickle(file_path)
                if sample_save_problems:
                    save_problems = True
                    writer.write(str(file_found).rsplit('_', 2)[0])

                os.remove(file_path)

    return save_problems


def run_concatenate_fastq_files(args):
    threads = args.threads
    indir = os.path.abspath(args.indir)
    outdir = os.path.abspath(args.outdir)

    check_create_directory(outdir)

    if indir == outdir:
        sys.exit('ERROR: input directory and output directory cannot be the same!')

    # Start logger
    _, _ = start_logger(outdir)

    fastq_files = gather_files(indir, args.fastq_ends)
    sample_files = organize_sample_files(fastq_files, args.sampleNameDelimiter, args.fastq_ends)
    saving_problems = save_sample_fastq(sample_files, outdir, threads, args.fastq_ends)

    if saving_problems:
        sys.exit('ERROR: some errors occur when saving the combined files! Search for ERROR to find what happened.')

    print('\n' + 'concatenate_fastq_files.py is FINISHED')


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 concatenate_fastq_files.py"')

    parser = argparse.ArgumentParser(prog='python3 concatenate_fastq_files.py',
                                     description="Concatenate compressed fastq files from the same sample",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-i', '--indir', type=str, metavar='/path/to/input/directory/',
                                 help='Path for input directory containing the fastq files or folders with the fastq'
                                      ' files inside',
                                 required=True, default='./')
    parser_required.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                 help='Path for output directory', required=True)

    parser_optional = parser.add_argument_group('Facultative options')
    parser_optional.add_argument('-s', '--sampleNameDelimiter', type=str, metavar='_',
                                 help='First left character in sample name delimiting sample name and other sequencing'
                                      ' information',
                                 required=False, default='_')
    parser_optional.add_argument('-j', '--threads', metavar='N', type=int, help='Number of threads to be used',
                                 required=False, default=1)
    parser.add_argument('--fastq_ends', nargs=2, type=str, metavar=('_R1_001.fastq.gz', '_R2_001.fastq.gz'),
                        help='By default, ena_submission.py searches for pair-end fastq files ending with'
                             ' "_R1_001.fastq.gz" and "_R2_001.fastq.gz". If your fastq files end differently, you can'
                             ' provide two strings containing the end of fastq files names (for example, "_1.fastq.gz"'
                             ' and "_2.fastq.gz")',
                        required=False, default=['_R1_001.fastq.gz', '_R2_001.fastq.gz'])

    parser.set_defaults(func=run_concatenate_fastq_files)

    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
