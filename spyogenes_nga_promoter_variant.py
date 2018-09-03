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
import shlex
import subprocess
import sys
from threading import Timer

from Bio.Blast import NCBIXML


version = '0.1'


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
    return run_successfully, stdout.decode("utf-8"), stderr.decode("utf-8")


def check_db_exists(db_path):
    """
    Check if Blast DB exists

    Parameters
    ----------
    db_path : str
        Path to Blast DB files, e.g. /output/blast_db.nucl.sequences.fasta

    Returns
    -------
    db_exists : bool
        Tells if the Blast DB already exists
    original_file : bool
        Tells if the original fasta file from which the Blast DB was produced is present
    """
    db_exists = False
    files = [f for f in os.listdir(os.path.dirname(db_path))
             if not f.startswith('.')
             and os.path.isfile(os.path.join(os.path.dirname(db_path), f))]
    counter = 0
    original_file = False
    for file_found in files:
        if file_found == os.path.basename(db_path):
            original_file = True
        elif os.path.splitext(file_found)[0] == os.path.basename(db_path):
            counter += 1
    if counter > 0:
        db_exists = True
    return db_exists, original_file


def create_blast_db(db_sequences, db_output, db_type):
    """
    Creates a Blast DB

    Parameters
    ----------
    db_sequences : str
        Path to fasta file containing the sequences from which Blast DB will be created, e.g. /input/sequences.fasta
    db_output : str
        Path to Blast output DB files, e.g. /output/blast_db.nucl.sequences.fasta
    db_type : str
        Blast DB type. Can only be 'nucl' or 'prot'

    Returns
    -------
    run_successfully : bool
        Tells if the Blast DB was successfully created
    """

    run_successfully = False

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}). Use one of the following: nucl, prot'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(db_output)):
        os.makedirs(os.path.dirname(db_output))
    else:
        db_exists, original_file = check_db_exists(db_output)
        if db_exists and original_file:
            print('Blast DB already found at {db_output}'
                  ' for {db_sequences}'.format(db_output=os.path.dirname(db_output),
                                               db_sequences=db_sequences,
                                               file_found=os.path.basename(db_sequences)))

            run_successfully = True
        elif db_exists and not original_file:
            print('The original fasta file ({db_sequences}) from which the Blast DB was produced is not present. Make'
                  ' sure it is found in {db_dir} and it is'
                  ' named {original_file_name}'.format(db_dir=os.path.dirname(db_output), db_sequences=db_sequences,
                                                       original_file_name=os.path.basename(db_output)))
        elif not db_exists and original_file:
            print('The original fasta file ({db_sequences}) from which the Blast DB was supposed to be produced is'
                  ' present in Blast DB directory ({db_dir}), but no Blast DB files were found'
                  ' there.'.format(db_sequences=db_sequences, db_dir=os.path.dirname(db_output)))
        else:
            run_successfully, _, _ = run_command_popen_communicate(['makeblastdb', '-parse_seqids', '-dbtype',
                                                                    db_type, '-in', db_sequences, '-out', db_output],
                                                                   False, None, True)
            if run_successfully:
                from shutil import copyfile
                copyfile(db_sequences, db_output)

    return run_successfully


def run_blast_command(query_file, blast_db, db_type, blast_output, threads=1):
    """
    Run Blast: blastn or blastp

    Parameters
    ----------
    query_file : str
        Path to fasta file containing the query sequences, e.g. /input/queries.fasta
    blast_db : str
        Path to Blast DB files, e.g. /input/blast_db.nucl.sequences.fasta
    db_type : str
        Blast DB type. Can only be 'nucl' or 'prot'
    blast_output : str
        Path to Blast output tabular file
    threads : int
        Number of CPUs to use during Blast

    Returns
    -------
    run_successfully : bool
        Tells if the Blast DB was successfully created
    """

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}).'
             ' Use one of the following: "nucl", "prot"'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(blast_output)):
        os.makedirs(os.path.dirname(blast_output))

    command = ['', '-query', query_file, '-db', blast_db, '-out', blast_output, '-outfmt', '5', '-dust', 'no',
               '-culling_limit', '1', '-num_threads', str(threads)]

    if db_type == 'nucl':
        command[0] = 'blastn -task blastn'
    else:
        command[0] = 'blastp -task blastp'

    run_successfully, _, _ = run_command_popen_communicate(command, False, None, True)

    return run_successfully


def run_blast(blast_db_path, outdir, blast_type, query_fasta_file):
    """
    Parse Blast output

    Parameters
    ----------
    blast_db_path : str
        Path to Blast DB files, e.g. /input/blast_db.nucl.sequences.fasta. If Blast DB is not found, it will be created
        under the outdir folder
    outdir : str
        Path to output directory
    blast_type : str
        Blast DB type. Can only be 'nucl' or 'prot'
    query_fasta_file : str
        Path to fasta file containing the query sequences, e.g. /input/queries.fasta

    Returns
    -------
    folders_2_remove : list
        List of folders that can be removed at the end
    blast_results : dict
        Dictionary with Blast results cleaned and almost already formated for parse_results.py: subject as key and
        values similar to ReMatCh results
    blast_db_path : str
        Path to Blast DB used
    """

    folders_2_remove = []

    # Check Blast DB
    db_exists, original_file = check_db_exists(blast_db_path)
    if not db_exists:
        blast_db = os.path.join(outdir, 'blast_db', '{blast_DB}'.format(blast_DB=os.path.basename(blast_db_path)))
        folders_2_remove.append(os.path.dirname(blast_db))

        if not os.path.isdir(os.path.dirname(blast_db)):
            os.makedirs(os.path.dirname(blast_db))

        db_exists = create_blast_db(blast_db_path, blast_db, blast_type)
        if db_exists:
            blast_db_path = str(blast_db)
            original_file = True
    elif db_exists and not original_file:
        sys.exit('Original fasta file from which the Blast DB was produced ({blast_db_path}) is missing'
                 ' from {blast_db_dir}'.format(blast_db_path=os.path.basename(blast_db_path),
                                               blast_db_dir=os.path.dirname(blast_db_path)))

    # Run Blast
    blast_output = os.path.join(outdir, 'blast_out', 'results.xml')
    if db_exists and original_file:
        folders_2_remove.append(os.path.dirname(blast_output))

        if not os.path.isdir(os.path.dirname(blast_output)):
            os.makedirs(os.path.dirname(blast_output))
        run_successfully = run_blast_command(query_fasta_file, blast_db_path, blast_type, blast_output)

        if not run_successfully:
            sys.exit('Blast was not run successfully')
    else:
        sys.exit('It was not found any Blast DB and/or the original fasta file from which the Blast DB was produced')

    return blast_output


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 spyogenes_nga_promoter_variant.py"')

    parser = argparse.ArgumentParser(prog='spyogenes_nga_promoter_variant.py',
                                     description='Determines Streptococcus pyogenes nga gene promoter variant using'
                                                 ' fasta sequences',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-r', '--reference', nargs=1, type=argparse.FileType('r'), required=True,
                                 metavar='/path/to/reference_sequence.fasta',
                                 help='Fasta file containing nga gene and the promoter sequence'
                                      ' ("spyogenes_nga_promoter_variant.CP000017_nga_promoter.fasta").')
    parser_required.add_argument('-f', '--fasta', nargs=1, type=argparse.FileType('r'),
                                 metavar='/path/to/query/assembly_file.fasta',
                                 help='Path to fasta file containing the query sequences from which the'
                                      ' nga gene promoter variant should be assessed',
                                 required=True)

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the information will be stored',
                                         required=False, default='.')
    parser_optional_general.add_argument('-t', '--type', choices=['nucl', 'prot'], type=str, metavar='nucl',
                                         help='Blast DB type (available options: %(choices)s)')

    args = parser.parse_args()

    reference = os.path.abspath(args.reference[0].name)
    fasta = os.path.abspath(args.fasta[0].name)
    outdir = os.path.abspath(args.outdir)
    blast_type = args.type

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    blast_output = run_blast(blast_db_path=fasta, outdir=outdir, blast_type=blast_type, query_fasta_file=reference)

    blast_results = {}
    with open(blast_output, 'rt') as results:
        records = NCBIXML.parse(results)
        for query in records:
            # query_name = query.query
            # query_length = query.query_length
            best_hit = None
            for alignment in query.alignments:
                # subject_name = alignment.accession
                # subject_length = alignment.length
                for hsp in alignment.hsps:
                    if best_hit is None:
                        best_hit = (alignment.accession, hsp)
                    else:
                        to_change = False

                        if best_hit[1].align_length < hsp.align_length:
                            to_change = True
                        elif best_hit[1].align_length == hsp.align_length and \
                                best_hit[1].identities <= hsp.identities:
                            if best_hit[1].identities < hsp.identities:
                                to_change = True
                            elif best_hit[1].identities == hsp.identities and \
                                    best_hit[1].expect >= hsp.expect:  # expect for e-value
                                if best_hit[1].expect > hsp.expect:
                                    to_change = True
                                elif best_hit[1].expect == hsp.expect and \
                                        best_hit[1].gaps > hsp.gaps:
                                    to_change = True

                        if to_change:
                            best_hit = (alignment.accession, hsp)

            blast_results[query.query] = (query.query_length, best_hit)

    variants = {'TTGACAAAAACTATTTGCTAATGTATAGT': 1,
                'TTGACAAAGACTATTTGTTAATGTATAGT': 2,
                'TTGACAAAAACTAGTTGTTAATGTATAGT': 3}
    min_gene_coverage = 90

    result = ('No results found', None, None, None, None)

    for query, info in blast_results.items():
        query_length = info[0]
        subject_name = info[1][0]
        hsp = info[1][1]
        if hsp is not None:
            if hsp.frame[0] < 0:
                sys.exit('Unexpected reference (query) frame')
            else:
                if hsp.align_length / float(query_length) * 100 >= min_gene_coverage:
                    var_found = None
                    for var_seq, var_id in variants.items():
                        if var_seq in hsp.query:
                            if hsp.query.find(var_seq) + len(var_seq) + 1 < 45:
                                var_found = var_id

                    if var_found is not None:
                        result = ('Variant known', var_found, subject_name, hsp.sbjct_start, hsp.sbjct_end)
                    else:
                        result = ('New variant', None, subject_name, hsp.sbjct_start, hsp.sbjct_end)
                else:
                    result = ('Partial sequence found', None, None, None, None)

    with open(os.path.join(outdir, 'spyogenes_nga_promoter_variant.result.tab'), 'wt') as writer:
        writer.write('\t'.join(['#result', 'variant', 'sequence', 'start_position', 'end_position']) + '\n' +
                     '\t'.join(map(str, result)) + '\n')

    # print(best_hit.align_length, best_hit.query_start, best_hit.query_end, best_hit.sbjct_start,
    #       best_hit.sbjct_end, best_hit.frame)
    # print("****Alignment****")
    # print(best_hit.query[0:75] + "...")
    # print(best_hit.match[0:75] + "...")
    # print(best_hit.sbjct[0:75] + "...")


if __name__ == "__main__":
    main()
