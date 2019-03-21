#!/usr/bin/env python
"""Run whole workflow from a multi-fast5 file to fully embeded original fast5 formats"""
########################################################################
# File: run_embed_fast5.py
#
# Authors: Andrew Bailey
#
# History: 03/20/17
########################################################################

import sys
import h5py
import pandas as pd
import numpy as np
import subprocess
from argparse import ArgumentParser
from embed.split_multi_read import *
from py3helpers.multiprocess import *
from py3helpers.utils import list_dir


def multiprocess_generate_embedded_reads(multi_fast5_dir, out_dir, worker_count=2):
    """Generate embedded reads from multi_fast5 reads

    :return: list of lists of random numbers with max integer
    """
    test_args = {"out_dir": out_dir}
    service = BasicService(multi_fast5_wrapper)
    total, failure, messages, output = run_service(service.run, list_dir(multi_fast5_dir, ext="fast5"),
                                                   test_args, ["multi_fast5"], worker_count)
    return output


def generate_embedded_reads(multi_fast5_dir, out_dir):
    """Generate individual reads from combined reads

    :return: tuple (total_processed_files, total_errors)
    """
    total_n_processed = 0
    total_error = 0
    for multi_fast5 in list_dir(multi_fast5_dir, ext="fast5"):
        print("Processing: {}".format(multi_fast5))
        n_processed, n_error = multi_fast5_wrapper(multi_fast5, out_dir)
        total_error += n_error
        total_n_processed += n_processed
        print("Processed {} reads".format(total_n_processed))

    return total_n_processed, total_error


def generate_embedded_read(multi_fast5, out_dir, fastq_path=None):
    """Generate embedded individual reads from combined reads

    :param fastq_path: boolean option to write a fastq file for the extracted reads
    :param multi_fast5: path to multi fast5 file
    :param out_dir: output directory
    :return: tuple (total_processed_files, total_errors)
    """
    print("Processing: {}".format(multi_fast5))
    n_processed = 0
    n_error = 0
    try:
        mf5h = MultiFast5(multi_fast5)
        n_processed, n_error = mf5h.write_individual_fast5s(out_dir, fastq_path)

    except KeyError:
        pass
    print("Processed {} reads".format(n_processed))
    return n_processed, n_error


def call_nanopolish_index(nanopolish_dir, output_dir, fastq):
    """Call nanopolish index on some files"""
    nanopolish_path = os.path.join(nanopolish_dir, "nanopolish")
    assert os.path.exists(nanopolish_path), "Nanopolish does not exist in directory {}".format(nanopolish_dir)
    nanopolish_command = nanopolish_path + " index -d {output_dir} {fastq}".format(output_dir=output_dir, fastq=fastq)

    try:
        command = nanopolish_command.split()
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        output, errors = proc.communicate()
        errors = errors.decode().splitlines()
        for x in errors:
            print(x)
        output = output.decode().splitlines()
        for x in output:
            print(x)

    except Exception as e:
        print("[run_embed_fast5] exception ({}) running nanopolish extract: {}".format(type(e), e))
        raise e

    return True


def call_embed_main(main_cpp_dir, fastq):
    """Call embed on all files"""
    main_cpp_path = os.path.join(main_cpp_dir, "main_cpp")
    assert os.path.exists(main_cpp_path), "main_cpp_path does not exist in directory {}".format(main_cpp_dir)
    embed_command = main_cpp_path + " embed -r {fastq}".format(fastq=fastq)

    try:
        command = embed_command.split()
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        output, errors = proc.communicate()
        errors = errors.decode().splitlines()
        for x in errors:
            print(x)
        output = output.decode().splitlines()
        for x in output:
            print(x)
    except Exception as e:
        print("[run_embed_fast5] exception ({}) running nanopolish extract: {}".format(type(e), e))
        raise e

    return True


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # parsers for running the full pipeline
    # required arguments
    parser.add_argument('--fast5', '-f', action='store',
                        dest='fast5', required=True, type=str, default=None,
                        help="path to multi read fast5 file")

    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str,
                        help="Directory to place all new originally formatted fast5 files")

    parser.add_argument('--main_cpp_dir', '-m', action='store',
                        dest='main_cpp_dir', required=True, type=str,
                        help="Directory where the main_cpp executable is located")

    parser.add_argument('--nanopolish_dir', '-n', action='store',
                        dest='nanopolish_dir', required=True, type=str,
                        help="Directory where the nanopolish executable is located")

    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=2, type=int, help="number of jobs to run in parallel")

    args = parser.parse_args()
    return args


def main():
    # parse args
    start = timer()

    args = parse_args()

    fastq_path = os.path.join(args.output_dir, os.path.basename(args.fast5).split(".")[0]+".fastq")
    generate_embedded_read(args.fast5, args.output_dir, fastq_path)
    success = call_nanopolish_index(args.nanopolish_dir, args.output_dir, fastq_path)
    if success:
        call_embed_main(args.main_cpp_dir, fastq_path)

    stop = timer()
    print("Running Time = {} seconds".format(stop - start))


if __name__ == "__main__":
    main()
