#!/usr/bin/env python
"""Embed signal reads with fastq and event table"""
########################################################################
# File: embed_fast5s
#
# Authors: Andrew Bailey
#
# History: 04/03/19
########################################################################

from argparse import ArgumentParser
from embed.embedding_helpers import *
from embed.fast5 import *
from py3helpers.utils import time_it


# TODO update documentation and add option for adding fastq and/or event table
def parse_args():
    parser = ArgumentParser(description=__doc__)

    # parsers for running the full pipeline
    # required arguments
    parser.add_argument('--fast5_dir', '-f', action='store',
                        dest='fast5_dir', required=True, type=str, default=None,
                        help="path to signal fast5 files")
    parser.add_argument('--fastq', '-q', action='store',
                        dest='fastq', required=True, type=str, default=None,
                        help="path to fastq file")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=2, type=int, help="number of jobs to run in parallel")
    parser.add_argument('--debug', '-d', action='store_true', dest='debug', required=False,
                        help="Option to not multiprocess reads to check for errors")
    parser.add_argument('--embed_build_dir', '-m', action='store',
                        dest='embed_build_dir', required=True, type=str,
                        help="Directory where the main_cpp executable is located")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str,
                        help="Output directory for edited fastq if it is in wrong format")
    parser.add_argument('--no_fastq', action='store_true',
                        dest='no_fastq', required=False,
                        help="If set this will not embed fastqs into the files")
    parser.add_argument('--no_events', action='store_true',
                        dest='no_events', required=False,
                        help="If set this will not embed event tables into the files")
    parser.add_argument('--no_index', action='store_true',
                        dest='no_index', required=False,
                        help="If set this will not generate index files")

    args = parser.parse_args()
    return args


def main():
    # parse args
    args = parse_args()
    # some fastqs have bad input formats which messes up nanopolish parsing
    if not check_fast5_extension_on_fastq(args.fastq):
        output_fastq = os.path.join(args.output_dir, os.path.basename(args.fastq))
        edit_fastq_input(args.fastq, output_fastq)
        args.fastq = output_fastq

    # get full paths to files
    if os.path.exists(args.fast5_dir):
        args.fast5_dir = os.path.abspath(args.fast5_dir)
    if os.path.exists(args.fastq):
        args.fastq = os.path.abspath(args.fastq)

    # index fastq files to fast5s
    if not args.no_index:
        call_nanopolish_index(args.fast5_dir, args.fastq)

    assert os.path.exists(args.fastq+".index"), "Indexing did not work. Please use nanopolish to index files"
    assert os.path.exists(args.fastq+".index.readdb"), "Indexing did not work. Please use nanopolish to index files"

    # will not embed fastqs if set
    if not args.no_fastq:
        if args.debug:
            n_processed = embed_fastqs_into_fast5s(args.fastq, [args.fast5_dir])
        else:
            n_processed = multiprocess_embed_fastqs_into_fast5s(args.fastq, [args.fast5_dir],
                                                                worker_count=args.nb_jobs)
        print("[embed_fast5s] {} reads with Fastqs embedded".format(n_processed))

    # will not embed events if set
    if not args.no_events:
        # sometimes nanopolish readdb does not find the fast5 files
        os.chdir(args.fast5_dir)
        call_embed_main(args.fastq)


if __name__ == "__main__":
    _, run_time = time_it(main)
    print("Running Time = {} seconds".format(run_time))
