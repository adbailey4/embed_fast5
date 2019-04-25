#!/usr/bin/env python
"""Split multi fast5 file into multiple fast5s"""
########################################################################
# File: split_multi_fast5
#
# Authors: Andrew Bailey
#
# History: 04/03/19
########################################################################

from argparse import ArgumentParser
from embed.split_multi_read import *
from py3helpers.utils import list_dir, time_it


def parse_args():
    parser = ArgumentParser(description=__doc__)

    # parsers for running the full pipeline
    # required arguments
    parser.add_argument('--multi_fast5_dir', '-f', action='store',
                        dest='multi_fast5_dir', required=True, type=str, default=None,
                        help="Path to multi read fast5 directory")

    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str,
                        help="Directory to place all new originally formatted fast5 files")

    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False,
                        default=2, type=int, help="number of jobs to run in parallel")

    parser.add_argument('--debug', '-d', action='store_true', dest='debug', required=False,
                        help="Option to not multiprocess reads to check for errors")

    parser.add_argument('--delete_multi', action='store_true', dest='delete_multi', required=False,
                        help="Option to delete the multi_fast5 file after processing")

    args = parser.parse_args()
    return args


def main():
    # parse args
    args = parse_args()

    if args.debug:
        total_processed = 0
        total_errors = 0
        for multi_fast5_file in list_dir(args.multi_fast5_dir, ext="fast5"):

            fast5_specific_output_dir = os.path.join(args.output_dir, os.path.basename(multi_fast5_file).split(".")[0])
            if not os.path.exists(fast5_specific_output_dir):
                os.mkdir(fast5_specific_output_dir)
            mf5h = MultiFast5(multi_fast5_file)
            n_processed, n_error = mf5h.write_individual_fast5s(fast5_specific_output_dir)
            total_processed += n_processed
            total_errors += total_errors
            if args.delete_multi:
                os.remove(multi_fast5_file)
    else:
        total_processed, total_errors = multiprocess_generate_individual_reads(args.multi_fast5_dir, args.output_dir,
                                                                               worker_count=args.nb_jobs,
                                                                               delete_multi=args.delete_multi)

    print("Total Processed: {}\nTotal Errors: {}".format(total_processed, total_errors))


if __name__ == "__main__":
    _, run_time = time_it(main)
    print("Running Time = {} seconds".format(run_time))
