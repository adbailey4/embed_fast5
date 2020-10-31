#!/usr/bin/env python3
"""
    Place unit tests for python wrapped c++ functions
"""
########################################################################
# File: test_embedding_helpers.py
#  executable: test_embedding_helpers.py
#
# Author: Andrew Bailey
# History: 07/30/2019 Created
########################################################################

import unittest
import tempfile
import os
import filecmp
from embed import bindings
from py3helpers.utils import list_dir, captured_output, count_lines_in_file


class TestPythonWrappers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestPythonWrappers, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.positions_file = os.path.join(cls.HOME,
                                          "tests/test_files/rRNA_test_files/16S_final_branch_points.positions")
        cls.rna_signal_files = os.path.join(cls.HOME,
                                            "tests/test_files/rRNA_test_files/rRNA_signal_files")
        cls.correct_per_path = os.path.join(cls.HOME,
                                            "tests/test_files/rRNA_test_files/test_output_dir/per_path_counts.tsv")
        cls.correct_per_read = os.path.join(cls.HOME,
                                            "tests/test_files/rRNA_test_files/test_output_dir/per_read_calls.tsv")
        cls.assignment_dir = os.path.join(cls.HOME, "tests/test_files/assignment_files")
        cls.alignment_dir = os.path.join(cls.HOME, "tests/test_files/alignment_files")

    def test_LoadVariantPaths(self):
        with captured_output() as (_, _):
            lvp = bindings.LoadVariantPaths(self.positions_file, self.rna_signal_files, True, 100)
            with tempfile.TemporaryDirectory() as temp_dir:
                output_per_path = os.path.join(temp_dir, "per_path_counts.tsv")
                output_per_read = os.path.join(temp_dir, "per_read_calls.tsv")
                lvp.write_per_path_counts(output_per_path)
                lvp.write_per_read_calls(output_per_read)
                with open(self.correct_per_path, 'r') as fh, open(output_per_path) as fh2:
                    self.assertSetEqual(set(fh.readlines()), set(fh2.readlines()))
                with open(self.correct_per_read, 'r') as fh, open(output_per_read) as fh2:
                    self.assertSetEqual(set(fh.readlines()), set(fh2.readlines()))

    def test_generate_master_kmer_table(self):
        with captured_output() as (_, _):
            heap_size = 1000
            alphabet = "ATGCE"
            n_threads = 2
            with tempfile.TemporaryDirectory() as temp_dir:
                output_path = os.path.join(temp_dir, "out_file.tsv")
                log_file_path = os.path.join(temp_dir, "log_file.tsv")

                bindings.generate_master_kmer_table(list_dir(self.assignment_dir, ext="tsv"), output_path,
                                                    log_file_path,
                                                    heap_size,
                                                    alphabet, min_prob=0, n_threads=n_threads)
                self.assertEqual(15626, count_lines_in_file(log_file_path))
                self.assertEqual(17350, count_lines_in_file(output_path))
                bindings.generate_master_kmer_table(
                    [os.path.join(self.alignment_dir, "7d31de25-8c15-46d8-a08c-3d5043258c89.sm.forward.tsv")],
                    output_path, log_file_path, heap_size, "ATGCEF", n_threads=n_threads)
                self.assertEqual(7777, count_lines_in_file(log_file_path))
                self.assertEqual(1565, count_lines_in_file(output_path))

                bindings.generate_master_kmer_table(
                    [os.path.join(self.alignment_dir, "7d31de25-8c15-46d8-a08c-3d5043258c89.sm.forward.tsv")],
                    output_path, log_file_path, heap_size, "ATGCEF", min_prob=0.5, n_threads=n_threads,
                    full=False)
                self.assertEqual(7777, count_lines_in_file(log_file_path))
                self.assertEqual(417, count_lines_in_file(output_path))

            # self.assertRaises(RuntimeError, bindings.generate_master_kmer_table,
            #                       list_dir(self.alignment_dir, ext="tsv"),
            #                       output_path, heap_size, "ATGCEF", n_threads)


if __name__ == '__main__':
    unittest.main()
