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
from py3helpers.utils import list_dir, captured_output


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

                self.assertTrue(filecmp.cmp(self.correct_per_path, output_per_path))
                self.assertTrue(filecmp.cmp(self.correct_per_read, output_per_read))

    def test_generate_master_kmer_table(self):
        with captured_output() as (_, _):
            heap_size = 1000
            alphabet = "ATGCE"
            n_threads = 2
            with tempfile.TemporaryDirectory() as temp_dir:
                output_path = os.path.join(temp_dir, "out_file.tsv")
                out_file = bindings.generate_master_kmer_table(list_dir(self.assignment_dir, ext="tsv"), output_path,
                                                               heap_size,
                                                               alphabet, n_threads)
                self.assertEqual(192835, os.stat(out_file).st_size)

                out_file2 = bindings.generate_master_kmer_table(list_dir(self.alignment_dir, ext="tsv")[0:1], output_path,
                                                                heap_size,
                                                                "ATGCEF", n_threads=n_threads)
                self.assertEqual(80496, os.stat(out_file2).st_size)
                # self.assertRaises(RuntimeError, bindings.generate_master_kmer_table,
                #                       list_dir(self.alignment_dir, ext="tsv"),
                #                       output_path, heap_size, "ATGCEF", n_threads)


if __name__ == '__main__':
    unittest.main()
