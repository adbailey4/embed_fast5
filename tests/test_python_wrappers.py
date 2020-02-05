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


class TestPythonWrappers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestPythonWrappers, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])

    def test_add(self):
        self.assertEqual(3, bindings.add(1, 2))

    def test_subtract(self):
        self.assertEqual(1, bindings.subtract(2, 1))

    def test_LoadVariantPaths(self):
        positions_file = os.path.join(self.HOME,
                                      "tests/test_files/rRNA_test_files/16S_final_branch_points.positions")
        rna_signal_files = os.path.join(self.HOME,
                                        "tests/test_files/rRNA_test_files/rRNA_signal_files")
        correct_per_path = os.path.join(self.HOME,
                                        "tests/test_files/rRNA_test_files/test_output_dir/per_path_counts.tsv")
        correct_per_read = os.path.join(self.HOME,
                                        "tests/test_files/rRNA_test_files/test_output_dir/per_read_calls.tsv")
        lvp = bindings.LoadVariantPaths(positions_file, rna_signal_files, True, 100)
        with tempfile.TemporaryDirectory() as temp_dir:
            output_per_path = os.path.join(temp_dir, "per_path_counts.tsv")
            output_per_read = os.path.join(temp_dir, "per_read_calls.tsv")
            lvp.write_per_path_counts(output_per_path)
            lvp.write_per_read_calls(output_per_read)

            self.assertTrue(filecmp.cmp(correct_per_path, output_per_path))
            self.assertTrue(filecmp.cmp(correct_per_read, output_per_read))

    def test_generate_master_assignment_table(self):
        assignment_dir = os.path.join(self.HOME,"tests/test_files/assignment_files")
        heap_size = 1000
        alphabet = "ATGC"
        n_threads = 2
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = temp_dir
            a = bindings.generate_master_assignment_table(assignment_dir, output_dir, heap_size, alphabet, n_threads)
            size = os.stat(a).st_size
            self.assertEqual(438960, size)


if __name__ == '__main__':
    unittest.main()
