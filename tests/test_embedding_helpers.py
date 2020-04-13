#!/usr/bin/env python3
"""
    Place unit tests for the embedding_helpers.py library functions
"""
########################################################################
# File: test_embedding_helpers.py
#  executable: test_embedding_helpers.py
#
# Author: Andrew Bailey
# History: 04/23/2019 Created
########################################################################

import unittest
import tempfile
import shutil
import sys
import os
import subprocess
from embed.embedding_helpers import *
from embed.fast5 import Fast5
from py3helpers.utils import list_dir, captured_output


class TestSplitMultiRead(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestSplitMultiRead, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.test_files_dir = os.path.join(cls.HOME, "tests/test_files/test_readdb/")
        cls.read_db = os.path.join(cls.HOME, "tests/test_files/test_readdb/new_individual.fastq.index.readdb")
        cls.fastq = os.path.join(cls.HOME, "tests/test_files/test_readdb/new_individual.fastq")
        cls.fast5 = os.path.join(cls.HOME,
                                 "tests/test_files/test_readdb/read_002f9702-c19e-48c2-8e72-9021adbd4a48.fast5")

        cls.r94_fastq = os.path.join(cls.HOME, "tests/test_files/r94_tests/small_sample.fastq")
        cls.r94_dir = os.path.join(cls.HOME, "tests/test_files/r94_tests")

    def test_parse_readdb(self):
        with captured_output() as (_, _):
            total_reads = 0
            for read_id, fast5_path in parse_readdb(self.read_db, [self.test_files_dir]):
                total_reads += 1
                self.assertTrue(os.path.exists(fast5_path))
            self.assertEqual(1, total_reads)

    def test_call_nanopolish_index(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                fastq_file = os.path.join(tempdir, "test.fastq")
                shutil.copyfile(self.fastq, fastq_file)
                fast5_file = os.path.join(tempdir, "test.fast5")
                shutil.copyfile(self.fast5, fast5_file)

                works = call_nanopolish_index(tempdir, fastq_file)
                self.assertTrue(works)
                output = list_dir(tempdir)
                self.assertEqual(len(output), 6)

    def test_call_embed_main(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                fastq_file = os.path.join(tempdir, "test.fastq")
                shutil.copyfile(self.fastq, fastq_file)
                fast5_file = os.path.join(tempdir, "test.fast5")
                shutil.copyfile(self.fast5, fast5_file)
                # check if fast5 file has data
                fh = Fast5(fast5_file)
                self.assertRaises(IndexError, fh.get_basecall_data)
                self.assertRaises(IndexError, fh.get_fastq)
                fh.close()

                works = call_nanopolish_index(tempdir, fastq_file)
                self.assertTrue(works)
                works = call_embed_main(fastq_file)
                self.assertTrue(works)
                output = list_dir(tempdir)
                self.assertEqual(len(output), 6)
                # check if fast5 file has data
                fh = Fast5(fast5_file)
                data = fh.get_basecall_data()
                self.assertEqual(len(data), 20561)

    def test_check_fast5_extension_on_fastq(self):
        with captured_output() as (_, _):
            check_passed = check_fast5_extension_on_fastq(self.fastq)
            self.assertTrue(check_passed)
            check_passed = check_fast5_extension_on_fastq(self.r94_fastq)
            self.assertFalse(check_passed)

    def test_edit_fastq_input(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                # confirm that the r94 fastq failed test
                check_passed = check_fast5_extension_on_fastq(self.r94_fastq)
                self.assertFalse(check_passed)
                # create new fastq with correct output
                output_fastq = os.path.join(tempdir, os.path.basename(self.r94_fastq))
                out_fastq = edit_fastq_input(self.r94_fastq, output_fastq)
                self.assertTrue(os.path.exists(out_fastq))
                # make sure that new fastq passes
                check_passed = check_fast5_extension_on_fastq(out_fastq)
                self.assertTrue(check_passed)

    def test_embed_fastqs_into_fast5s(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                fastq_file = os.path.join(tempdir, "test.fastq")
                shutil.copyfile(self.fastq, fastq_file)
                readdb_file = os.path.join(tempdir, "test.fastq.index.readdb")
                shutil.copyfile(self.read_db, readdb_file)
                fast5_file = os.path.join(tempdir, "read_002f9702-c19e-48c2-8e72-9021adbd4a48.fast5")
                shutil.copyfile(self.fast5, fast5_file)
                # check if fast5 file has data
                fh = Fast5(fast5_file)
                self.assertRaises(IndexError, fh.get_fastq)
                fh.close()

                works = embed_fastqs_into_fast5s(fastq_file, [tempdir], fastq_readdb=readdb_file)
                self.assertEqual(works, 1)

                fh = Fast5(fast5_file)
                data = fh.get_fastq()
                self.assertEqual(len(data), 22041)
                fh.close()

    def test_multiprocess_embed_fastqs_into_fast5s(self):
        with captured_output() as (_, _):
            with tempfile.TemporaryDirectory() as tempdir:
                fastq_file = os.path.join(tempdir, "test.fastq")
                shutil.copyfile(self.fastq, fastq_file)
                readdb_file = os.path.join(tempdir, "test.fastq.index.readdb")
                shutil.copyfile(self.read_db, readdb_file)
                fast5_file = os.path.join(tempdir, "read_002f9702-c19e-48c2-8e72-9021adbd4a48.fast5")
                shutil.copyfile(self.fast5, fast5_file)
                # check if fast5 file has data
                fh = Fast5(fast5_file)
                self.assertRaises(IndexError, fh.get_fastq)
                fh.close()

                works = multiprocess_embed_fastqs_into_fast5s(fastq_file, [tempdir], fastq_readdb=readdb_file)
                self.assertEqual(works, 1)

                fh = Fast5(fast5_file)
                data = fh.get_fastq()
                self.assertEqual(len(data), 22041)
                fh.close()

    def test_embed_fast5_script(self):
        with captured_output() as (stdout, stderr):
            with tempfile.TemporaryDirectory() as tempdir:
                temp_r94_dir = os.path.join(tempdir, "r94_tests")
                shutil.copytree(self.r94_dir, temp_r94_dir)
                fast5_dir = os.path.join(temp_r94_dir, "fast5")
                fastq = os.path.join(temp_r94_dir, "small_sample.fastq")

                test_fast5_file = "DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_88738_ch1_read1464_strand.fast5"
                fast5_file_path = os.path.join(fast5_dir, test_fast5_file)

                fh = Fast5(fast5_file_path)
                self.assertRaises(IndexError, fh.get_fastq)
                self.assertRaises(IndexError, fh.get_basecall_data)
                fh.close()

                embed_command = "embed_fast5s.py --fast5_dir {} --fastq {} " \
                                "--jobs {} --output_dir {}".format(fast5_dir,
                                                                   fastq,
                                                                   1,
                                                                   tempdir)
                command = embed_command.split()
                proc = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                proc.communicate()

                fh = Fast5(fast5_file_path)
                fastq_string = fh.get_fastq()
                events = fh.get_basecall_data()
                self.assertEqual(len(events), 1700)
                self.assertEqual(len(fastq_string), 2046)
                fh.close()


if __name__ == '__main__':
    unittest.main()
