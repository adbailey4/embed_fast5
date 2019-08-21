#!/usr/bin/env python3
"""
    Place unit tests for splitting mutli read fast5 files into individual files
"""
########################################################################
# File: test_split_multi_read.py
#  executable: test_split_multi_read.py
#
# Author: Andrew Bailey
# History: 12/18/2019 Created
########################################################################

import unittest
import tempfile
import shutil
import os
from embed.split_multi_read import *
from py3helpers.utils import list_dir


class TestSplitMultiRead(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestSplitMultiRead, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])

        cls.test_create_file = os.path.join(cls.HOME, "tests/test_files/test_create_file.fast5")


        multi_fast5_with_fastq = os.path.join(cls.HOME, "tests/test_files/multi_read_test.fast5")
        multi_fast5_without_fastq = os.path.join(cls.HOME, "tests/test_files/multi_read_test_no_fastq.fast5")

        with tempfile.TemporaryDirectory() as tempdir:
            cls.test_file_with_fastq = os.path.join(tempdir, "test_file.fast5")
            shutil.copyfile(multi_fast5_with_fastq, cls.test_file_with_fastq)

            cls.test_file_with_out_fastq = os.path.join(tempdir, "test_file2.fast5")
            shutil.copyfile(multi_fast5_without_fastq, cls.test_file_with_out_fastq)

            cls.fast5handle_with_fastq = MultiFast5(cls.test_file_with_fastq, 'r+')
            cls.fast5handle_no_fastq = MultiFast5(cls.test_file_with_out_fastq, 'r+')
        cls.new_fast5handle = NewFast5File(cls.test_create_file)

    def test_multiprocess_generate_individual_reads(self):
        with tempfile.TemporaryDirectory() as tempdir:
            fast5_file_dir = os.path.join(self.HOME, "tests/test_files/")
            total_n_processed, total_error = multiprocess_generate_individual_reads(fast5_file_dir, tempdir,
                                                                                    worker_count=2)
            self.assertEqual(total_n_processed, 10)
            self.assertEqual(total_error, 10)
            all_files = list_dir(os.path.join(tempdir, "multi_read_test_no_fastq"), ext="fast5")
            self.assertEqual(5, len(all_files))
            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))
            all_files = list_dir(os.path.join(tempdir, "multi_read_test"), ext="fast5")
            self.assertEqual(5, len(all_files))
            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))

    def test_generate_individual_reads(self):
        with tempfile.TemporaryDirectory() as tempdir:
            fast5_file_dir = os.path.join(self.HOME, "tests/test_files/")
            total_n_processed, total_error = generate_individual_reads(fast5_file_dir, tempdir)
            self.assertEqual(total_n_processed, 10)
            self.assertEqual(total_error, 10)
            all_files = list_dir(os.path.join(tempdir, "multi_read_test_no_fastq"), ext="fast5")
            self.assertEqual(5, len(all_files))
            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))
            all_files = list_dir(os.path.join(tempdir, "multi_read_test"), ext="fast5")
            self.assertEqual(5, len(all_files))
            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))

    def test_write_individual_fast5s(self):
        with tempfile.TemporaryDirectory() as tempdir:
            self.fast5handle_with_fastq.write_individual_fast5s(tempdir)
            all_files = list_dir(tempdir)
            self.assertEqual(5, len(all_files))
            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))

        with tempfile.TemporaryDirectory() as tempdir:
            self.fast5handle_with_fastq.write_individual_fast5s(tempdir, write_fastq_file=os.path.join(tempdir, "test.fastq"))
            all_files = list_dir(tempdir)
            self.assertEqual(6, len(all_files))

            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))

        with tempfile.TemporaryDirectory() as tempdir:
            self.fast5handle_no_fastq.write_individual_fast5s(tempdir)
            all_files = list_dir(tempdir)
            self.assertEqual(5, len(all_files))
            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))

        with tempfile.TemporaryDirectory() as tempdir:
            self.fast5handle_no_fastq.write_individual_fast5s(tempdir, write_fastq_file=os.path.join(tempdir, "test.fastq"))
            all_files = list_dir(tempdir)
            self.assertEqual(6, len(all_files))

            self.assertEqual(5, sum([int(os.path.basename(x).startswith("read")) for x in all_files]))

    def test_get_basecall_1d_attributes(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        bc_attrs = self.fast5handle_with_fastq.get_basecall_1d_attributes(get_id)
        self.assertEqual(3, len(bc_attrs.keys()))

    def test_get_fastq(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        fastq = self.fast5handle_with_fastq.get_fastq(get_id)
        self.assertEqual(22042, len(fastq))

    def test_get_signal(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        signal = self.fast5handle_with_fastq.get_signal(get_id)
        self.assertEqual(684, signal[0])
        self.assertEqual(578, signal[-1])

    def test_get_signal_attrbutes(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        signal = self.fast5handle_with_fastq.get_signal_attrbutes(get_id)
        self.assertTrue("duration" in signal.keys())
        self.assertTrue("read_id" in signal.keys())
        self.assertTrue("start_mux" in signal.keys())
        self.assertTrue("read_number" in signal.keys())
        self.assertTrue("start_time" in signal.keys())

    def test_get_channel_id_attrbutes(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        ch_attrs = self.fast5handle_with_fastq.get_channel_id_attrbutes(get_id)
        self.assertSetEqual({"digitisation", "offset", "range", "sampling_rate", "channel_number"},
                            set(ch_attrs.keys()))

    def test_get_context_tags_attrbutes(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        ct_attrs = self.fast5handle_with_fastq.get_context_tags_attrbutes(get_id)
        self.assertEqual(15, len(ct_attrs.keys()))

    def test_get_tracking_id_attrbutes(self):
        get_id = "read_002f9702-c19e-48c2-8e72-9021adbd4a48"
        ti_attrs = self.fast5handle_with_fastq.get_tracking_id_attrbutes(get_id)
        self.assertEqual(32, len(ti_attrs.keys()))

    def test_write_signal(self):
        attributes = {"duration": 241066, "median_before": 245.39717658547795,
                      "read_id": "104aa826-0469-4781-8c88-d6f88049f09c",
                      "read_number": 26,
                      "start_mux": 1,
                      "start_time": 224439}
        self.new_fast5handle.write_signal(np.asarray([10, 20, 30]), attributes)
        self.assertTrue(self.new_fast5handle[self.new_fast5handle.__default_raw_read__], [10, 20, 30])

    def test_write_1d_fastq(self):
        attributes = {"name": "MinKNOW-Live-Basecalling",
                      "time_stamp": "2019-01-23T01:36:38Z",
                      "version": "3.1.18"}
        self.new_fast5handle.write_1d_fastq("@Some_fastq string \n and more data", attributes)
        self.assertTrue(self.new_fast5handle['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'],
                        "@Some_fastq string \n and more data")

    def test_write_unique_global_key(self):
        context_tags = {"channel_number": 458,
                        "digitisation": 8192.0,
                        "offset": 8.0,
                        "range": 1469.85,
                        "sampling_rate": 4000.0}
        self.new_fast5handle.write_unique_global_key(context_tags, context_tags, context_tags)
        self.assertTrue(self.new_fast5handle[self.new_fast5handle.__default_global_key__], context_tags)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.test_create_file)


if __name__ == '__main__':
    unittest.main()
