#!/usr/bin/env python
"""Split a multi-fast5 into multiple reads """
########################################################################
# File: split_multi_read.py
#
# Authors: Andrew Bailey
#
# History: 03/19/17
########################################################################

import os
import h5py
import pandas as pd
import numpy as np
from py3helpers.multiprocess import *
from py3helpers.utils import list_dir


def multiprocess_generate_individual_reads(multi_fast5_dir, out_dir, worker_count=2, delete_multi=False):
    """Generate individual reads from combined reads

    :return: list of lists of random numbers with max integer
    """
    test_args = {"out_dir": out_dir,
                 "delete_multi": delete_multi}
    service = BasicService(multi_fast5_wrapper, service_name="multiprocess_generate_individual_reads")
    total, failure, messages, output = run_service(service.run, list_dir(multi_fast5_dir, ext="fast5"),
                                                   test_args, ["multi_fast5"], worker_count)
    total_n_processed = 0
    total_error = 0
    for n_processed, n_error in output:
        total_error += n_error
        total_n_processed += n_processed
    return total_n_processed, total_error


def multi_fast5_wrapper(multi_fast5, out_dir, delete_multi=False):
    """Wrap MultiFast5 write_individual_fast5s into a single function call
    :param multi_fast5: multi_fast5
    :param out_dir: path to empty output directory
    :param delete_multi: boolean option to delete old multi fast5 after processing
    :return: (number of reads processed, number of errors)
    """
    n_processed = 0
    n_error = 0
    try:
        fast5_specific_output_dir = os.path.join(out_dir, os.path.basename(multi_fast5).split(".")[0])
        if not os.path.exists(fast5_specific_output_dir):
            os.mkdir(fast5_specific_output_dir)
        mf5h = MultiFast5(multi_fast5)
        n_processed, n_error = mf5h.write_individual_fast5s(fast5_specific_output_dir)
        if delete_multi:
            os.remove(multi_fast5)
    except KeyError:
        pass
    return n_processed, n_error


def generate_individual_reads(multi_fast5_dir, out_dir):
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


def read_in_sequencing_summary(summary_path):
    """Read in the sequencing summary file from the guppy basecaller and r9.4.1 sequencing chemistry"""
    data = pd.read_csv(summary_path, sep="\t")
    return data


class MultiFast5(h5py.File):
    """Class for grabbing data from fast5 file with multiple reads. """
    __analysis__ = 'Analyses'
    __basecall_1d__ = 'Basecall_1D_000'
    __basecall_template__ = 'BaseCalled_template'

    __fastq__ = 'Fastq'
    __full_path_fastq__ = 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
    __raw__ = "Raw"
    __signal__ = 'Signal'
    __full_path_signal__ = 'Raw/Signal'

    __channel_id__ = "channel_id"
    __context_tags__ = "context_tags"
    __tracking_id__ = "tracking_id"

    def __init__(self, fname, read='r'):
        super(MultiFast5, self).__init__(fname, read)

        self.fname = fname
        self.reads = list(self.keys())

    def write_individual_fast5s(self, out_dir, write_fastq_file=False):
        """Write out individual fast5 files from the multi fast5 file

        :param out_dir: path to directory to write all the files
        """
        n_processed = 0
        n_errors = 0
        fh = None
        # option to write fastq data
        if write_fastq_file is not False:
            fh = open(write_fastq_file, "w")
        for read_id_group in self.reads:
            try:
                output_path = os.path.join(out_dir, read_id_group+".fast5")
                # Make sure everything is in the right place first
                read_id_has_basecall_data = False
                try:
                    fastq = self.get_fastq(read_id_group)
                    basecall_attrs = self.get_basecall_1d_attributes(read_id_group)
                    read_id_has_basecall_data = True
                except KeyError:
                    print("No Fastq data in {}".format(read_id_group))

                signal = self.get_signal(read_id_group)
                signal_attrs = self.get_signal_attrbutes(read_id_group)
                channel_id_attrs = self.get_channel_id_attrbutes(read_id_group)
                context_tags_attrs = self.get_context_tags_attrbutes(read_id_group)
                tracking_id_attrs = self.get_tracking_id_attrbutes(read_id_group)
                # then create a new file
                nf5h = NewFast5File(output_path)
                if read_id_has_basecall_data:
                    nf5h.write_1d_fastq(fastq, basecall_attrs)
                nf5h.write_signal(signal, signal_attrs)
                nf5h.write_unique_global_key(channel_id_attrs,
                                             context_tags_attrs,
                                             tracking_id_attrs)
                if write_fastq_file is not False:
                    if read_id_has_basecall_data:
                        fh.write(fastq.decode())

                n_processed += 1
            except KeyError:
                n_errors += 1
        #   clean up file
        if write_fastq_file is not False:
            fh.close()
        return n_processed, n_errors

    def get_basecall_1d_attributes(self, read_id_group):
        """Get basecall 1d string from read id"""
        location = self._join_path(read_id_group, self.__analysis__, self.__basecall_1d__)
        return dict(self[location].attrs)

    def get_fastq(self, read_id_group):
        """Get fastq string from read id"""
        location = self._join_path(read_id_group, self.__full_path_fastq__)
        return self[location][()]

    def get_signal(self, read_id_group):
        """Get signal from read id"""
        location = self._join_path(read_id_group, self.__full_path_signal__)
        return self[location][()]

    def get_signal_attrbutes(self, read_id_group):
        """Get signal from read id"""
        location = self._join_path(read_id_group, self.__raw__)
        return dict(self[location].attrs)

    def get_channel_id_attrbutes(self, read_id_group):
        """Get signal from read id"""
        location = self._join_path(read_id_group, self.__channel_id__)
        return dict(self[location].attrs)

    def get_context_tags_attrbutes(self, read_id_group):
        """Get signal from read id"""
        location = self._join_path(read_id_group, self.__context_tags__)
        return dict(self[location].attrs)

    def get_tracking_id_attrbutes(self, read_id_group):
        """Get signal from read id"""
        location = self._join_path(read_id_group, self.__tracking_id__)
        return dict(self[location].attrs)

    @staticmethod
    def _join_path(*args):
        return '/'.join(args)


class NewFast5File(h5py.File):
    """Class for creating and writing a new fast5 file"""
    __base_analysis__ = '/Analyses'
    __default_raw_read__ = '/Raw/Reads/Read_1'
    __default_basecall_1d__ = "Basecall_1D_000"
    __default_basecalled_template__ = "BaseCalled_template"
    __default_global_key__ = "UniqueGlobalKey"
    __channel_id__ = "channel_id"
    __context_tags__ = "context_tags"
    __tracking_id__ = "tracking_id"

    def __init__(self, fname):
        super(NewFast5File, self).__init__(fname, 'w')
        self.fname = fname
        self.analysis_grp = self.create_group(self.__base_analysis__)
        self.raw_grp = self.create_group(self.__default_raw_read__)

    def write_signal(self, signal, attributes):
        """Write out signal with attributes to new fast5 file

        :param signal: raw signal data
        :param attributes: attributes to write for signal data
        """
        self.raw_grp.create_dataset("Signal", data=signal, dtype=np.dtype('i2'))
        attrs = self.raw_grp.attrs
        for k, v in attributes.items():
            attrs[k] = v

    def write_1d_fastq(self, fastq, attributes):
        """Write fastq with attributes to new fast5 file

        :param fastq: fastq string
        :param attributes: attributes to write for fastq data
        """
        basecall_1d_group = self.analysis_grp.create_group(self.__default_basecall_1d__)
        template_data_group = basecall_1d_group.create_group(self.__default_basecalled_template__)
        dt = h5py.special_dtype(vlen=bytes)
        template_data_group.create_dataset("Fastq", data=fastq, dtype=dt)
        attrs = basecall_1d_group.attrs
        for k, v in attributes.items():
            attrs[k] = v

    def write_unique_global_key(self, channel_id, context_tags, tracking_id):
        """Write to Unique Global Key group with the attributes of channel_id, context_tags and tracking_id"""

        global_key_group = self.create_group(self.__default_global_key__)
        self.write_channel_id(global_key_group, channel_id)
        self.write_tracking_id(global_key_group, tracking_id)
        self.write_context_tags(global_key_group, context_tags)

    def write_channel_id(self, global_key_group, channel_id):
        """Write to Unique Global Key group with the attributes of channel_id"""
        channel_id_group = global_key_group.create_group(self.__channel_id__)
        attrs = channel_id_group.attrs
        for k, v in channel_id.items():
            attrs[k] = v

    def write_tracking_id(self, global_key_group, tracking_id):
        """Write to Unique Global Key group with the attributes of tracking_id"""
        tracking_id_group = global_key_group.create_group(self.__tracking_id__)
        attrs = tracking_id_group.attrs
        for k, v in tracking_id.items():
            attrs[k] = v

    def write_context_tags(self, global_key_group, context_tags):
        """Write to Unique Global Key group with the attributes of context_tags"""
        context_tags_group = global_key_group.create_group(self.__context_tags__)
        attrs = context_tags_group.attrs
        for k, v in context_tags.items():
            attrs[k] = v
