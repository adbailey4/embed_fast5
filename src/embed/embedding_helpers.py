import os
import pysam
import numpy as np
from contextlib import closing
import subprocess
from Bio import SeqIO
from embed.fast5 import Fast5
from py3helpers.utils import get_all_sub_directories
from py3helpers.multiprocess import *


def parse_readdb(readdb, directories=None):
    """Parse readdb file

    :param readdb: path to readdb file
    :param directories: path to directories of where the reads are
    :return yields read_id, fast5_path
    """
    assert readdb.endswith("readdb"), "readdb file must end with .readdb: {}".format(readdb)
    with open(readdb, 'r') as fh:
        while True:
            # read line
            line = fh.readline()
            if not line:
                break
            split_line = line.split()
            if len(split_line) == 2:
                if directories:
                    for dir_path in directories:
                        full_path = os.path.join(dir_path, split_line[1])
                        if os.path.exists(full_path):
                            yield split_line[0], full_path
                else:
                    full_path = os.path.abspath(split_line[1])
                    if os.path.exists(full_path):
                        yield split_line[0], full_path


def call_nanopolish_index(fast5_dir, fastq):
    """Call nanopolish index on some files"""
    fast5_dir = os.path.abspath(fast5_dir)
    fastq = os.path.abspath(fastq)
    nanopolish_command = "embed_main index -d {fast5_dir} {fastq}".format(fast5_dir=fast5_dir, fastq=fastq)
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


def call_embed_main(fastq):
    """Call embed on all files"""
    embed_command = "embed_main embed -r {fastq}".format(fastq=fastq)
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


def edit_fastq_input(fastq_file, out_fastq):
    command = "sed s/strand/strand.fast5/ {}".format(fastq_file)
    with open(out_fastq, 'w') as fh:
        try:
            command = command.split()
            proc = subprocess.Popen(command, stdout=fh, stderr=subprocess.PIPE)
            output, errors = proc.communicate()
            errors = errors.decode().splitlines()
            for x in errors:
                print(x)
        except Exception as e:
            print("[edit_fastq_input] exception ({}) running edit_fastq_input: {}".format(type(e), e))
            raise e
    return out_fastq


def check_fast5_extension_on_fastq(fastq_file):
    """Check if fastq file is formatted correctly"""
    record = next(SeqIO.parse(fastq_file, "fastq"))
    if record.description.endswith("strand"):
        return False
    return True


def embed_fastqs_into_fast5s(fastq_file, fast5_dirs, fastq_readdb=None):
    """Embed fastqs into fast5 file using readdb file
    :param fastq_file: path to fastq file
    :param fast5_dirs: list of fast5 directories
    :param fastq_readdb: path to readdb. If not set it assumes the readdb is in fastq directory
    :return: number of reads processed
    """
    if fastq_readdb is None:
        fastq_readdb = fastq_file+".index.readdb"
    record_dict = SeqIO.index(fastq_file, "fastq")
    counter = 0
    for read_id, fast5_path in parse_readdb(fastq_readdb, fast5_dirs):
        fastq_record = record_dict[read_id]
        f5_handle = Fast5(fast5_path, 'r+')
        f5_handle.set_fastq("Analyses/Basecall_1D_000", fastq_record.format("fastq").rstrip().encode('ascii'),
                            overwrite=True)
        f5_handle.close()
        counter += 1
    record_dict.close()
    return counter


def multiprocess_embed_fastqs_into_fast5s(fastq_file, fast5_dirs, fastq_readdb=None, worker_count=2):
    """Multiprocess embed fastqs into fast5 file using readdb file
    :param fastq_file: path to fastq file
    :param fast5_dirs: list of fast5 directories
    :param fastq_readdb: path to readdb. If not set it assumes the readdb is in fastq directory
    :return: number of reads processed
    """
    if fastq_readdb is None:
        fastq_readdb = fastq_file+".index.readdb"
    record_dict = SeqIO.index(fastq_file, "fastq")
    all_data = []
    for read_id, fast5_path in parse_readdb(fastq_readdb, fast5_dirs):
        fastq_string = record_dict[read_id].format("fastq").rstrip().encode('ascii')
        all_data.append((fast5_path, fastq_string))
    record_dict.close()

    service = BasicService(embed_fastqs_wrapper, service_name="multiprocess_embed_fastqs")
    total, failure, messages, output = run_service(service.run, all_data,
                                                   {}, ["fast5_path", "fastq_string"], worker_count=worker_count)
    return sum(output)


def embed_fastqs_wrapper(fast5_path, fastq_string):
    """Multiprocess embed fastqs into fast5 file using readdb file
    :param fast5_path: path to fast5 file
    :param fastq_string: fastq in string format
    :return: number of reads processed
    """
    f5_handle = Fast5(fast5_path, 'r+')
    f5_handle.set_fastq("Analyses/Basecall_1D_000", fastq_string,
                        overwrite=True)
    f5_handle.close()
    return True


def parse_seq_summary(seq_summary, directories):
    """Parse seq_summary file

    :param seq_summary: path to seq_summary file
    :param directories: path to directories of where the reads are
    """
    assert seq_summary.endswith("tsv"), "seq_summary file must end with .tsv: {}".format(seq_summary)
    with open(seq_summary, 'r') as fh:
        for line in fh:
            split_line = line.split()
            for dir_path in directories:
                full_path = os.path.join(dir_path, split_line[0])
                if os.path.exists(full_path):
                    yield split_line[1], full_path


def parse_read_name_map_file(read_map, directories, recursive=False):
    """Parse either a seq summary file or a readdb file
    :param read_map: either a readdb file or sequencing summary file
    :param directories: check all the directories for the fast5 path
    :param recursive: boolean option to check the sub directories of input directories
    """
    if read_map.endswith("readdb"):
        name_index = 0
        path_index = 1
    else:
        name_index = 1
        path_index = 0
    for dir_path in directories:
        assert os.path.isdir(dir_path), "Path provided does not exist or is not a directory: {}".format(dir_path)
    with open(read_map, 'r') as fh:
        for line in fh:
            split_line = line.split()
            if len(split_line) == 2:
                for dir_path in directories:
                    if recursive:
                        directories2 = get_all_sub_directories(dir_path)
                        for dir_path2 in directories2:
                            full_path = os.path.join(dir_path2, split_line[path_index])
                            if os.path.exists(full_path):
                                yield split_line[name_index], os.path.abspath(full_path)
                    else:
                        full_path = os.path.join(dir_path, split_line[path_index])
                        if os.path.exists(full_path):
                            yield split_line[name_index], os.path.abspath(full_path)


def filter_reads(alignment_file, readdb, read_dirs, quality_threshold=7, recursive=False, trim=False):
    """Filter fast5 files based on a quality threshold and if there is an alignment
    :param alignment_file: bam aligment file
    :param readdb: readdb or sequence summary file
    :param read_dirs: list of directories
    :param quality_threshold: phred quality score min threshold for passing
    :param recursive: search directories recursively for more fast5 dirs
    :param trim: number of bases to analyze
    """
    assert alignment_file.endswith("bam"), "Alignment file must be in BAM format: {}".format(alignment_file)
    # grab aligned segment
    if trim:
        assert isinstance(trim, int), "Trim needs to be an integer: {}".format(trim)
    else:
        trim = np.inf
    n_bases = 0
    n_files = 0
    with closing(pysam.AlignmentFile(alignment_file, 'rb')) as bamfile:
        name_indexed = pysam.IndexedReads(bamfile)
        name_indexed.build()
        for name, fast5 in parse_read_name_map_file(readdb, read_dirs, recursive=recursive):
            try:
                if trim < n_bases:
                    print("Filtered {} files for {} bases".format(n_files, n_bases))
                    break
                iterator = name_indexed.find(name)
                for aligned_segment in iterator:
                    if aligned_segment.is_secondary or aligned_segment.is_unmapped \
                            or aligned_segment.is_supplementary or aligned_segment.has_tag("SA"):
                        continue
                    # get data and sanity check
                    if aligned_segment.query_qualities is not None:
                        if np.mean(aligned_segment.query_qualities) < quality_threshold:
                            continue
                    n_files += 1
                    n_bases += aligned_segment.query_length
                    yield fast5, aligned_segment
            except KeyError:
                print("Found no alignments for {}".format(fast5))


def create_new_directories_for_filter_reads(in_dir, out_dir):
    """Copy directory structure and return the input directory and output directory for each interal dir
    :param in_dir: top input directory to
    :param out_dir: top of the output directories
    """
    for sub_in_dir in get_all_sub_directories(in_dir):
        # make dir if it doesnt exist
        sub_out_dir = os.path.join(out_dir, os.path.basename(sub_in_dir))
        if not os.path.isdir(sub_out_dir):
            os.mkdir(sub_out_dir)
        yield sub_in_dir, sub_out_dir


def find_fast5s_from_ids_readdb(readdb, read_ids, read_dirs, recursive=False):
    """Find the corresponding fast5 files given readids"""
    for name, fast5 in parse_read_name_map_file(readdb, read_dirs, recursive=recursive):
        if name.split("_")[0] in read_ids:
            yield name, fast5
