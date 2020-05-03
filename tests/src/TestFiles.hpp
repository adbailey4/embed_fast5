//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_TESTFILES_HPP_
#define EMBED_FAST5_TESTS_SRC_TESTFILES_HPP_

// boost
#include <boost/filesystem.hpp>
using namespace boost::filesystem;

// Keep track of a bunch of paths
#define ORIGINAL_FAST51 "/tests/test_files/fast5s/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_mux_scan_MA_821_R9_4_NA12878_11_17_16_95723_ch458_read26_strand.fast5"
#define EMPTY_FAST51 "/tests/test_files/fast5s/empty_tester.fast5"
#define SIGNAL_FAST51 "/tests/test_files/fast5s/just_signal.fast5"
#define NO_FAST51 "/tests/test_files/fast5s/non_existent.fast5"
#define NO_EVENT1 "/tests/test_files/new_individual_files/read_002f9702-c19e-48c2-8e72-9021adbd4a48.fast5"
#define READ_DB1 "/tests/test_files/new_individual_files/new_individual.fastq"
#define R94_FAST51 "tests/test_files/r94_tests/fast5"
#define R94_FASTQ1 "tests/test_files/r94_tests/small_sample.fastq"
#define R94_TEST_DIR1 "tests/test_files/r94_tests"
#define POSITIONS_FILE1 "tests/test_files/positions_tests/CCWGG_ecoli_k12_mg1655.positions"
#define TEST_POSITIONS_FILE2 "tests/test_files/positions_tests/short_test.positions"
#define ALIGNMENT_FILE1 "tests/test_files/positions_tests/5cc86bac-79fd-4897-8631-8f1c55954a45.sm.backward.tsv"
#define ALIGNMENT_DIR1 "tests/test_files/positions_tests/"
#define CORRECT_OUTPUT1 "tests/test_files/positions_tests/correct_outputs/"
#define ASSIGNMENT_FILE1 "tests/test_files/assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv"
#define ASSIGNMENT_DIR1 "tests/test_files/assignment_files"
#define ALIGNMENT_FILE_MOD1 "tests/test_files/alignment_files/7d31de25-8c15-46d8-a08c-3d5043258c89.sm.forward.tsv"
#define TEST_FILES1 "tests/test_files/"
#define READ_DB_DIR1 "/tests/test_files/new_individual_files"
#define RRNA_SIGNAL_FILES1 "tests/test_files/rRNA_test_files/rRNA_signal_files"
#define RRNA_TEST_FILES1 "tests/test_files/rRNA_test_files"
#define RRNA_TEST_VARIANTS1 "tests/test_files/rRNA_test_files/rRNA_signal_files/1a424c03-7f13-4565-b991-3f7d9c36bf86.sm.forward.tsv"
#define DNA_TEST_VARIANTS1 "tests/test_files/alignment_files/c53bec1d-8cd7-43d0-8e40-e5e363fa9fca.sm.backward.tsv"
#define TOP_KMERS_ASSIGNMENT1 "tests/test_files/test_top_kmers/builtAssignment2.tsv"
#define TOP_KMERS_ALIGNMENT1 "tests/test_files/test_top_kmers/builtAssignment1.tsv"
#define TOP_KMERS_ALIGNMENT2 "tests/test_files/test_top_kmers/builtAssignment3.tsv"
#define AMBIG_MOD_FILE1 "tests/test_files/test_ambig_model/mod_variants.model"

namespace test_files {
path HOME = "This is not a path";
path ORIGINAL_FAST5 = ORIGINAL_FAST51;
path EMPTY_FAST5 = EMPTY_FAST51;
path SIGNAL_FAST5 = SIGNAL_FAST51;
path NO_EVENT = NO_EVENT1;
path NO_FAST5 = NO_FAST51;
path READ_DB = READ_DB1;
path R94_FAST5 = R94_FAST51;
path R94_FASTQ = R94_FASTQ1;
path R94_TEST_DIR = R94_TEST_DIR1;
path POSITIONS_FILE = POSITIONS_FILE1;
path ALIGNMENT_FILE = ALIGNMENT_FILE1;
path ALIGNMENT_FILE_MOD = ALIGNMENT_FILE_MOD1;
path ALIGNMENT_DIR = ALIGNMENT_DIR1;
path CORRECT_OUTPUT = CORRECT_OUTPUT1;
path ASSIGNMENT_FILE = ASSIGNMENT_FILE1;
path TEST_FILES = TEST_FILES1;
path ASSIGNMENT_DIR = ASSIGNMENT_DIR1;
path READ_DB_DIR = READ_DB_DIR1;
path TEST_POSITIONS_FILE = TEST_POSITIONS_FILE2;
path RRNA_SIGNAL_FILES = RRNA_SIGNAL_FILES1;
path RRNA_TEST_FILES = RRNA_TEST_FILES1;
path RRNA_TEST_VARIANTS = RRNA_TEST_VARIANTS1;
path DNA_TEST_VARIANTS = DNA_TEST_VARIANTS1;
path TOP_KMERS_ASSIGNMENT = TOP_KMERS_ASSIGNMENT1;
path TOP_KMERS_ALIGNMENT = TOP_KMERS_ALIGNMENT1;
path FILTERED_TOP_KMERS_ALIGNMENT = TOP_KMERS_ALIGNMENT2;
path AMBIG_MOD_FILE = AMBIG_MOD_FILE1;
}

#endif //EMBED_FAST5_TESTS_SRC_TESTFILES_HPP_
