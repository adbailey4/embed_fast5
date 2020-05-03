//
// Created by Andrew Bailey on 03/15/19.
//

// embed tests
#include "BinaryFileTests.hpp"
#include "ConcurrentQueueTests.hpp"
#include "EmbedUtilsTests.hpp"
#include "Fast5Tests.hpp"
#include "FileTests.hpp"
#include "FilterAlignmentsTests.hpp"
#include "MarginalizeVariantsTests.hpp"
#include "MaxKmersTests.hpp"
#include "TopKmersTests.hpp"
#include "VariantPathTests.hpp"
// boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>
// Standard Libray
#include <iostream>

using namespace test_files;

int main(int argc, char **argv) {
  H5Eset_auto(0, nullptr, nullptr);
  HOME = argv[1];
  cout << HOME << '\n';
  ORIGINAL_FAST5 = HOME / ORIGINAL_FAST5;
  EMPTY_FAST5 = HOME / EMPTY_FAST5;
  SIGNAL_FAST5 = HOME / SIGNAL_FAST5;
  NO_EVENT = HOME / NO_EVENT;
  NO_FAST5 = HOME / NO_FAST5;
  READ_DB = HOME / READ_DB;
  R94_FAST5 = HOME / R94_FAST5;
  R94_FASTQ = HOME / R94_FASTQ;
  R94_TEST_DIR = HOME / R94_TEST_DIR;
  POSITIONS_FILE = HOME / POSITIONS_FILE;
  ALIGNMENT_FILE = HOME / ALIGNMENT_FILE;
  ALIGNMENT_DIR = HOME / ALIGNMENT_DIR;
  CORRECT_OUTPUT = HOME / CORRECT_OUTPUT;
  ASSIGNMENT_FILE = HOME / ASSIGNMENT_FILE;
  ASSIGNMENT_DIR = HOME / ASSIGNMENT_DIR;
  ALIGNMENT_FILE_MOD = HOME / ALIGNMENT_FILE_MOD;
  TEST_FILES = HOME / TEST_FILES;
  READ_DB_DIR = HOME / READ_DB_DIR;
  TEST_POSITIONS_FILE = HOME / TEST_POSITIONS_FILE;
  RRNA_SIGNAL_FILES = HOME / RRNA_SIGNAL_FILES;
  RRNA_TEST_FILES = HOME / RRNA_TEST_FILES;
  RRNA_TEST_VARIANTS = HOME / RRNA_TEST_VARIANTS;
  DNA_TEST_VARIANTS = HOME / DNA_TEST_VARIANTS;
  TOP_KMERS_ASSIGNMENT = HOME / TOP_KMERS_ASSIGNMENT;
  TOP_KMERS_ALIGNMENT = HOME / TOP_KMERS_ALIGNMENT;
  FILTERED_TOP_KMERS_ALIGNMENT = HOME / FILTERED_TOP_KMERS_ALIGNMENT;
  AMBIG_MOD_FILE = HOME / AMBIG_MOD_FILE;
  ::testing::InitGoogleMock(&argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
