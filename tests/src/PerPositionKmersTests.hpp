//
// Created by Andrew Bailey on 5/3/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_

// embed source
#include "PerPositonKmers.hpp"
#include "EmbedUtils.hpp"
#include "TestFiles.hpp"

// boost
#include <boost/filesystem.hpp>

// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace boost::filesystem;
using namespace std;
using namespace embed_utils;
using namespace test_files;

TEST (PerPositionKmersTests, test_process_alignment) {
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "test.event";
  uint64_t num_locks = 1000;
  ReferenceHandler reference(PUC_REFERENCE.string());
  PerPositonKmers ppk(num_locks, test_file, reference);
  path alignment_file = PUC_5MER_ALIGNMENTS/"03274a9a-0eab-422e-ace7-b35fd3a0f48c.sm.forward.tsv";
  AlignmentFile af(alignment_file.string());
  ppk.process_alignment(af);
  EXPECT_EQ("hello!", "hello!");
}


#endif //EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
