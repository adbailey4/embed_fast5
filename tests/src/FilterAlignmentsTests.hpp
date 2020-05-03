//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_FILTERALIGNMENTSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_FILTERALIGNMENTSTESTS_HPP_

// embed lib
#include "FilterAlignments.hpp"
#include "EmbedUtils.hpp"
//embed test files
#include "TestFiles.hpp"
// boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace test_files;
using namespace embed_utils;
using namespace boost::filesystem;
using ::testing::ElementsAreArray;

TEST (FilterAlignmentTests, test_filter_alignment_files) {
  Redirect a(true, false);
  path input_dir = temp_directory_path() / "input" ;
  path output_dir = temp_directory_path() / "output" ;
  copyDir(ALIGNMENT_DIR, input_dir);
  string empty;
  string in_dir(input_dir.string());
  const string& pos_file(POSITIONS_FILE.string());
  string out_dir(output_dir.string());
  //    omp_set_num_threads(2); // Use 2 threads
  filter_alignment_files(in_dir, pos_file, out_dir, empty);

  directory_iterator end_itr;
//    Get all tsvs to process
  for (directory_iterator itr(output_dir); itr != end_itr; ++itr) {

    EXPECT_TRUE(compare_files(itr->path(), CORRECT_OUTPUT / itr->path().filename()));

  }
  remove_all(input_dir);
  remove_all(output_dir);
}



#endif //EMBED_FAST5_TESTS_SRC_FILTERALIGNMENTSTESTS_HPP_
