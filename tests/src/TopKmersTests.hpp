//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_TOPKMERSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_TOPKMERSTESTS_HPP_

// embed lib
#include "TopKmers.hpp"
// embed test files
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


TEST (TopKmersTests, test_generate_master_kmer_table_assignments){
  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  create_directory(tempdir);
  string outpath = tempdir.string();
  string alphabet = "ACTGlmnop";
  path assignment_file = TEST_FILES / "assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv";
  vector<string> data = {assignment_file.string()};
  path expected_log_file =  outpath / "log_file.tsv";
  string log_file = expected_log_file.string();

  path expected_output_file =  outpath / "builtAssignment.tsv";
  string out_file = expected_output_file.string();
  generate_master_kmer_table<AssignmentFile, eventkmer>(data, out_file, log_file, alphabet, 1000, 0, 2, false);
  EXPECT_EQ(lines_in_file(TOP_KMERS_ASSIGNMENT), lines_in_file(expected_output_file));
}

TEST (TopKmersTests, test_generate_master_kmer_table_alignments){
  Redirect a(true, true);
  testing::FLAGS_gtest_death_test_style="threadsafe";
  path tempdir = temp_directory_path() / "temp";
  create_directory(tempdir);
  string outpath = tempdir.string();
  string alphabet = "ACTGP";
  path alignment_file = TEST_FILES / "alignment_files/c53bec1d-8cd7-43d0-8e40-e5e363fa9fca.sm.backward.tsv";
  vector<string> data = {alignment_file.string()};
  path expected_log_file =  outpath / "log_file.tsv";
  string log_file = expected_log_file.string();

  ASSERT_THROW({(generate_master_kmer_table<AlignmentFile, FullSaEvent>)(data, outpath, log_file, alphabet, 1000, 0, 2, false); } , AssertionFailureException);
  alphabet = "ACTGE";
  path expected_output_file =  outpath / "builtAssignment.tsv";
  string out_file = expected_output_file.string();
  generate_master_kmer_table<AlignmentFile, FullSaEvent>(data, out_file, log_file, alphabet, 1000, 0, 2, false);
  EXPECT_EQ(lines_in_file(TOP_KMERS_ALIGNMENT), lines_in_file(expected_output_file));

  alphabet = "ACTGE";
  expected_output_file =  outpath / "builtAssignment.tsv";

  out_file = expected_output_file.string();
  generate_master_kmer_table<AlignmentFile, FullSaEvent>(data, out_file, log_file, alphabet, 1000, 0.5, 2, false);
  EXPECT_EQ(lines_in_file(FILTERED_TOP_KMERS_ALIGNMENT), lines_in_file(expected_output_file));
}

TEST (TopKmersTests, test_generate_master_kmer_table_wrapper){
  Redirect a(true, true);
  testing::FLAGS_gtest_death_test_style="threadsafe";
  path tempdir = temp_directory_path() / "temp";
  create_directory(tempdir);
  string outpath = tempdir.string();
  string alphabet = "ACTGlmnop";
  path assignment_file = TEST_FILES / "assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv";
  vector<string> data = {assignment_file.string()};
  path expected_output_file =  outpath / "built.backward.tsv";
  path log_path =  outpath / "log_file.tsv";
  string log_file = log_path.string();
  string out_file = expected_output_file.string();
  generate_master_kmer_table_wrapper(data, out_file, log_file, 1000, alphabet, 0, 2, false, false);
  EXPECT_EQ(lines_in_file(TOP_KMERS_ASSIGNMENT), lines_in_file(expected_output_file));

  alphabet = "ACTGP";
  path alignment_file = TEST_FILES / "alignment_files/c53bec1d-8cd7-43d0-8e40-e5e363fa9fca.sm.backward.tsv";
  data = {alignment_file.string()};
  ASSERT_THROW({ generate_master_kmer_table_wrapper(data, out_file, log_file, 1000, alphabet, 0, 2, false, false); } ,
  AssertionFailureException);
  alphabet = "ACTGE";
  generate_master_kmer_table_wrapper(data, out_file, log_file, 1000, alphabet, 0, 2, false, true);
  EXPECT_EQ(lines_in_file(TOP_KMERS_ALIGNMENT), lines_in_file(expected_output_file));
  cout << "test" << '\n';
  AlignmentFile af(expected_output_file.string());
  EXPECT_EQ(af.strand, "-");
}

TEST (AssignmentFileTests, test_iterate_assignment) {
  Redirect a(true, true);
  path input_dir = temp_directory_path() / "input" ;
  path output_dir = temp_directory_path() / "output" ;
  AssignmentFile af(ASSIGNMENT_FILE.string());
  float answer = 83.7093;
  for (auto &x: af.iterate()){
    EXPECT_FLOAT_EQ(x.descaled_event_mean, answer);
    break;
  }
}

TEST (AssignmentFileTests, test_get_k) {
  Redirect a(true, true);
  path input_dir = temp_directory_path() / "input" ;
  path output_dir = temp_directory_path() / "output" ;
  AssignmentFile af(ASSIGNMENT_FILE.string());
  int64_t answer = af.get_k();
  EXPECT_EQ(6, answer);
}

#endif //EMBED_FAST5_TESTS_SRC_TOPKMERSTESTS_HPP_
