//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_VARIANTPATHTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_VARIANTPATHTESTS_HPP_


// embed lib
#include "VariantPath.hpp"
#include "LoadVariantPaths.hpp"
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


TEST (VariantPathTests, test_load_positions_file){
  Redirect a(true, true);
  VariantPath vp;
  vp.load_positions_file(TEST_POSITIONS_FILE.string());
  vector<uint64_t> path_multiplier_answer = {1, 3, 12};
  ASSERT_THAT(path_multiplier_answer, ElementsAreArray(vp.path_multiplier["gi_ecoli+"]));
  EXPECT_EQ(3, vp.num_positions["gi_ecoli+"]);
  EXPECT_EQ(48, vp.num_ids["gi_ecoli+"]);
  EXPECT_EQ("LP", vp.all_variant_chars);
  EXPECT_EQ(0, vp.position_to_path_index["gi_ecoli+"][419]);
  EXPECT_EQ(1, vp.position_to_path_index["gi_ecoli+"][507]);
  EXPECT_EQ(2, vp.position_to_path_index["gi_ecoli+"][1089]);
}

TEST (VariantPathTests, test_id_2_path_2_id){
  Redirect a(true, true);

  VariantPath vp(TEST_POSITIONS_FILE.string());
  uint64_t n_ids = vp.num_ids["gi_ecoli+"];
  for (uint64_t i = 0; i < n_ids; ++i){
//    cout << i << " (";
//    vector<uint64_t> path = vp.id_to_path("gi_ecoli+", i);
//    for (auto &j: path){
//      cout << j;
//    }
//    cout << ")\n";
    EXPECT_EQ(i, vp.path_to_id("gi_ecoli+", vp.id_to_path("gi_ecoli+", i)));
  }
}

TEST (VariantPathTests, test_variant_call_to_path_to_id){
  Redirect a(true, true);
//  create some variant calls
  vector<VariantCall> some_calls(3);
  VariantCall vc("gi_ecoli", "+", 419, "CE");
  vc.normalized_probs = {0.2, 0.8};
  some_calls[0] = vc;
  VariantCall vc1("gi_ecoli", "+", 507, "CEO");
  vc1.normalized_probs = {0.2, 0.7, 0.1};
  some_calls[1] = vc1;
  VariantCall vc2("gi_ecoli", "+", 1089, "CEO");
  vc2.normalized_probs = {0.2, 0.2, 0.8};
  some_calls[2] = vc2;
//  tests
  VariantPath vp(TEST_POSITIONS_FILE.string());
  vector<uint64_t> id_answer = {2, 2, 3};
  ASSERT_THAT(id_answer, ElementsAreArray(vp.variant_call_to_path("gi_ecoli+", some_calls)));
  EXPECT_EQ(44, vp.variant_call_to_id("gi_ecoli+", some_calls));
}

TEST (VariantPathTests, test_path_to_bases){
  Redirect a(true, true);
  VariantPath vp(TEST_POSITIONS_FILE.string());
  vector<uint64_t> id_answer = {1, 1, 1};
  EXPECT_EQ("CCC", vp.path_to_bases("gi_ecoli+", id_answer));
  id_answer = {0, 0, 0};
  EXPECT_EQ("---", vp.path_to_bases("gi_ecoli+", id_answer));
  id_answer = {1, 2, 2};
  EXPECT_EQ("CEE", vp.path_to_bases("gi_ecoli+", id_answer));
  id_answer = {0, 1, 3};
  EXPECT_EQ("-CO", vp.path_to_bases("gi_ecoli+", id_answer));
  id_answer = {0, 1, 3, 3};
  ASSERT_THROW(vp.path_to_bases("gi_ecoli+", id_answer), AssertionFailureException);
}

TEST (VariantPathTests, test_load_variants){
  Redirect a(true, true);
  path positions_file = RRNA_TEST_FILES/"/16S_final_branch_points.positions";
  LoadVariantPaths lvp(positions_file.string(), RRNA_SIGNAL_FILES.string(), true, 10);
  path tempdir = temp_directory_path() / "temp";
  create_directory(tempdir);
  path output_per_read = tempdir/"per_read_calls.tsv";
  path correct_per_read = RRNA_TEST_FILES/"test_output_dir/per_read_calls.tsv";
  lvp.write_per_read_calls(output_per_read.string());
  EXPECT_EQ(number_of_columns(correct_per_read), number_of_columns(output_per_read));

  path output_per_path = tempdir/"per_path_counts.tsv";
  lvp.write_per_path_counts(output_per_path.string());
  path correct_per_path = RRNA_TEST_FILES/"test_output_dir/per_path_counts.tsv";
  EXPECT_TRUE(compare_files(correct_per_path, output_per_path));
}

#endif //EMBED_FAST5_TESTS_SRC_VARIANTPATHTESTS_HPP_
