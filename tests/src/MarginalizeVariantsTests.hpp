//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_MARGINALIZEVARIANTSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_MARGINALIZEVARIANTSTESTS_HPP_

// embed lib
#include "MarginalizeVariants.hpp"
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

TEST (MarginalizeVariantsTests, test_load_variants) {
  Redirect a(true, true);
  vector<VariantCall> data;
  VariantCall vc1("test", "+", 10, "AT");
  vc1.normalized_probs = {0.1, 0.9};
  VariantCall vc2("test", "+", 14, "GT");
  vc2.normalized_probs = {0.1, 0.9};
  VariantCall vc3("test", "+", 10, "AT");
  vc3.normalized_probs = {0.9, 0.1};

  data.push_back(vc1);
  data.push_back(vc2);
  data.push_back(vc3);

  MarginalizeVariants mf;
  mf.load_variants(&data);
  EXPECT_EQ(2, mf.per_genomic_position["test"].first.size());
  EXPECT_EQ(2, mf.per_genomic_position["test"].first[10].coverage);
  EXPECT_EQ(2, mf.per_genomic_position["test"].first[10].coverage);
  EXPECT_EQ(1, mf.per_genomic_position["test"].first[14].coverage);
  EXPECT_EQ(0, mf.per_genomic_position["test"].first[23].coverage);
  EXPECT_EQ(-1, mf.per_genomic_position["test"].first[23].start);
  EXPECT_EQ(1, mf.per_genomic_position["test"].first[10].hits[0]);
  EXPECT_EQ(1, mf.per_genomic_position["test"].first[10].hits[1]);
  EXPECT_EQ(1, mf.per_genomic_position["test"].first[14].hits[1]);
  EXPECT_EQ(0, mf.per_genomic_position["test"].first[14].hits[0]);

}

TEST (MarginalizeVariantsTests, test_write_to_file) {
  Redirect a(true, true);
  vector<VariantCall> data;
  VariantCall vc1("test", "+", 10, "AT");
  vc1.normalized_probs = {0.1, 0.9};
  VariantCall vc2("test", "+", 14, "GT");
  vc2.normalized_probs = {0.1, 0.9};
  VariantCall vc3("test", "+", 10, "AT");
  vc3.normalized_probs = {0.9, 0.1};

  data.push_back(vc1);
  data.push_back(vc2);
  data.push_back(vc3);

  MarginalizeVariants mf;
  mf.load_variants(&data);
  path bed_file = TEST_FILES / "bed_files/test.bed";
  path bed_file2 = temp_directory_path() / "test2.bed";
  mf.write_to_file(bed_file2);
  EXPECT_TRUE(compare_files(bed_file, bed_file2));
}

#endif //EMBED_FAST5_TESTS_SRC_MARGINALIZEVARIANTSTESTS_HPP_
