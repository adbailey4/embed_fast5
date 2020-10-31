//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_FILETESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_FILETESTS_HPP_

// embed lib
#include "TopKmers.hpp"
// embed test files
#include "TestFiles.hpp"
// boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>
//std lib
#include <numeric>

using namespace test_files;
using namespace embed_utils;
using namespace boost::filesystem;


TEST (PositionsFileTests, test_load) {
  Redirect a(true, true);
  PositionsFile pf(POSITIONS_FILE.string(), 5);
  string contig = "gi_ecoli+";
  EXPECT_TRUE(contains(pf.m_data[contig], 419));
}

TEST (PositionsFileTests, test_iterate) {
  Redirect a(true, true);
  PositionsFile pf(POSITIONS_FILE.string());
  for (auto &position: pf.iterate()) {
    EXPECT_EQ("gi_ecoli", position.contig);
  }
  PositionsFile pf2(POSITIONS_FILE.string());

  positions_coro::pull_type data = pf2.iterate();
  PositionLine something = data.get();
  PositionLine true_something("gi_ecoli", 419, "+", "C", "E");
  EXPECT_EQ(true_something, something);
}

TEST (PositionsFileTests, test_load_positions_map) {
  Redirect a(true, true);
  PositionsFile pf(POSITIONS_FILE.string());
  pf.load_positions_map();
  EXPECT_EQ(make_tuple("C", "E"), pf.positions_map["gi_ecoli+"][419]);
}

TEST (PositionsFileTests, test_is_in) {
  Redirect a(true, true);
  PositionsFile pf(POSITIONS_FILE.string(), 5);
  string contig = "gi_ecoli+";
  EXPECT_TRUE(pf.is_in(contig, 419));
  EXPECT_TRUE(pf.is_in(contig, 415));
  EXPECT_FALSE(pf.is_in(contig, 414));
}

TEST (AlignmentFileTests, test_get_strand) {
  Redirect a(true, true);
  AlignmentFile af(ALIGNMENT_FILE.string());
  af.get_strand();
  EXPECT_EQ("-", af.strand);
}

TEST (AlignmentFileTests, test_iterate) {
  Redirect a(true, true);
  AlignmentFile af(ALIGNMENT_FILE.string());
  for (auto &event: af.iterate()){
    EXPECT_EQ("gi_ecoli", event.contig);
  }
}

TEST (AlignmentFileTests, test_filter_by_ref_bases) {
  Redirect a(true, true);
  AlignmentFile af(ALIGNMENT_FILE.string());
  string bases = "G";
  for (auto &event: af.filter_by_ref_bases(bases)){
    EXPECT_TRUE(event.reference_kmer.find("G")!=string::npos);
  }
}

TEST (AlignmentFileTests, test_get_variant_calls) {
  Redirect a(true, true);
  AlignmentFile af(ALIGNMENT_FILE_MOD.string(), true);
  string bases = "f";
  std::map<string, string> ambig_bases = create_ambig_bases();
  vector<VariantCall> data = af.get_variant_calls(bases, &ambig_bases);
  for (auto i: data){
    EXPECT_EQ("+", i.strand);
    EXPECT_EQ("rna_fake", i.contig);
    EXPECT_EQ("AF", i.bases);
    EXPECT_FLOAT_EQ(1.0, std::accumulate(i.normalized_probs.begin(),
        i.normalized_probs.end(), 0.0));
  }
  AlignmentFile af2(RRNA_TEST_VARIANTS.string(), true);
  string bases2 = "Y";
  vector<VariantCall> data2 = af2.get_variant_calls(bases2, &ambig_bases);
  for (auto i: data2){
    EXPECT_EQ("+", i.strand);
    EXPECT_EQ("ecoli_MRE600", i.contig);
    EXPECT_EQ("CT", i.bases);
    EXPECT_FLOAT_EQ(1.0, std::accumulate(i.normalized_probs.begin(),
        i.normalized_probs.end(), 0.0));
  }
  bases2 = "K";
  vector<VariantCall> data3 = af2.get_variant_calls(bases2, &ambig_bases);
  for (auto i: data3){
    EXPECT_EQ("+", i.strand);
    EXPECT_EQ("ecoli_MRE600", i.contig);
    EXPECT_EQ("GT", i.bases);
    EXPECT_FLOAT_EQ(1.0, std::accumulate(i.normalized_probs.begin(),
        i.normalized_probs.end(), 0.0));
  }

  AlignmentFile af3(DNA_TEST_VARIANTS.string(), false);
  bases2 = "P";
  vector<VariantCall> data4 = af3.get_variant_calls(bases2, &ambig_bases);
  for (auto i: data4){
    EXPECT_EQ("-", i.strand);
    EXPECT_EQ("gi_ecoli", i.contig);
    EXPECT_EQ("CE", i.bases);
    float number = std::accumulate(i.normalized_probs.begin(),
                                   i.normalized_probs.end(), 0.0);
    if (isnan(number)){
      cout << "nan?" << "\n";
    }
    EXPECT_FLOAT_EQ(1.0, std::accumulate(i.normalized_probs.begin(),
        i.normalized_probs.end(), 0.0));
  }
}

TEST (AlignmentFileTests, test_filter) {
  Redirect a(true, true);
  path tempdir = temp_directory_path()/ "temp";
  create_directory(tempdir);
  path test_output = tempdir / "test_output.assignment.tsv";
  AlignmentFile af(ALIGNMENT_FILE.string());
  PositionsFile pf(POSITIONS_FILE.string(), 5);
  af.filter_by_positions(&pf, test_output, "");
  std::ifstream in_file(test_output.c_str());
  if (in_file.good()) {
    // read the file
    std::string line;
    while (getline(in_file, line)) {
      std::vector<std::string> fields = embed_utils::split_string(line, '\t');
      std::size_t found = fields[0].find('C');
      EXPECT_TRUE(found!=std::string::npos);
    }
  }
  in_file.close();
  remove_all(tempdir);
}

#endif //EMBED_FAST5_TESTS_SRC_FILETESTS_HPP_
