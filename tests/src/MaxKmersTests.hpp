//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_MAXKMERSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_MAXKMERSTESTS_HPP_

// embed lib
#include "MaxKmers.hpp"
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

TEST (MaxKmersTests, test_max_kmer_initialization) {
  Redirect a(true, true);
  ASSERT_THROW(MaxKmers<eventkmer> mk(10, "ATGC", 0, 0), std::exception);
}

TEST (MaxKmersTests, test_create_kmers) {
  Redirect a(true, true);
  MaxKmers<eventkmer> mk(10, "ATGC", 5, 0);
  string alphabet = "ATGC";
  vector<string> kmers = mk.create_kmers(alphabet, 5);
  EXPECT_EQ(kmers.size(), 1024);
  kmers = mk.create_kmers(alphabet, 6);
  EXPECT_EQ(kmers.size(), 4096);
}

TEST (MaxKmersTests, test_get_kmer_index) {
  Redirect a(true, true);
  MaxKmers<eventkmer> mk(10, "ATGC", 5, 0);
  string poly_a = "AAAAA";
  int poly_a_index = mk.get_kmer_index(poly_a);
  EXPECT_EQ(poly_a_index, 0);
  string poly_t = "TTTTT";
  int poly_t_index = mk.get_kmer_index(poly_t);
  EXPECT_EQ(poly_t_index, 1023);
}

TEST (MaxKmersTests, test_get_index_kmer) {
  Redirect a(true, true);
  MaxKmers<eventkmer> mk(10, "ATGC", 5, 0);
  size_t index = 0;
  string poly_a = mk.get_index_kmer(index);
  EXPECT_EQ(poly_a, "AAAAA");
  index = 1023;
  string poly_t = mk.get_index_kmer(index);
  EXPECT_EQ(poly_t, "TTTTT");
}

TEST (MaxKmersTests, test_add_to_heap) {
  Redirect a(true, true);
  MaxKmers<eventkmer> mk(10, "ATGC", 5, 0);
  eventkmer my_kmer = eventkmer("AAAAA", 10, "t", 1.2);
  mk.add_to_heap(my_kmer);
  eventkmer my_kmer2 = eventkmer("AAAAA", 10, "t", 1.4);
#pragma omp parallel for shared(mk, my_kmer2) default(none)
  for (int i = 0; i < 10; i++){
    mk.add_to_heap(my_kmer2);
  }
  auto match = mk.kmer_queues[0].top();
  EXPECT_EQ(match.path_kmer, my_kmer2.path_kmer);
  EXPECT_EQ(match.strand, my_kmer2.strand);
  EXPECT_EQ(match.descaled_event_mean, my_kmer2.descaled_event_mean);
  EXPECT_FLOAT_EQ(match.posterior_probability, my_kmer2.posterior_probability);

  MaxKmers<FullSaEvent> mk2(10, "ATGC", 5, 0);
  FullSaEvent my_kmer3("a", 1, "string reference_kmer", "string read_file", "t",
                       10, 20, 20, 20, "string aligned_kmer",
                       20, 20, 1.2, 10,
                       3, "AAAAA");
  mk2.add_to_heap(my_kmer3);
  FullSaEvent my_kmer4("a", 1, "string reference_kmer", "string read_file", "t",
                       10, 20, 20, 20, "string aligned_kmer",
                       20, 20, 1.4, 10,
                       3, "AAAAA");
#pragma omp parallel for shared(mk2, my_kmer4) default(none)
  for (int i = 0; i < 10; i++){
    mk2.add_to_heap(my_kmer4);
  }
  auto match2 = mk.kmer_queues[0].top();
  EXPECT_EQ(match2.path_kmer, my_kmer4.path_kmer);
  EXPECT_EQ(match2.strand, my_kmer4.strand);
  EXPECT_EQ(match2.descaled_event_mean, my_kmer4.descaled_event_mean);
  EXPECT_FLOAT_EQ(match2.posterior_probability, my_kmer4.posterior_probability);
}

TEST (MaxKmersTests, test_write_to_file) {
  Redirect a(true, true);
  MaxKmers<eventkmer> mk(10, "ATGC", 5, 0);
  for (int i = 10; i > 0; i--){
    eventkmer my_kmer = eventkmer("AAAAA", 10, "t", i);
    mk.add_to_heap(my_kmer);
  }
  if (exists(NO_FAST5)){
    remove(NO_FAST5);
  }
  mk.write_to_file(NO_FAST5, true);
  EXPECT_TRUE(is_regular_file(NO_FAST5));
  remove(NO_FAST5);

  path log_file = TEST_FILES / "log_file.txt";

  mk.write_to_file(NO_FAST5, log_file, true);
  EXPECT_TRUE(is_regular_file(NO_FAST5));
  EXPECT_TRUE(is_regular_file(log_file));
  remove(NO_FAST5);
  remove(log_file);
}

#endif //EMBED_FAST5_TESTS_SRC_MAXKMERSTESTS_HPP_
