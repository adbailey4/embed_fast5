//
// Created by Andrew Bailey on 5/23/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_BASEKMERTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_BASEKMERTESTS_HPP_

// embed source
#include "BaseKmer.hpp"
#include "TestFiles.hpp"
//boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace boost::filesystem;
using namespace std;
using namespace embed_utils;
using namespace test_files;


TEST (BaseKmerTests, test_Event) {
  Redirect a(true, true);
  Event e(1, 2);
  EXPECT_EQ(1, e.descaled_event_mean);
  EXPECT_EQ(2, e.posterior_probability);
}

TEST (BaseKmerTests, test_PosKmer) {
  Redirect a(true, true);
  Event e(1, 2);
  Event e2(1.1, 2.2);
  PosKmer k("ATGCC", 1);
  k.add_event(e);
  k.add_event(e2);
  k.add_event(1.02, 2.1);
  EXPECT_EQ("ATGCC", k.kmer);
  EXPECT_FLOAT_EQ(2.2, k.events.front().posterior_probability);
  EXPECT_EQ(1, k.num_events());
  PosKmer k2("ATGCC", 2);
  k2.add_event(e);
  k2.add_event(e2);
  k2.add_event(1.02, 2.1);
  EXPECT_EQ("ATGCC", k2.kmer);
  EXPECT_FLOAT_EQ(2.2, k2.events.back().posterior_probability);
  EXPECT_EQ(2, k2.num_events());
}

TEST (BaseKmerTests, test_PosKmer_kde) {
  Redirect a(true, true);
  PosKmer k("ATGCC");
  k.add_event(1, 1);
  k.add_event(2, 0.5);
  k.add_event(3, 1);
  k.add_event(4, 0.5);
  k.add_event(5, 1);
  k.add_event(6, 0.5);
  k.add_event(7, 1);
  k.add_event(8, 0.5);
  k.add_event(9, 1);
  k.add_event(10, 1);
  EXPECT_EQ(10, k.num_events());
  ASSERT_THROW(k.get_hist(0,10,10), AssertionFailureException);
  vector<uint64_t> hist{1,1,1,1,1,1,1,1,1,1};
  ASSERT_THAT(hist, ElementsAreArray(k.get_hist(0,10.1,10)));
  ASSERT_THAT(hist, ElementsAreArray(k.get_hist(0,10.1,10, 0.5)));
  hist = {1,0,1,0,1,0,1,0,1,1};
  ASSERT_THAT(hist, ElementsAreArray(k.get_hist(0,10.1,10, 0.51)));
  hist = {0,1,1,1,1,1,1,1,1,1,1};
  ASSERT_THAT(hist, ElementsAreArray(k.get_hist(0,11,11)));
}

TEST (BaseKmerTests, test_Position) {
  Redirect a(true, true);
  string kmer = "ATGCC";
  Event e(1, 2);
  Event e2(1.1, 2.2);
  Event e3 = e;
  shared_ptr<PosKmer> k = make_shared<PosKmer>(kmer, 1);
  k->add_event(e);
  k->add_event(e2);
  shared_ptr<PosKmer> k2 = k;
  Position p(1);
  p.add_kmer(move(k));
  ASSERT_THROW({(p.add_kmer)(k2);}, AssertionFailureException);
  p.soft_add_kmer_event(kmer, e3);
  EXPECT_EQ(kmer, p.get_pos_kmer(kmer)->kmer);
  EXPECT_EQ(1, p.position);
  EXPECT_FLOAT_EQ(2.2, p.get_pos_kmer(kmer)->events.front().posterior_probability);
}

TEST (BaseKmerTests, test_ContigStrand) {
  Redirect a(true, true);
  string kmer = "ATGCC";
  Event e(1, 2);
  Event e2(1.1, 2.2);
  shared_ptr<PosKmer> k = make_shared<PosKmer>(kmer, 1);
  k->add_event(e);
  k->add_event(e2);
  string contig = "asd";
  string strand = "+";
  uint64_t num_positions = 10;
  ContigStrand cs(contig, strand, num_positions, "t");
  EXPECT_EQ(num_positions, cs.positions.capacity());
  EXPECT_EQ("asd", cs.get_contig());
  EXPECT_EQ("+", cs.get_strand());
  for (uint64_t i=0; i < num_positions; i++){
    EXPECT_EQ(i, cs.positions[i].position);
  }
  uint64_t pos = 1;
  cs.add_kmer(pos, k);
  EXPECT_FLOAT_EQ(2.2, cs.positions[pos].get_pos_kmer(kmer)->events.front().posterior_probability);
  EXPECT_FLOAT_EQ(2.2, cs.get_position(pos).get_pos_kmer(kmer)->events.front().posterior_probability);
  float c = 2;
  float b = 3.3;
  cs.add_event(pos, kmer, c, b);
  EXPECT_FLOAT_EQ(3.3, cs.get_position(pos).get_pos_kmer(kmer)->events.front().posterior_probability);
}

TEST (BaseKmerTests, test_Kmer) {
  Redirect a(true, true);
  Event e(1, 2);
  Event e2(1.1, 2.2);
  shared_ptr<PosKmer> k = make_shared<PosKmer>("ATGCC", 1);

  k->add_event(e);
  k->add_event(e2);
  k->add_event(1.02, 2.1);
  EXPECT_EQ("ATGCC", k->kmer);
  EXPECT_FLOAT_EQ(2.2, k->events.front().posterior_probability);
  EXPECT_EQ(1, k->num_events());
  Kmer kmer("ATGCC");
  string contig_strand = "asdf+t";
  uint64_t pos = 1;
  kmer.add_pos_kmer(contig_strand, pos, k);
  EXPECT_EQ(contig_strand, get<0>(kmer.contig_positions[0]));
  EXPECT_EQ(pos, get<1>(kmer.contig_positions[0]));
  EXPECT_EQ(k, kmer.get_pos_kmer(contig_strand, pos));
  EXPECT_EQ(true, kmer.has_pos_kmer(contig_strand, pos));
  EXPECT_EQ(false, kmer.has_pos_kmer(contig_strand, 2));
  kmer.soft_add_pos_kmer(contig_strand, pos, k);
  EXPECT_EQ(1, kmer.contig_positions.size());
  EXPECT_EQ(1, kmer.pos_kmer_map.size());
  ContigStrandPosition csp = kmer.split_pos_kmer_map_key("asdf+t1234");
  EXPECT_EQ("asdf", csp.contig);
  EXPECT_EQ("+", csp.strand);
  EXPECT_EQ("t", csp.nanopore_strand);
  EXPECT_EQ(1234, csp.position);

}

TEST (BaseKmerTests, test_ByKmer_kde) {
  Redirect a(true, true);
  shared_ptr<PosKmer> k = make_shared<PosKmer>("ATGCC");
  k->add_event(1, 1);
  k->add_event(2, 0.5);
  k->add_event(3, 1);
  k->add_event(4, 0.5);
  k->add_event(5, 1);
  k->add_event(6, 0.5);
  k->add_event(7, 1);
  k->add_event(8, 0.5);
  k->add_event(9, 1);
  k->add_event(10, 1);
  shared_ptr<PosKmer> k2 = make_shared<PosKmer>("ATGCC");
  k2->add_event(1, 1);
  k2->add_event(2, 0.5);
  k2->add_event(3, 1);
  k2->add_event(4, 0.5);
  k2->add_event(5, 1);
  k2->add_event(6, 0.5);
  k2->add_event(7, 1);
  k2->add_event(8, 0.5);
  k2->add_event(9, 1);
  k2->add_event(10, 1);
  Kmer kmer("ATGCC");
  kmer.add_pos_kmer("asdf", 1, k);
  kmer.add_pos_kmer("asdf", 2, k2);
  ASSERT_THROW(kmer.get_hist(0,10,10), AssertionFailureException);
  vector<uint64_t> hist{2,2,2,2,2,2,2,2,2,2};
  ASSERT_THAT(hist, ElementsAreArray(kmer.get_hist(0,10.1,10)));
  ASSERT_THAT(hist, ElementsAreArray(kmer.get_hist(0,10.1,10, 0.5)));
  hist = {2,0,2,0,2,0,2,0,2,2};
  ASSERT_THAT(hist, ElementsAreArray(kmer.get_hist(0,10.1,10, 0.51)));
  hist = {0,2,2,2,2,2,2,2,2,2,2};
  ASSERT_THAT(hist, ElementsAreArray(kmer.get_hist(0,11,11)));
}

TEST (BaseKmerTests, test_ByKmer) {
  Redirect a(true, true);
  Event e(1, 2);
  Event e2(1.1, 2.2);
  shared_ptr<PosKmer> k = make_shared<PosKmer>("ATGCC", 1);
  k->add_event(e);
  k->add_event(e2);
  k->add_event(1.02, 2.1);
  EXPECT_EQ("ATGCC", k->kmer);
  EXPECT_FLOAT_EQ(2.2, k->events.front().posterior_probability);
  EXPECT_EQ(1, k->num_events());
  string contig_strand = "asdf";
  uint64_t pos = 1;
  uint64_t kmer_length = 5;
  ByKmer by_kmer({'A', 'C', 'T', 'G'}, kmer_length);
  for (auto &k: all_string_permutations("ACTG", kmer_length)) {
    by_kmer.get_kmer(k);
  }
  ASSERT_THROW(by_kmer.get_kmer("AAFAA"), AssertionFailureException);
  by_kmer.add_kmer_ptr(contig_strand, pos, k);
  EXPECT_EQ(1, by_kmer.get_kmer("ATGCC").pos_kmer_map.size());
  by_kmer.add_kmer_ptr(contig_strand, pos, k);
  EXPECT_EQ(1, by_kmer.get_kmer("ATGCC").pos_kmer_map.size());
}


#endif //EMBED_FAST5_TESTS_SRC_BASEKMERTESTS_HPP_
