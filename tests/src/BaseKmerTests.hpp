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
  EXPECT_EQ(kmer, p.get_kmer(kmer)->kmer);
  EXPECT_EQ(1, p.position);
  EXPECT_FLOAT_EQ(2.2, p.get_kmer(kmer)->events.front().posterior_probability);
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
  EXPECT_FLOAT_EQ(2.2, cs.positions[pos].get_kmer(kmer)->events.front().posterior_probability);
  EXPECT_FLOAT_EQ(2.2, cs.get_position(pos).get_kmer(kmer)->events.front().posterior_probability);
  float c = 2;
  float b = 3.3;
  cs.add_event(pos, kmer, c, b);
  EXPECT_FLOAT_EQ(3.3, cs.get_position(pos).get_kmer(kmer)->events.front().posterior_probability);
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
  string contig_strand = "asdf";
  uint64_t pos = 1;
  kmer.add_pos_kmer(contig_strand, pos, k);
  EXPECT_EQ(contig_strand, kmer.contig_strands[0]);
  EXPECT_EQ(pos, kmer.positions[0]);
  EXPECT_EQ(k, kmer.kmers[0]);
  EXPECT_EQ(0, kmer.find_index(contig_strand, pos));
  EXPECT_EQ(-1, kmer.find_index(contig_strand, 2));
  kmer.soft_add_pos_kmer(contig_strand, pos, k);
  EXPECT_EQ(1, kmer.contig_strands.size());
  EXPECT_EQ(1, kmer.positions.size());
  EXPECT_EQ(1, kmer.kmers.size());

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
  EXPECT_EQ(1, by_kmer.get_kmer("ATGCC").kmers.size());
  by_kmer.add_kmer_ptr(contig_strand, pos, k);
  EXPECT_EQ(1, by_kmer.get_kmer("ATGCC").kmers.size());
}


#endif //EMBED_FAST5_TESTS_SRC_BASEKMERTESTS_HPP_
