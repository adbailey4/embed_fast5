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

TEST (BaseKmerTests, test_Kmer) {
  Redirect a(true, true);
  Event e(1, 2);
  Event e2(1.1, 2.2);
  Kmer k("ATGCC", 1);
  k.add_event(e);
  k.add_event(e2);
  k.add_event(1.02, 2.1);
  EXPECT_EQ("ATGCC", k.kmer);
  EXPECT_FLOAT_EQ(2.2, k.events.front().posterior_probability);
  EXPECT_EQ(1, k.num_events());
  Kmer k2("ATGCC", 2);
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
  Kmer k(kmer, 1);
  k.add_event(e);
  k.add_event(e2);
  Kmer k2;
  k2 = k;
  Position p(1);
//  p.kmers.emplace(kmer, move(k));
//  cout << "emplaced"  << '\n';
  p.add_kmer(k);
  cout << "add kmer"  << '\n';
//  p.kmers.insert({"ATGCC", k});
//  cout << "inserted"  << '\n';
//  p.add_kmer(move(k));
  ASSERT_THROW({(p.add_kmer)(k2);}, AssertionFailureException);
  p.soft_add_kmer_event(kmer, e3);
  EXPECT_EQ(kmer, p.get_kmer(kmer).kmer);
  EXPECT_EQ(1, p.position);
  EXPECT_FLOAT_EQ(2.2, p.get_kmer(kmer).events.front().posterior_probability);
}

TEST (BaseKmerTests, test_ContigStrand) {
  Redirect a(true, true);
  string kmer = "ATGCC";
  Event e(1, 2);
  Event e2(1.1, 2.2);
  Kmer k(kmer, 1);
  k.add_event(e);
  k.add_event(e2);
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
  EXPECT_FLOAT_EQ(2.2, cs.positions[pos].kmers[kmer].events.front().posterior_probability);
  EXPECT_FLOAT_EQ(2.2, cs.get_position(pos).kmers[kmer].events.front().posterior_probability);
  float c = 2;
  float b = 3.3;
  cs.add_event(pos, kmer, c, b);
  EXPECT_FLOAT_EQ(3.3, cs.get_position(pos).kmers[kmer].events.front().posterior_probability);
}



#endif //EMBED_FAST5_TESTS_SRC_BASEKMERTESTS_HPP_
