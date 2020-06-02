//
// Created by Andrew Bailey on 5/23/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_BINARYEVENTTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_BINARYEVENTTESTS_HPP_

// embed source
#include "BinaryEventReader.hpp"
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


TEST (BinaryEventTests, test_read_and_write) {
//  Redirect a(true, true);
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  create data structure
  string kmer = "ATGCC";
  string contig = "asd";
  string strand = "+";
  uint64_t num_positions = 10;
//  PosKmer k(kmer, 2);
  shared_ptr<PosKmer> k = make_shared<PosKmer>(kmer, 2);
  k->add_event(1, 1);
  k->add_event(2, .5);
  shared_ptr<PosKmer> k2 = k;
  ContigStrand cs(contig, strand, num_positions, "t");
  uint64_t pos = 1;
  cs.add_kmer(pos, k);
//  write data structure
  std::set<char> alphabet{'A', 'G', 'C', 'T'};
  BinaryEventWriter bew(test_file, alphabet, 5);
  bew.write_contig_strand(cs);
  bew.write_indexes();
  bew.close();
  string contig_strand = "asd+t";
  BinaryEventReader ber(test_file.string());
  EXPECT_EQ("ACGT", ber.alphabet_string);
  ASSERT_THAT(alphabet, ElementsAreArray(ber.alphabet));
  EXPECT_EQ(5, ber.kmer_length);

  EXPECT_EQ(1, ber.indexes.size());
  EXPECT_EQ(3, ber.indexes[contig_strand].contig_string_length);
  EXPECT_EQ(strand, ber.indexes[contig_strand].strand);
  EXPECT_EQ("t", ber.indexes[contig_strand].nanopore_strand);
  EXPECT_EQ(contig, ber.indexes[contig_strand].contig);
  EXPECT_EQ(num_positions, ber.indexes[contig_strand].num_positions);
  EXPECT_EQ(1, ber.indexes[contig_strand].num_written_positions);
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes.size());
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes[kmer]->name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes[kmer]->name_length);

  shared_ptr<PosKmer> kmer_struct = ber.get_position_kmer(kmer, "asd", "+", pos, "t");
  EXPECT_EQ(kmer, kmer_struct->kmer);
  ASSERT_THAT(k2->events, ElementsAreArray(kmer_struct->events));
}

TEST (BinaryEventTests, test_create_kmer_map) {
  Redirect a(true, true);
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  create data structure
  string kmer = "ATGCC";
  string contig = "asd";
  string strand = "+";
  uint64_t num_positions = 10;
  shared_ptr<PosKmer> k = make_shared<PosKmer>(kmer, 2);
  k->add_event(1, 1);
  k->add_event(2, .5);
  shared_ptr<PosKmer> k2 = k;
  ContigStrand cs(contig, strand, num_positions, "t");
  uint64_t pos = 1;
  cs.add_kmer(pos, k);
  cs.add_kmer(2, k2);

//  write data structure
  std::set<char> alphabet{'A', 'G', 'C', 'T'};
  BinaryEventWriter bew(test_file, alphabet, 5);
  bew.write_contig_strand(cs);
  bew.write_indexes();
  bew.close();
  string contig_strand = "asd+t";
  BinaryEventReader ber(test_file.string());
  EXPECT_EQ(2, ber.kmer_map.get_kmer_index(kmer).kmer_index_ptrs.size());
  Kmer kmer_struct("ATGCC");
  ber.populate_kmer(kmer_struct);
  EXPECT_EQ(2, kmer_struct.kmers.size());
}


#endif //EMBED_FAST5_TESTS_SRC_BINARYEVENTTESTS_HPP_
