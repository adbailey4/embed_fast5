//
// Created by Andrew Bailey on 5/3/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_

// embed source
#include "PerPositionKmers.hpp"
#include "TestFiles.hpp"
#include "BinaryEventWriter.hpp"
#include "BinaryEventReader.hpp"
//boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace boost::filesystem;
using namespace std;
using namespace embed_utils;
using namespace test_files;

TEST (PerPositionKmersTests, test_Event) {
  Redirect a(true, true);
  Event e(1, 2);
  EXPECT_EQ(1, e.descaled_event_mean);
  EXPECT_EQ(2, e.posterior_probability);
}

TEST (PerPositionKmersTests, test_Kmer) {
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
}

TEST (PerPositionKmersTests, test_Position) {
  Redirect a(true, true);
  string kmer = "ATGCC";
  Event e(1, 2);
  Event e2(1.1, 2.2);
  Kmer k(kmer, 1);
  k.add_event(e);
  k.add_event(e2);
  Position p(1);
  p.add_kmer(k);
  ASSERT_THROW({(p.add_kmer)(k);}, AssertionFailureException);
  p.soft_add_kmer_event(kmer, e);
  EXPECT_EQ(kmer, p.get_kmer(kmer).kmer);
  EXPECT_EQ(1, p.position);
  EXPECT_FLOAT_EQ(2.2, k.events.front().posterior_probability);
}


TEST (PerPositionKmersTests, test_ContigStrand) {
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
  uint64_t pos_offset = 2;
  ContigStrand cs(contig, strand, num_positions, pos_offset);
  EXPECT_EQ(num_positions, cs.positions.capacity());
  EXPECT_EQ("asd", cs.get_contig());
  EXPECT_EQ("+", cs.get_strand());
  uint64_t pos = 1;
  cs.add_kmer(pos, k);
  EXPECT_FLOAT_EQ(2.2, cs.positions[pos+pos_offset].kmers[kmer].events.front().posterior_probability);
  EXPECT_FLOAT_EQ(2.2, cs.get_position(pos).kmers[kmer].events.front().posterior_probability);
  float c = 2;
  float b = 3.3;
  cs.add_event(pos, kmer, c, b);
  EXPECT_FLOAT_EQ(3.3, cs.get_position(pos).kmers[kmer].events.front().posterior_probability);
}

TEST (PerPositionKmersTests, test_initialization) {
  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "test.event";
  uint64_t num_locks = 1000;
  ReferenceHandler reference(PUC_REFERENCE.string());
  PerPositionKmers ppk(num_locks, test_file, reference, true);
  string contig_strand = "pUC19+c";
  uint64_t length = reference.get_chromosome_sequence_length("pUC19");
  EXPECT_EQ("pUC19", ppk.data[contig_strand].contig);
  EXPECT_EQ(length, ppk.data[contig_strand].num_positions);
}

TEST (PerPositionKmersTests, test_process_alignment) {
  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "test.event";
  uint64_t num_locks = 1000;
  ReferenceHandler reference(PUC_REFERENCE.string());
  PerPositionKmers ppk(num_locks, test_file, reference);
  path alignment_file = PUC_5MER_ALIGNMENTS/"03274a9a-0eab-422e-ace7-b35fd3a0f48c.sm.forward.tsv";
  AlignmentFile af(alignment_file.string());
  ppk.process_alignment(af);
  string contig_strand = "pUC19+c";
  uint64_t position = 1770;
  string kmer = "ATTGA";
  EXPECT_EQ(3, ppk.data[contig_strand].get_position(position).get_kmer(kmer).num_events());
  position = 2681;
  EXPECT_EQ(1, ppk.data[contig_strand].get_position(position).get_kmer(kmer).num_events());
}

TEST (BinaryEventTests, test_process_alignment) {
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
  uint64_t pos_offset = 2;
  Kmer k(kmer, 2);
  k.add_event(1, 1);
  k.add_event(2, .5);
  ContigStrand cs(contig, strand, num_positions, pos_offset);
  uint64_t pos = 1;
  cs.add_kmer(pos, k);
//  write data structure
  BinaryEventWriter bew(test_file);
  bew.write_contig_strand(cs);
  bew.write_indexes();
  bew.close();
  string contig_strand = "asd+t";
  BinaryEventReader ber(test_file.string());
  EXPECT_EQ(1, ber.indexes.size());
  EXPECT_EQ(3, ber.indexes[contig_strand].contig_string_length);
  EXPECT_EQ(strand, ber.indexes[contig_strand].strand);
  EXPECT_EQ("t", ber.indexes[contig_strand].nanopore_strand);
  EXPECT_EQ(contig, ber.indexes[contig_strand].contig);
  EXPECT_EQ(num_positions, ber.indexes[contig_strand].num_positions);
  EXPECT_EQ(1, ber.indexes[contig_strand].num_written_positions);
  EXPECT_EQ(pos_offset, ber.indexes[contig_strand].pos_offset);
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes.size());
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes[kmer].name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes[kmer].name_length);
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[0].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[0].kmer_indexes[kmer].name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[0].kmer_indexes[kmer].name_length);

//  Kmer kmer_struct;
//  ber.get_kmer(kmer_struct, kmer, contig, strand, pos);
//  EXPECT_EQ(kmer, kmer_struct.kmer);
//  ASSERT_THAT(k.events, ElementsAreArray(kmer_struct.events));
}


#endif //EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
