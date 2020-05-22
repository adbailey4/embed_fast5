//
// Created by Andrew Bailey on 5/3/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_

// embed source
#include "PerPositionKmers.hpp"
#include "TestFiles.hpp"
#include "BinaryEventReader.hpp"
#include "SplitByRefPosition.hpp"

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
  Kmer k2("ATGCC", 2);
  k2.add_event(e);
  k2.add_event(e2);
  k2.add_event(1.02, 2.1);
  EXPECT_EQ("ATGCC", k2.kmer);
  EXPECT_FLOAT_EQ(2.2, k2.events.back().posterior_probability);
  EXPECT_EQ(2, k2.num_events());
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

TEST (PerPositionKmersTests, test_initialization) {
  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "test.event";
  uint64_t num_locks = 1000;
  ReferenceHandler reference(PUC_REFERENCE.string());
  PerPositionKmers ppk(num_locks, reference, true);
  string contig_strand = "pUC19+c";
  uint64_t length = reference.get_chromosome_sequence_length("pUC19");
  EXPECT_EQ("pUC19", ppk.data.data[contig_strand].contig);
  EXPECT_EQ(length, ppk.data.data[contig_strand].num_positions);
}

TEST (PerPositionKmersTests, test_process_alignment) {
  Redirect a(true, true);
  uint64_t num_locks = 1000;
  ReferenceHandler reference(PUC_REFERENCE.string());
  PerPositionKmers ppk(num_locks, reference, true);
  path alignment_file = PUC_5MER_ALIGNMENTS/"03274a9a-0eab-422e-ace7-b35fd3a0f48c.sm.forward.tsv";
  AlignmentFile af(alignment_file.string());
  ppk.process_alignment(af);
  string contig_strand = "pUC19+c";
  uint64_t position = 1770;
  string kmer = "ATTGA";
  EXPECT_EQ(3, ppk.data.get_kmer("pUC19", "+", "c", position, kmer).num_events());
  position = 2681;
  EXPECT_EQ(1, ppk.data.get_kmer("pUC19", "+", "c", position, kmer).num_events());
}

TEST (BinaryEventTests, test_read_and_write) {
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
  Kmer k(kmer, 2);
  k.add_event(1, 1);
  k.add_event(2, .5);
  ContigStrand cs(contig, strand, num_positions, "t");
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
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes.size());
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes[kmer].name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[pos].kmer_indexes[kmer].name_length);

  Kmer kmer_struct;
  ber.get_kmer(kmer_struct, kmer, contig, strand, pos);
  EXPECT_EQ(kmer, kmer_struct.kmer);
  ASSERT_THAT(k.events, ElementsAreArray(kmer_struct.events));
}

TEST (PerPositionKmersTests, test_write_to_file) {
  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
  uint64_t num_locks = 1000;
  ReferenceHandler reference(PUC_REFERENCE.string());
  PerPositionKmers ppk(num_locks, reference, true);
  path alignment_file = PUC_5MER_ALIGNMENTS/"03274a9a-0eab-422e-ace7-b35fd3a0f48c.sm.forward.tsv";
  AlignmentFile af(alignment_file.string());
  ppk.process_alignment(af);
  ppk.write_to_file(test_file);

  BinaryEventReader ber(test_file.string());

  string contig_strand = "pUC19+c";
  uint64_t position = 1770;
  string kmer = "ATTGA";
  EXPECT_EQ(4, ber.indexes.size());
  EXPECT_EQ(5, ber.indexes[contig_strand].contig_string_length);
  EXPECT_EQ("+", ber.indexes[contig_strand].strand);
  EXPECT_EQ("c", ber.indexes[contig_strand].nanopore_strand);
  EXPECT_EQ("pUC19", ber.indexes[contig_strand].contig);
  EXPECT_EQ(2686, ber.indexes[contig_strand].num_positions);
  EXPECT_EQ(2514, ber.indexes[contig_strand].num_written_positions);
  EXPECT_EQ(2514, ber.indexes[contig_strand].position_indexes.size());
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[position].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[position].kmer_indexes[kmer].name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[position].kmer_indexes[kmer].name_length);
  vector<string> kmers = ber.get_position_index("pUC19", "+", "c", position).get_kmers();
  ASSERT_THAT(kmers, ElementsAreArray(ber.indexes[contig_strand].position_indexes[position].get_kmers()));
  Kmer kmer_struct;
  ber.get_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(kmer, kmer_struct.kmer);
  EXPECT_EQ(3, kmer_struct.num_events());
  position = 2681;
  ber.get_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(1, kmer_struct.num_events());
}

TEST (PerPositionKmersTests, test_split_signal_align_by_ref_position) {
//  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
  string sa_input_dir = PUC_5MER_ALIGNMENTS.string();
  string output_file_path = test_file.string();
  string reference = PUC_REFERENCE.string();
  uint64_t num_locks = 1000;
  uint64_t n_threads = 4;
  bool verbose = false;
  bool rna = false;
  bool two_d = true;
  split_signal_align_by_ref_position(sa_input_dir, output_file_path, reference,
      num_locks, n_threads, verbose, rna, two_d);

  BinaryEventReader ber(test_file.string());

  string contig_strand = "pUC19+c";
  uint64_t position = 1770;
  string kmer = "ATTGA";
  EXPECT_EQ(4, ber.indexes.size());
  EXPECT_EQ(5, ber.indexes[contig_strand].contig_string_length);
  EXPECT_EQ("+", ber.indexes[contig_strand].strand);
  EXPECT_EQ("c", ber.indexes[contig_strand].nanopore_strand);
  EXPECT_EQ("pUC19", ber.indexes[contig_strand].contig);
  EXPECT_EQ(2686, ber.indexes[contig_strand].num_positions);
  EXPECT_EQ(2682, ber.indexes[contig_strand].num_written_positions);
  EXPECT_EQ(2682, ber.indexes[contig_strand].position_indexes.size());
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[position].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[position].kmer_indexes[kmer].name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[position].kmer_indexes[kmer].name_length);
  vector<string> kmers = ber.get_position_index("pUC19", "+", "c", position).get_kmers();
  ASSERT_THAT(kmers, ElementsAreArray(ber.indexes[contig_strand].position_indexes[position].get_kmers()));
  Kmer kmer_struct;
  ber.get_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(kmer, kmer_struct.kmer);
  EXPECT_EQ(22, kmer_struct.num_events());
  position = 2681;
  ber.get_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(5, kmer_struct.num_events());
}


#endif //EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
