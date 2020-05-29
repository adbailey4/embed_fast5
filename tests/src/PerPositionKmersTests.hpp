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
  EXPECT_EQ(3, ppk.data.get_position_kmer("pUC19", "+", "c", position, kmer).num_events());
  position = 2681;
  EXPECT_EQ(1, ppk.data.get_position_kmer("pUC19", "+", "c", position, kmer).num_events());
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
  ber.get_position_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(kmer, kmer_struct.kmer);
  EXPECT_EQ(3, kmer_struct.num_events());
  position = 2681;
  ber.get_position_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
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

  BinaryEventReader ber(output_file_path);

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
  ber.get_position_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(kmer, kmer_struct.kmer);
  EXPECT_EQ(22, kmer_struct.num_events());
  position = 2681;
  ber.get_position_kmer(kmer_struct, kmer, "pUC19", "+", position, "c");
  EXPECT_EQ(5, kmer_struct.num_events());
}

TEST (PerPositionKmersTests, test_rna_reads) {
//  Redirect a(true, true);
  path tempdir = temp_directory_path() / "temp";
  path test_file = tempdir / "rna_test.event";
  if (exists(test_file)){
    remove(test_file);
  }
  string sa_input_dir = RRNA_SIGNAL_FILES.string();
  string output_file_path = test_file.string();
  string reference = ECOLI_16S_REFERENCE.string();

  uint64_t num_locks = 1000;
  uint64_t n_threads = 4;
  bool verbose = false;
  bool rna = true;
  bool two_d = true;
  split_signal_align_by_ref_position(sa_input_dir, output_file_path, reference,
                                     num_locks, n_threads, verbose, rna, two_d);

  BinaryEventReader ber(output_file_path);
  string contig_strand = "ecoli_MRE600+t";
  uint64_t position = 200;
  string kmer = "AGGGG";
  EXPECT_EQ(4, ber.indexes.size());
  EXPECT_EQ(12, ber.indexes[contig_strand].contig_string_length);
  EXPECT_EQ("+", ber.indexes[contig_strand].strand);
  EXPECT_EQ("t", ber.indexes[contig_strand].nanopore_strand);
  EXPECT_EQ("ecoli_MRE600", ber.indexes[contig_strand].contig);
  EXPECT_EQ(1542, ber.indexes[contig_strand].num_positions);
  EXPECT_EQ(1527, ber.indexes[contig_strand].num_written_positions);
  EXPECT_EQ(1527, ber.indexes[contig_strand].position_indexes.size());
  EXPECT_EQ(1, ber.indexes[contig_strand].position_indexes[position].kmer_indexes.size());
  EXPECT_EQ(kmer, ber.indexes[contig_strand].position_indexes[position].kmer_indexes[kmer].name);
  EXPECT_EQ(5, ber.indexes[contig_strand].position_indexes[position].kmer_indexes[kmer].name_length);
  vector<string> kmers = ber.get_position_index("ecoli_MRE600", "+", "t", position).get_kmers();
  ASSERT_THAT(kmers, ElementsAreArray(ber.indexes[contig_strand].position_indexes[position].get_kmers()));
  Kmer kmer_struct;
  ber.get_position_kmer(kmer_struct, kmer, "ecoli_MRE600", "+", position, "t");
  EXPECT_EQ(kmer, kmer_struct.kmer);
  EXPECT_EQ(16, kmer_struct.num_events());
  position = 144;
  ber.get_position_kmer(kmer_struct, kmer, "ecoli_MRE600", "+", position, "t");
  EXPECT_EQ(9, kmer_struct.num_events());
}


#endif //EMBED_FAST5_TESTS_SRC_PERPOSITIONKMERSTESTS_HPP_
