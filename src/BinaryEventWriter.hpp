//
// Created by Andrew Bailey on 5/2/20.
// Most of this code came from Ryan Lorig-Roach. https://github.com/rlorigro/runlength_analysis_cpp
//

#ifndef EMBED_FAST5_SRC_BINARYEVENTWRITER_HPP_
#define EMBED_FAST5_SRC_BINARYEVENTWRITER_HPP_

// embed libs
#include "BaseKmer.hpp"
#include "BinaryIO.hpp"
// boost libs
#include <boost/filesystem.hpp>
// std libs
#include <string>
#include <vector>
#include <ostream>

using namespace std;
using namespace embed_utils;
using namespace boost::filesystem;

class KmerIndex {
 public:
  /// Attributes ///
  string name;
  uint64_t name_length;
  uint64_t sequence_byte_index;
  uint64_t sequence_length;
};

class PositionIndex {
 public:
  /// Attributes ///
  uint64_t position;
  unordered_map<string, KmerIndex> kmer_indexes;
  uint64_t num_kmers;

  KmerIndex& get_kmer_index(const string& kmer){
    auto found = kmer_indexes.find(kmer);
    if (found != kmer_indexes.end()) {
      // Found it
      return found->second;
    } else {
      throw runtime_error(kmer + " was not found in kmer index map");
    }
  }

  vector<string> get_kmers(){
    vector<string> v;
    for (auto i : kmer_indexes) {
      v.push_back(i.first);
    }
    return v;
  }

};

class ContigStrandIndex {
 public:
  /// Attributes ///
  string contig;
  uint64_t contig_string_length;
  string strand;
  string nanopore_strand;
  uint64_t num_positions;
  unordered_map<uint64_t, PositionIndex> position_indexes;
  uint64_t num_written_positions;

  PositionIndex& get_position_index(const uint64_t& position){
    auto found = position_indexes.find(position);
    if (found != position_indexes.end()) {
      // Found it
      return found->second;
    } else {
      throw runtime_error(to_string(position) + " was not found in position index map");
    }
  }
  KmerIndex& get_kmer_index(const uint64_t& position, const string& kmer) {
    return get_position_index(position).get_kmer_index(kmer);
  }
};


//ostream& operator<<(ostream& s, KmerIndex& index);


class BinaryEventWriter {
 public:
  /// Attributes ///
  path sequence_file_path;
  std::ofstream sequence_file;
  // When writing the binary file, this vector is appended, so the position of each sequence is stored
  vector<ContigStrandIndex> contig_strand_indexes;

  PositionIndex write_position(Position& position){
    PositionIndex index;
    // Store the number of kmers
    index.num_kmers = position.kmers.size();
    // Store the name of this sequence
    index.position = position.position;

    for (auto &kmer: position.kmers){
      index.kmer_indexes.insert(std::make_pair(kmer.first, this->write_kmer(kmer.second)));
    }
    return index;
  }

  KmerIndex write_kmer(Kmer& kmer){
    if (kmer.events.empty()){
      throw runtime_error("ERROR: empty sequence provided to BinaryEventWriter: " + kmer.kmer);
    }
    KmerIndex index;
    // Add sequence start position to index
    index.sequence_byte_index = this->sequence_file.tellp();
    // Store the length of this sequence
    index.sequence_length = kmer.events.size();
    // Store the name of this sequence
    index.name = kmer.kmer;
    write_vector_to_binary(this->sequence_file, kmer.events);
    // Append index object to vector
    return index;
  }


  void write_kmer_index(KmerIndex& index){
    // Where is the sequence
    write_value_to_binary(this->sequence_file, index.sequence_byte_index);
    // How long is the sequence
    write_value_to_binary(this->sequence_file, index.sequence_length);
    // How long is the name of the sequence
    write_value_to_binary(this->sequence_file, index.name.size());
    // What is the name
    write_string_to_binary(this->sequence_file, index.name);
  }

  void write_position_index(PositionIndex& index){
    // write position
    write_value_to_binary(this->sequence_file, index.position);
    // How many kmers
    write_value_to_binary(this->sequence_file, index.num_kmers);
    for (auto &kmer_index: index.kmer_indexes){
      write_kmer_index(kmer_index.second);
    }
  }

  void write_contig_strand_index(ContigStrandIndex& index){
    // write contig name length
    write_value_to_binary(this->sequence_file, index.contig_string_length);
    // Contig, strand and nanopore strand
    write_string_to_binary(this->sequence_file, index.contig);
    write_string_to_binary(this->sequence_file, index.strand);
    write_string_to_binary(this->sequence_file, index.nanopore_strand);
    write_value_to_binary(this->sequence_file, index.num_positions);
    write_value_to_binary(this->sequence_file, index.num_written_positions);
    for (auto &pos_index: index.position_indexes){
      write_position_index(pos_index.second);
    }
  }

 public:
  /// Methods ///
  BinaryEventWriter(path file_path) {
    this->sequence_file_path = file_path;
    // Ensure that the output directory exists
    create_directories(this->sequence_file_path.parent_path());
    this->sequence_file = std::ofstream(this->sequence_file_path.c_str(), std::ofstream::binary);
    throw_assert(this->sequence_file.is_open(), "ERROR: could not open file " + file_path.string());
  }

  ~BinaryEventWriter() {
    this->close();
  }

  void close(){
    sequence_file.close();
  }

  void write_contig_strand(ContigStrand& contig){
    ContigStrandIndex index;
    // Store contig strand attributes
    index.contig = contig.contig;
    index.contig_string_length = contig.contig.length();
    index.strand = contig.strand;
    index.nanopore_strand = contig.nanopore_strand;
    index.num_positions = contig.num_positions;
    index.num_written_positions = 0;

    for (auto &position: contig.positions){
      if (position.has_data){
        index.position_indexes.insert(std::make_pair(position.position, this->write_position(position)));
        index.num_written_positions += 1;
      }
    }
    contig_strand_indexes.push_back(index);
  }

  void write_indexes(){
    // Store the current file byte index so the beginning of the INDEX table can be located later
    uint64_t indexes_start_position = this->sequence_file.tellp();

    // Iterate all the indexes, write them to the file
    for (auto& index: this->contig_strand_indexes){
      write_contig_strand_index(index);
    }
    // Store the current file byte index so the beginning of the CHANNEL table can be located later
//    uint64_t channel_metadata_start_position = this->sequence_file.tellp();
    // Write the pointer to the beginning of the index table
    write_value_to_binary(this->sequence_file, indexes_start_position);
  }
};

#endif //EMBED_FAST5_SRC_BINARYEVENTWRITER_HPP_
