//
// Created by Andrew Bailey on 5/3/20.
// Most of this code came from Ryan Lorig-Roach. https://github.com/rlorigro/runlength_analysis_cpp
//

#ifndef EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_
#define EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_

// embed libs
#include "BinaryEventWriter.hpp"
#include <utility>

using namespace std;
using namespace embed_utils;
using namespace boost::filesystem;

using namespace std;

class BinaryEventReader {
 public:

  /// Methods ///

  // Initialize the class with a file path
  BinaryEventReader(string file_path) {
    this->initialize(file_path);
  }
  BinaryEventReader() {}
  ~BinaryEventReader() {
    this->close();
  }

  void initialize(const string& file_path){
    this->initialized = true;
    this->sequence_file_path = file_path;

    // Open the input file.
    this->sequence_file_descriptor = ::open(file_path.c_str(), O_RDONLY);

    // Verify it is working
    if(this->sequence_file_descriptor == -1) {
      throw runtime_error("ERROR: could not read " + file_path);
    }

    // Find file size in bytes
    this->file_length = lseek(this->sequence_file_descriptor, 0, SEEK_END);

    // Initialize remaining parameters using the file footer data
    this->read_footer();

    // Read table of contents, needed for indexed reading
    this->read_indexes();
  }

  void close() {
    ::close(sequence_file_descriptor);
  }

  void get_kmer(Kmer& kmer_struct, const string& kmer, const string& contig_name, const string& strand,
      const uint64_t& position, const string nanopore_strand="t"){
    KmerIndex& ki = this->get_kmer_index(contig_name, strand, nanopore_strand, position, kmer);
    kmer_struct.kmer = ki.name;
    off_t byte_index = ki.sequence_byte_index;
    uint64_t sequence_length = ki.sequence_length;
    kmer_struct.events.reserve(sequence_length);
    pread_vector_from_binary(this->sequence_file_descriptor, kmer_struct.events, sequence_length, byte_index);
  }

  void get_position(Position& position_struct, const string& contig_name, const string& strand,
                const uint64_t& position, const string nanopore_strand="t"){
    PositionIndex& ki = this->get_position_index(contig_name, strand, nanopore_strand, position);
    vector<string> kmers= ki.get_kmers();
    for (auto kmer: kmers){
      Kmer k;
      get_kmer(k, kmer, contig_name, strand, position, nanopore_strand);
      position_struct.add_kmer(k);
    }
  }

  ContigStrandIndex& get_contig_index(const string& contig, const string& strand, const string& nanopore_strand){
    string contig_strand = contig+strand+nanopore_strand;
    auto found = indexes.find(contig_strand);
    if (found != indexes.end()) {
      // Found it
      return found->second;
    } else {
      throw runtime_error(contig_strand + " was not found in indexes");
    }
  }

  PositionIndex& get_position_index(const string& contig, const string& strand, const string& nanopore_strand, const uint64_t& position){
    return this->get_contig_index(contig, strand, nanopore_strand).get_position_index(position);
  }

  KmerIndex& get_kmer_index(const string& contig, const string& strand, const string& nanopore_strand, const uint64_t& position, const string& kmer){
    return this->get_position_index(contig, strand, nanopore_strand, position).get_kmer_index(kmer);
  }



  /// Attributes ///
  unordered_map<string, ContigStrandIndex> indexes;
  bool initialized = false;

 private:

  /// Attributes ///
  string sequence_file_path;
  int sequence_file_descriptor;

  uint64_t indexes_start_position;
  off_t file_length;

  /// Methods ///
  void read_indexes(){
    off_t byte_index = off_t(this->indexes_start_position);
    while (byte_index > 0 and uint64_t(byte_index) < (this->file_length - 1*sizeof(uint64_t))){
      ContigStrandIndex index_element;
      this->read_contig_strand_index_entry(index_element, byte_index);
      auto ret = this->indexes.emplace(index_element.contig + index_element.strand + index_element.nanopore_strand, move(index_element));
      throw_assert(ret.second,
          "ERROR: possible duplicate read name (" + ret.first->first + ") found in events file: " +
          this->sequence_file_path)
    }
  }

  void read_footer(){
    off_t byte_index = off_t(this->file_length - 1*sizeof(uint64_t));
    pread_value_from_binary(this->sequence_file_descriptor, this->indexes_start_position, byte_index);
  }
  
  void read_kmer_index_entry(KmerIndex& index_element, off_t& byte_index){
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_byte_index, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_length, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.name_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.name, index_element.name_length, byte_index);
  }

  void read_position_index_entry(PositionIndex& index_element, off_t& byte_index){
    pread_value_from_binary(this->sequence_file_descriptor, index_element.position, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.num_kmers, byte_index);
    index_element.kmer_indexes.reserve(index_element.num_kmers);
    for (uint64_t i=0; i < index_element.num_kmers; i++){
      KmerIndex ki;
      this->read_kmer_index_entry(ki, byte_index);
      auto ret = index_element.kmer_indexes.emplace(ki.name, move(ki));
      throw_assert(ret.second,
                   "ERROR: possible duplicate kmer (" + ret.first->first + ") at position (" +
                   to_string(index_element.position) + ") found in events file: " +
                       this->sequence_file_path)
    }
  }

  void read_contig_strand_index_entry(ContigStrandIndex& index_element, off_t& byte_index){
    // read string length, Contig, strand and nanopore strand
    pread_value_from_binary(this->sequence_file_descriptor, index_element.contig_string_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.contig, index_element.contig_string_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.strand, 1, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.nanopore_strand, 1, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.num_positions, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.num_written_positions, byte_index);
    index_element.position_indexes.reserve(index_element.num_written_positions);
    for (uint64_t i=0; i < index_element.num_written_positions; i++){
      PositionIndex pi;
      this->read_position_index_entry(pi, byte_index);
      auto ret = index_element.position_indexes.emplace(pi.position, move(pi));
      throw_assert(ret.second,
          "ERROR: possible duplicate position (" + to_string(ret.first->first) + ") found in events file: " +
          this->sequence_file_path)
    }
  }

};

#endif //EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_
