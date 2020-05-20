//
// Created by Andrew Bailey on 5/3/20.
// Most of this code came from Ryan Lorig-Roach. https://github.com/rlorigro/runlength_analysis_cpp
//

#ifndef EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_
#define EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_

// embed libs
#include "BinaryEventWriter.hpp"

using namespace std;
using namespace embed_utils;
using namespace boost::filesystem;

using namespace std;

class BinaryEventReader {
 public:

  /// Methods ///

  // Initialize the class with a file path
  BinaryEventReader(string file_path) {
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
  ~BinaryEventReader() {
    this->close();
  }

  void close() {
    ::close(sequence_file_descriptor);
  }

  void get_kmer(Kmer& kmer_struct, string kmer, string& contig_name, string strand, uint64_t position, string nanopore_strand="t"){
    KmerIndex& ki = indexes[contig_name+strand+nanopore_strand].position_indexes[position].kmer_indexes[kmer];
    kmer_struct.kmer = ki.name;
    off_t byte_index = ki.sequence_byte_index;
    uint64_t sequence_length = ki.sequence_length;
    pread_vector_from_binary(this->sequence_file_descriptor, kmer_struct.events, sequence_length, byte_index);
  }

  unordered_map<string, ContigStrandIndex> indexes;

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
      this->indexes.emplace(index_element.contig + index_element.strand + index_element.nanopore_strand, move(index_element));
      // Update the mapping of read names to their places in the vector of indexes
//    auto element = make_pair(index_element.name, this->indexes.size() - 1);
//    auto success = this->index_map.insert(move(element)).second;
//    if (not success){
//      throw runtime_error("ERROR: possible duplicate read name (" + index_element.name + ") found in runnie file: " + this->sequence_file_path);
//    }
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
      index_element.kmer_indexes.emplace(ki.name, move(ki));
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
      index_element.position_indexes.emplace(pi.position, move(pi));
    }
  }
};

#endif //EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_
