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

//ostream& operator<<(ostream& s, KmerIndex& index);


class BinaryEventWriter {
 private:
  /// Attributes ///
  path sequence_file_path;
  std::ofstream sequence_file;
  // When writing the binary file, this vector is appended, so the position of each sequence is stored
  vector<ContigStrandIndex> contig_strand_indexes;
  set<char> alphabet;
  uint64_t kmer_len;
  bool rna;
  bool two_d;

  PositionIndex write_position(Position& position){
    PositionIndex index;
    // Store the number of kmers
    index.num_kmers = position.num_kmers();
    // Store the name of this sequence
    index.position = position.position;
    for (auto &kmer_ptr: position.get_kmer_pointers()){
      index.kmer_indexes.insert(std::make_pair(kmer_ptr->kmer, this->write_kmer(kmer_ptr)));
    }
    return index;
  }

  shared_ptr<PosKmerIndex> write_kmer(shared_ptr<PosKmer> kmer){
    if (kmer->events.empty()){
      throw runtime_error("ERROR: empty sequence provided to BinaryEventWriter: " + kmer->kmer);
    }
    shared_ptr<PosKmerIndex> index = make_shared<PosKmerIndex>();
    // Add sequence start position to index
    index->sequence_byte_index = this->sequence_file.tellp();
    // Store the length of this sequence
    index->sequence_length = kmer->events.size();
    // Store the name of this sequence
    index->name = kmer->kmer;
    write_vector_to_binary(this->sequence_file, kmer->events);
    // Append index object to vector
    return index;
  }


  void write_kmer_index(PosKmerIndex& index){
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
      write_kmer_index(*kmer_index.second);
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
  BinaryEventWriter(path file_path,
                    const set<char> &alphabet,
                    const uint64_t &kmer_len,
                    bool rna,
                    bool two_d) :
      alphabet(move(alphabet)), kmer_len(move(kmer_len)), rna(move(rna)),
      two_d(move(two_d))
  {
    this->sequence_file_path = file_path;
    // Ensure that the output directory exists
    if (this->sequence_file_path.has_parent_path()){
      create_directories(this->sequence_file_path.parent_path());
    }
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
    write_value_to_binary(this->sequence_file, kmer_len);
    // How long is the name of the sequence
    string alphabet_str = char_set_to_string(alphabet);
    write_value_to_binary(this->sequence_file, alphabet_str.size());
    write_string_to_binary(this->sequence_file, alphabet_str);
    write_value_to_binary(this->sequence_file, rna);
    write_value_to_binary(this->sequence_file, two_d);

    // Iterate all the indexes, write them to the file
    for (auto& index: this->contig_strand_indexes){
      write_contig_strand_index(index);
    }
    // Write the pointer to the beginning of the indexes
    write_value_to_binary(this->sequence_file, indexes_start_position);
  }
};

#endif //EMBED_FAST5_SRC_BINARYEVENTWRITER_HPP_
