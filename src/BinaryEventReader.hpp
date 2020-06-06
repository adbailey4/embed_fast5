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
    if (this->initialized){
      ::close(sequence_file_descriptor);
    }
  }

  shared_ptr<PosKmer> get_position_kmer(const string& kmer, const string& contig_name, const string& strand,
                         const uint64_t& position, const string nanopore_strand= "t"){
    shared_ptr<PosKmerIndex> ki = this->get_pos_kmer_index(contig_name, strand, nanopore_strand, position, kmer);
    shared_ptr<PosKmer> kmer_struct = make_shared<PosKmer>();
    get_position_kmer(kmer_struct, ki);
    return kmer_struct;
  }

  void get_position_kmer(shared_ptr<PosKmer> kmer_struct, shared_ptr<PosKmerIndex> ki){
    kmer_struct->kmer = ki->name;
    off_t byte_index = ki->sequence_byte_index;
    uint64_t sequence_length = ki->sequence_length;
    kmer_struct->events.reserve(sequence_length);
    pread_vector_from_binary(this->sequence_file_descriptor, kmer_struct->events, sequence_length, byte_index);
  }


  void get_position(Position& position_struct, const string& contig_name, const string& strand,
                const uint64_t& position, const string nanopore_strand="t"){
    PositionIndex& pi = this->get_position_index(contig_name, strand, nanopore_strand, position);
    vector<string> kmers= pi.get_kmers();
    for (auto &kmer: kmers){
      if (!position_struct.has_kmer(kmer)){
        shared_ptr<PosKmer> k = get_position_kmer(kmer, contig_name, strand, position, nanopore_strand);
        position_struct.add_kmer(k);
      }
      position_struct.populated = true;
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

  shared_ptr<PosKmerIndex> get_pos_kmer_index(const string& contig, const string& strand, const string& nanopore_strand, const uint64_t& position, const string& kmer){
    return this->get_position_index(contig, strand, nanopore_strand, position).get_kmer_index(kmer);
  }

  KmerIndex& get_kmer_index(const string& kmer){
    return kmer_map.get_kmer_index(kmer);
  }

  void populate_kmer(Kmer& kmer){
    KmerIndex& kmer_index = this->get_kmer_index(kmer.kmer);
    uint64_t size = kmer_index.kmer_index_ptrs.size();
    kmer.pos_kmer_map.reserve(size);
//    kmer.positions.reserve(size);
//    kmer.contig_strands.reserve(size);
//    vector<shared_ptr<PosKmer>> kmer_struct(size);
    for (uint64_t i = 0; i < size; ++i){
      uint64_t& pos = kmer_index.positions[i];
      string& contig_strand = kmer_index.contig_strands[i];
      //      look for pos and contig strand
      if (!kmer.has_pos_kmer(contig_strand, pos)) {
//      if both are not found then
        shared_ptr<PosKmer> ptr = make_shared<PosKmer>();
        get_position_kmer(ptr, kmer_index.kmer_index_ptrs[i]);
        kmer.add_pos_kmer(contig_strand, pos, ptr);
      }
    }
  }


  /// Attributes ///
  unordered_map<string, ContigStrandIndex> indexes;
  ByKmerIndex kmer_map;
  uint64_t kmer_length;
  uint64_t alphabet_length;
  set<char> alphabet;
  string alphabet_string;
  bool rna = false;
  bool two_d = false;

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
    pread_value_from_binary(this->sequence_file_descriptor, this->kmer_length, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, this->alphabet_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, this->alphabet_string, this->alphabet_length, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, this->rna, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, this->two_d, byte_index);

    this->alphabet = string_to_char_set(this->alphabet_string);
//    read in all other indexes
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
  
  void read_kmer_index_entry(PosKmerIndex& index_element, off_t& byte_index){
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_byte_index, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_length, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.name_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.name, index_element.name_length, byte_index);
  }

  void read_position_index_entry(PositionIndex& index_element, off_t& byte_index, const string& contig_strand){
    pread_value_from_binary(this->sequence_file_descriptor, index_element.position, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.num_kmers, byte_index);
    index_element.kmer_indexes.reserve(index_element.num_kmers);
    for (uint64_t i=0; i < index_element.num_kmers; i++){
      std::shared_ptr<PosKmerIndex> p = std::make_shared<PosKmerIndex>();
      this->read_kmer_index_entry(*p, byte_index);
      kmer_map.add_kmer_index_ptr(contig_strand, index_element.position, p);
      auto ret = index_element.kmer_indexes.emplace(p->name, p);
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
      this->read_position_index_entry(pi, byte_index,
          index_element.contig+index_element.strand+index_element.nanopore_strand);
      auto ret = index_element.position_indexes.emplace(pi.position, move(pi));
      throw_assert(ret.second,
          "ERROR: possible duplicate position (" + to_string(ret.first->first) + ") found in events file: " +
          this->sequence_file_path)
    }
  }
};

#endif //EMBED_FAST5_SRC_BINARYEVENTREADER_HPP_
