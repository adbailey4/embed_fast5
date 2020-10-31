//
// Created by Andrew Bailey on 5/20/20.
//

#ifndef EMBED_FAST5_SRC_EVENTDATAHANDLER_HPP_
#define EMBED_FAST5_SRC_EVENTDATAHANDLER_HPP_

// embed libs
#include "ReferenceHandler.hpp"
#include "BinaryEventReader.hpp"

using namespace std;

class EventDataHandler {
 public:
  EventDataHandler(ReferenceHandler &reference,
                   set<char> alphabet,
                   uint64_t kmer_length,
                   bool two_d = false,
                   bool rna = false)
  {
    this->initialize(reference, alphabet, kmer_length, two_d, rna);
  }
  EventDataHandler(ReferenceHandler &reference, const string &event_file)
  {
    this->initialize(reference, event_file);
  }

  EventDataHandler() {}
  ~EventDataHandler() = default;

  set<char> get_alphabet(){
    return alphabet;
  }

  uint64_t get_kmer_length(){
    return kmer_length;
  }

  void initialize(ReferenceHandler &reference1,
                  set<char> alphabet1,
                  uint64_t kmer_length1,
                  bool two_d1 = false,
                  bool rna1 = false) {
    alphabet = alphabet1;
    kmer_length = kmer_length1;
    two_d = two_d1;
    rna = rna1;
    this->initialize_reference(reference1, two_d);
    this->initialize_by_kmer();
  }

  void initialize(ReferenceHandler &reference1, const string &event_file1){
    this->initialize_reader(event_file1);
    this->initialize_reference(reference1, two_d);
    this->initialize_by_kmer();
  }

  //  initialize internal data structure
  void initialize_reference(ReferenceHandler &reference, const bool& internal_two_d=false) {
    vector<string> chr_names = reference.get_chromosome_names();
    vector<string> strands{"+", "-"};
    vector<string> nanopore_strands{"t", "c"};
    for (auto &contig: chr_names){
      for (auto &strand: strands){
        if (internal_two_d) {
          for (auto &nanopore_strand: nanopore_strands) {
            data.emplace(std::piecewise_construct,
                        std::forward_as_tuple(contig + strand + nanopore_strand),
                        std::forward_as_tuple(contig, strand, reference.get_chromosome_sequence_length(contig), nanopore_strand));

          }
        } else {
          data.emplace(std::piecewise_construct,
                       std::forward_as_tuple(contig + strand + "t"),
                       std::forward_as_tuple(contig, strand, reference.get_chromosome_sequence_length(contig), "t"));
        }
      }
    }
  }

  void initialize_reader(const string& internal_event_file){
    event_file = internal_event_file;
    reader.initialize(event_file);
    alphabet = reader.alphabet;
    kmer_length = reader.kmer_length;
    two_d = reader.two_d;
    rna = reader.rna;
  }

  void initialize_by_kmer(){
    by_kmer_data.initialize_kmer_map(alphabet, kmer_length);
  }

  void add_kmer_event(const string& contig, const string& strand, const string& nanopore_strand, const uint64_t& reference_index,
                      const string& path_kmer, const float& descaled_event_mean, const float& posterior_probability){
    string contig_strand = contig+strand+nanopore_strand;
    Position& pos = data.at(contig_strand).get_position(reference_index);
    pos.soft_add_kmer_event(path_kmer, descaled_event_mean, posterior_probability);
    by_kmer_data.add_kmer_ptr(contig_strand, reference_index, pos.get_pos_kmer(path_kmer));
  }

  void write_to_file(path& output_file){
    throw_assert(kmer_length != (uint64_t)-1, "Kmer length must be set in order to write to file");
    throw_assert(alphabet != set<char>{}, "Alphabet must be set in order to write to file")
    BinaryEventWriter bew(output_file, alphabet, kmer_length, rna, two_d);
    for (auto &cs_pair: data){
      bew.write_contig_strand(cs_pair.second);
    }
    bew.write_indexes();
  }

  shared_ptr<PosKmer> get_position_kmer(const string& contig, const string& strand, const string& nanopore_strand,
                             const uint64_t& reference_index, const string& path_kmer){
    string contig_strand = contig+strand+nanopore_strand;
    throw_assert(data.find(contig_strand) != data.end(),
        "contig_strand: " + contig_strand + " is not in EventDataHandler.")
    Position &pos = data.at(contig_strand).get_position(reference_index);
    if (pos.has_kmer(path_kmer)) {
      return pos.get_pos_kmer(path_kmer);
    }
    return get_position_kmer_from_reader(contig, strand, nanopore_strand, reference_index, path_kmer);
  }

  Kmer& get_kmer(const string& path_kmer){
    Kmer& kmer_data = by_kmer_data.get_kmer(path_kmer);
    uint64_t expected_n_kmers = reader.get_kmer_index(path_kmer).positions.size();
    if (kmer_data.pos_kmer_map.size() != expected_n_kmers) {
      get_kmer_from_reader(kmer_data);
    }
    return kmer_data;
  }

  bool has_kmer(const string& path_kmer){
    if (by_kmer_data.has_kmer(path_kmer)){
      return reader.has_kmer_index(path_kmer);
    } else {
      return false;
    }
  }

  Position& get_position(const string& contig, const string& strand, const string& nanopore_strand,
      const uint64_t& reference_index){
    string contig_strand = contig+strand+nanopore_strand;
    throw_assert(data.find(contig_strand) != data.end(),
                 "contig_strand: " + contig_strand + " is not in EventDataHandler.")
    Position& pos = data.at(contig_strand).get_position(reference_index);
    if (reader.initialized & !pos.populated){
      reader.get_position(pos, contig, strand, reference_index, nanopore_strand);
      for (auto &k: pos.get_kmer_strings()){
        by_kmer_data.get_kmer(k).soft_add_pos_kmer(contig_strand, reference_index, pos.get_pos_kmer(k));
      }
    }
    return pos;
  }

  ContigStrand& get_contig_strand(const string& contig, const string& strand, const string& nanopore_strand){
    string contig_strand = contig+strand+nanopore_strand;
    throw_assert(data.find(contig_strand) != data.end(),"contig_strand: " + contig_strand + " is not in EventDataHandler.")
    ContigStrand& cs = data.at(contig_strand);
    for (uint64_t i=0; i < cs.num_positions; ++i) {
      get_position(contig, strand, nanopore_strand, i);
    }
    return cs;
  }

  void add_set_to_alphabet(const set<char>& kmer){
    alphabet.insert(kmer.begin(), kmer.end());
  }

  void add_kmer_to_alphabet(const string& kmer){
    alphabet.insert(kmer.begin(), kmer.end());
  }

  void set_alphabet(set<char> new_alphabet){
    alphabet = new_alphabet;
  }

  void set_kmer_length(uint64_t new_kmer_length){
    kmer_length = new_kmer_length;
  }

 private:
  set<char> alphabet = {};
  uint64_t kmer_length = -1;
  bool two_d = false;
  bool rna;
  string event_file;

  unordered_map<string, ContigStrand> data;
  ByKmer by_kmer_data;
  BinaryEventReader reader;


  shared_ptr<PosKmer> get_position_kmer_from_reader(const string& contig, const string& strand, const string& nanopore_strand,
                                                    const uint64_t& reference_index, const string& path_kmer){
    if (reader.initialized){
      shared_ptr<PosKmer> shared_ptr_pos_kmer = reader.get_position_kmer(path_kmer, contig, strand, reference_index, nanopore_strand);
      string contig_strand = contig+strand+nanopore_strand;
      data.at(contig_strand).get_position(reference_index).add_kmer(shared_ptr_pos_kmer);
      by_kmer_data.add_kmer_ptr(contig_strand, reference_index, shared_ptr_pos_kmer);
      return shared_ptr_pos_kmer;
    }
    // Not there
    throw runtime_error(contig+strand+nanopore_strand + " was not found in data and there is no reader initialized");

  }

  void get_kmer_from_reader(Kmer& kmer){
    if (reader.initialized){
      reader.populate_kmer(kmer);
      string contig_strand;
      uint64_t position;
      for (auto &k: kmer.contig_positions){
        contig_strand = get<0>(k);
        position = get<1>(k);
        Position& pos = data.at(contig_strand).get_position(position);
        if (!pos.has_kmer(kmer.kmer))
          pos.add_kmer(kmer.get_pos_kmer(contig_strand, position));
      }
    } else {
      throw runtime_error(kmer.kmer + " was not found in data and there is no reader initialized");
    }
  }
};

#endif //EMBED_FAST5_SRC_EVENTDATAHANDLER_HPP_
