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
  EventDataHandler(ReferenceHandler &reference, bool two_d=false) :
      two_d(move(two_d))
  {
    this->initialize_reference(reference, two_d);
  }
  EventDataHandler(const string& event_file, ReferenceHandler &reference, bool two_d=false) :
      two_d(move(two_d))
  {
    this->initialize_reference(reference, two_d);
    this->initialize_reader(event_file);

  }
  EventDataHandler() {}
  ~EventDataHandler() = default;

  string event_file;
  unordered_map<string, ContigStrand> data;
  bool two_d = false;
  BinaryEventReader reader;

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
    this->event_file = internal_event_file;
    reader.initialize(internal_event_file);
  }

  void add_kmer_event(const string& contig, const string& strand, const string& nanopore_strand, const uint64_t& reference_index,
                      const string& path_kmer, const float& descaled_event_mean, const float& posterior_probability){
    data.at(contig+strand+nanopore_strand).get_position(reference_index).soft_add_kmer_event(path_kmer, descaled_event_mean, posterior_probability);
  }

  void write_to_file(path& output_file){
    BinaryEventWriter bew(output_file);
    for (auto cs_pair: data){
      bew.write_contig_strand(cs_pair.second);
    }
    bew.write_indexes();
  }

  Kmer& get_kmer(const string& contig, const string& strand, const string& nanopore_strand,
      const uint64_t& reference_index, const string& path_kmer){
    string contig_strand = contig+strand+nanopore_strand;
    auto found = data.find(contig_strand);
    if (found != data.end()) {
      // Found it
      return data[contig_strand].get_position(reference_index).get_kmer(path_kmer);
    } else {
      if (reader.initialized){
        reader.get_kmer(data[contig_strand].get_position(reference_index).kmers[path_kmer],
            path_kmer, contig, strand, reference_index, nanopore_strand);
        return data[contig_strand].get_position(reference_index).get_kmer(path_kmer);
      }
      // Not there
      throw runtime_error(contig_strand + " was not found in data and there is no ");
    }
  }

  Position& get_position(const string& contig, const string& strand, const string& nanopore_strand,
      const uint64_t& reference_index){
    string contig_strand = contig+strand+nanopore_strand;
    auto found = data.find(contig_strand);
    if (found != data.end()) {
      // Found it
      return data[contig_strand].get_position(reference_index);
    } else {
      if (reader.initialized){
        reader.get_position(data[contig_strand].get_position(reference_index), contig, strand, reference_index, nanopore_strand);
        return data[contig_strand].get_position(reference_index);
      }
      // Not there
      throw runtime_error(contig_strand + " was not found in data and there is no ");
    }
  }


};

#endif //EMBED_FAST5_SRC_EVENTDATAHANDLER_HPP_
