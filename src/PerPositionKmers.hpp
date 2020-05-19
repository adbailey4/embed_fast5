//
// Created by Andrew Bailey on 2019-07-24.
//

#ifndef EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_
#define EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_

// embed source
#include "ReferenceHandler.hpp"
#include "AlignmentFile.hpp"
#include "EmbedUtils.hpp"
// boost lib
#include <boost/filesystem.hpp>
// std lib
#include <map>
#include <fstream>
#include <mutex>
#include <unordered_map>

using namespace std;
using namespace boost::filesystem;

/**
Simple data structure for kmer data

@param descaled_event_mean: descaled_event_mean
@param posterior_probability: posterior_probability
*/
struct Event
{
  float descaled_event_mean;
  float posterior_probability;
  Event(float descaled_event_mean, float posterior_probability) :
      descaled_event_mean(move(descaled_event_mean)), posterior_probability(move(posterior_probability))
  {}
  Event() {}
  ~Event() = default;
};

/**
Override < operator in order to make the priority queue a min heap
@param a: first eventkmer
@param b: second eventkmer

@return a.posterior_probability > b.posterior_probability
*/
static bool operator<(const Event& a, const Event& b) {
  return a.posterior_probability > b.posterior_probability;
}

/**
Override < operator in order to make the priority queue a min heap
@param a: first eventkmer
@param b: second eventkmer

@return a.posterior_probability > b.posterior_probability
*/
static bool operator==(const Event& a, const Event& b) {
  return a.posterior_probability == b.posterior_probability and a.descaled_event_mean == b.descaled_event_mean;
}


/**
Simple data structure for keeping track of events aligned to a kmer

@param kmer: string representation
@param max_kmers: limit number of events
*/
struct Kmer {
  string kmer;
  vector<Event> events;
  uint64_t max_events;
  bool max_events_set = false;
  uint64_t num_ignored_kmers;

  Kmer(string kmer) :
      kmer(move(kmer)) {}
  Kmer(string kmer, uint64_t max_events) :
      kmer(move(kmer)), max_events(move(max_events)), max_events_set(true) {
//    add one for extra when adding extras
    events.reserve(max_events + 1);
  }
  Kmer() {}
  ~Kmer() = default;

  /**
  Add an event to the priority queue

  @param descaled_event_mean: descaled_event_mean
  @param posterior_probability: posterior_probability
  */
  void add_event(float descaled_event_mean, float posterior_probability) {
    if (!max_events_set) {
      events.emplace_back(descaled_event_mean, posterior_probability);
      std::push_heap(events.begin(), events.end());
    } else {
      if (events.size() < this->max_events_set) {
//    add to queue if not at capacity
        events.emplace_back(descaled_event_mean, posterior_probability);
        std::push_heap(events.begin(), events.end());

//    if at capacity check to see if prob is greater than min
      } else if (events.front().posterior_probability < posterior_probability) {
        events.emplace_back(descaled_event_mean, posterior_probability);
        std::push_heap(events.begin(), events.end());

        while (events.size() > this->max_events_set) {
          std::pop_heap(events.begin(), events.end());
          events.pop_back();
        }
      } else {
        num_ignored_kmers += 1;
      }
    }
  }
  /**
  Add an event to the priority queue

  @param Event: event already built
  */
  void add_event(Event& event) {
    this->add_event(event.descaled_event_mean, event.posterior_probability);
  }
  /**
  Return the number of events aligned to this kmer

  @return: uint64_t number of events
  */

  uint64_t num_events(){
    return events.size();
  }

};

/**
Data structure for keeping track of kmers aligned to a position

@param position: reference position
*/
struct Position
{
  uint64_t position;
  unordered_map<string, Kmer> kmers;
  bool has_data = false;
  Position(uint64_t position) :
      position(move(position))
  {}
  Position() {}
  ~Position() = default;

  /**
  Add Kmer data structure to internal map via a move

  @param kmer: Kmer structure
  */
  void add_kmer(Kmer& kmer){
    throw_assert(kmers.find(kmer.kmer) == kmers.end(), "Kmer: " + kmer.kmer + " is already in Position.")
    kmers.insert({kmer.kmer, move(kmer)});
    has_data = true;
  }
  /**
  Add event to kmer or create and add kmer data structure

  @param kmer: kmer string
  @param event: Event data structure
  */
  void soft_add_kmer_event(string& kmer, Event& event){
    soft_add_kmer_event(kmer, event.descaled_event_mean, event.posterior_probability);
  }
  /**
  Add event to kmer or create and add kmer data structure

  @param kmer: kmer string
  @param event: Event data structure
  */
  void soft_add_kmer_event(string& kmer, float &descaled_event_mean, float &posterior_probability){
    auto search = kmers.find(kmer);
    if (search != kmers.end()) {
      kmers[kmer].add_event(descaled_event_mean, posterior_probability);
    } else {
      Kmer kmer1(kmer);
      kmer1.add_event(descaled_event_mean, posterior_probability);
      kmers.insert({kmer, move(kmer1)});
    }
    has_data = true;
  }

  /**
  Return reference to Kmer object

  @param kmer: kmer string
  */
  Kmer& get_kmer(string& kmer){
    throw_assert(kmers.find(kmer) != kmers.end(), "Kmer: " + kmer + " is not in Position.")
    return kmers[kmer];
  }
};

/**
Data structure for keeping track of Position data structures for a specific contig, strand and nanopore strand

@param contig: string name of contig
@param strand: "+" or "-" depending on strand of contig
@param num_positions: number of positions in contig
@param pos_offset: if position is based off of position offset
@param nanopore_strand: template or complement nanopore read ("t" or "c")
*/
struct ContigStrand
{
  string contig;
  string strand;
  vector<Position> positions;
  uint64_t num_positions;
  uint64_t pos_offset;
  string nanopore_strand;

  ContigStrand() {}

  ContigStrand(string contig, string strand, uint64_t num_positions, uint64_t pos_offset, string nanopore_strand="t") :
      contig(move(contig)), strand(move(strand)), num_positions(move(num_positions)),
      pos_offset(move(pos_offset)), nanopore_strand(move(nanopore_strand))
  {
    this->initialize_positions();
  }
  ~ContigStrand() = default;
  /**
  Initialize the number of positions in vector
  */
  void initialize_positions(){
    positions.resize(num_positions);
  }
  /**
  Return corresponding Position
  @param position: if position is based off of position offset
  */
  Position& get_position(uint64_t &position){
    return positions[position+pos_offset];
  }

  string get_contig(){
    return contig;
  }

  string get_strand(){
    return strand;
  }
  string get_nanopore_strand(){
    return nanopore_strand;
  }
  /**
  Move kmer to position
  @param position: 0 based position
  @param kmer: Kmer data structure
  */
  void add_kmer(uint64_t& position, Kmer& kmer){
    throw_assert(position+pos_offset < num_positions,
        "Position is out of range of initialized values: query (pos+offset): "
        + to_string(position+pos_offset) + " pos_offset: " + to_string(pos_offset) +
        " num_positions: "+ to_string(num_positions));
    positions[position+pos_offset].add_kmer(kmer);
  }
  /**
  Add event to kmer at a position
  @param position: 0 based position
  @param kmer: kmer string
  @param descaled_event_mean: descaled_event_mean
  @param posterior_probability: posterior_probability
  */
  void add_event(uint64_t& position, string& kmer, float &descaled_event_mean, float &posterior_probability){
    throw_assert(position+pos_offset < num_positions,
                 "Position is out of range of initialized values: query (pos+offset): "
                     + to_string(position+pos_offset) + " pos_offset: " + to_string(pos_offset) +
                     " num_positions: "+ to_string(num_positions));
    positions[position+pos_offset].soft_add_kmer_event(kmer, descaled_event_mean, posterior_probability);
  }

  void add_event(uint64_t& position, string& kmer, Event& event){
    throw_assert(position+pos_offset < num_positions,
                 "Position is out of range of initialized values: query (pos+offset): "
                     + to_string(position+pos_offset) + " pos_offset: " + to_string(pos_offset) +
                     " num_positions: "+ to_string(num_positions));
    positions[position+pos_offset].soft_add_kmer_event(kmer, event);
  }
};


/**
Class for handling processing of alignment files into the underlying ContigStrand data structure

@param num_locks: number of locks to use to protect position data
@param output_file: path to output file
@param reference: ReferenceHandler object to initialize data structure
@param two_d: option to initialize more data for possible complement or template signalalign data
*/
class PerPositionKmers {
 public:
//  constructors and deconstructors
  PerPositionKmers() {}
  PerPositionKmers(uint64_t num_locks, path output_file, ReferenceHandler &reference, bool two_d=false) :
      output_file(std::move(output_file)), num_locks(move(num_locks)), initialized_locks(true), two_d(move(two_d))
  {
    this->initialize_locks();
    this->initialize_reference(reference);
  }
  ~PerPositionKmers() = default;
//  data
  path output_file;
  uint64_t num_locks;
  bool initialized_locks;
  bool two_d;
  unordered_map<string, ContigStrand> data;
// methods
// initalize locks
  void initialize_locks() {
    locks = std::vector<std::mutex>(this->num_locks);
  }
//  initialize internal data structure
  void initialize_reference(ReferenceHandler &reference) {
    vector<string> chr_names = reference.get_chromosome_names();
    vector<string> strands{"+", "-"};
    uint64_t length;
    uint64_t pos_offset = 0;
    string contig_strand;
    for (auto &contig: chr_names){
      for (auto &strand: strands){
        if (two_d) {
          for (auto &nanopore_strand: {"t", "c"}) {
            contig_strand = contig + strand + nanopore_strand;
            length = reference.get_chromosome_sequence_length(contig);
            data.emplace(std::make_pair(contig_strand, ContigStrand(contig, strand, length, pos_offset, nanopore_strand)));
          }
        } else {
          contig_strand = contig + strand + "t";
          length = reference.get_chromosome_sequence_length(contig);
          data.emplace(std::make_pair(contig_strand, ContigStrand(contig, strand, length, pos_offset)));
        }
      }
    }
  }
  //  read in alignment file data
  void process_alignment(AlignmentFile &af) {
    string contig_strand = af.contig + af.strand;
    for (auto &event: af.iterate()){
      //    lock position
      std::unique_lock<std::mutex> lk(this->locks[event.reference_index % this->num_locks]);
      data[contig_strand+event.strand].get_position(event.reference_index).soft_add_kmer_event(event.path_kmer,
          (float&) event.descaled_event_mean,
          (float&) event.posterior_probability);
//    unlock position
      lk.unlock();
    }
  }
 private:
//  mutexes
  std::vector<mutex> locks;
};

#endif //EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_
