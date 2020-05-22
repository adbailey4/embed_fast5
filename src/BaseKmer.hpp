//
// Created by Andrew Bailey on 5/20/20.
//

#ifndef EMBED_FAST5_SRC_BASEKMER_HPP_
#define EMBED_FAST5_SRC_BASEKMER_HPP_

// embed lib
#include "EmbedUtils.hpp"
// stdlib
#include <unordered_map>

using namespace std;
using namespace embed_utils;

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
Create comparison function for min heap
@param a: first eventkmer
@param b: second eventkmer

@return a.posterior_probability > b.posterior_probability
*/
static bool event_greater_than(const Event& a, const Event& b) {
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
      std::push_heap(events.begin(), events.end(), event_greater_than);
    } else {
      if (events.size() < this->max_events) {
//    add to queue if not at capacity
        events.emplace_back(descaled_event_mean, posterior_probability);
        std::push_heap(events.begin(), events.end(), event_greater_than);

//    if at capacity check to see if prob is greater than min
      } else if (events.front().posterior_probability < posterior_probability) {
        events.emplace_back(descaled_event_mean, posterior_probability);
        std::push_heap(events.begin(), events.end(), event_greater_than);

        while (events.size() > this->max_events) {
          std::pop_heap(events.begin(), events.end(), event_greater_than);
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
  void soft_add_kmer_event(string kmer, Event& event){
    soft_add_kmer_event(move(kmer), event.descaled_event_mean, event.posterior_probability);
  }
  /**
  Add event to kmer or create and add kmer data structure

  @param kmer: kmer string
  @param event: Event data structure
  */
  void soft_add_kmer_event(string kmer, const float &descaled_event_mean, const float &posterior_probability){
    auto search = kmers.find(kmer);
    if (search != kmers.end()) {
      kmers[kmer].add_event(descaled_event_mean, posterior_probability);
    } else {
      Kmer kmer1(kmer);
      kmer1.add_event(descaled_event_mean, posterior_probability);
      kmers.insert({move(kmer), move(kmer1)});
    }
    has_data = true;
  }

  /**
  Return reference to Kmer object

  @param kmer: kmer string
  */
  Kmer& get_kmer(const string& kmer){
    throw_assert(kmers.find(kmer) != kmers.end(), "Kmer: " + kmer + " is not in Position.")
    return kmers[kmer];
  }

  vector<string> get_kmers(){
    vector<string> v;
    for (auto i : kmers) {
      v.push_back(i.first);
    }
    return v;
  }

};

/**
Data structure for keeping track of Position data structures for a specific contig, strand and nanopore strand

@param contig: string name of contig
@param strand: "+" or "-" depending on strand of contig
@param num_positions: number of positions in contig
@param nanopore_strand: template or complement nanopore read ("t" or "c")
*/
struct ContigStrand
{
  string contig;
  string strand;
  vector<Position> positions;
  uint64_t num_positions;
  string nanopore_strand;

  ContigStrand() {}

  ContigStrand(const string& contig, const string& strand, const uint64_t& num_positions, const string& nanopore_strand = "t") :
      contig(move(contig)), strand(move(strand)), num_positions(move(num_positions)),
      nanopore_strand(move(nanopore_strand))
  {
    this->initialize_positions();
  }
  ~ContigStrand() = default;
  /**
  Initialize the number of positions in vector
  */
  void initialize_positions(){
    positions.reserve(num_positions);
    for (uint64_t i=0; i < num_positions; i++) {
      positions.emplace_back(i);
    }
  }
  /**
  Return corresponding Position
  @param position: if position is based off of position offset
  */
  Position& get_position(const uint64_t &position){
    throw_assert(position < num_positions, "Reference Position " + to_string(position) +" not found in " + contig+nanopore_strand+strand)
    return positions.at(position);
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
    throw_assert(position < num_positions,
                 "Position is out of range of initialized values: query (pos): "
                     + to_string(position) + " num_positions: "+ to_string(num_positions));
    positions[position].add_kmer(kmer);
  }
  /**
  Add event to kmer at a position
  @param position: 0 based position
  @param kmer: kmer string
  @param descaled_event_mean: descaled_event_mean
  @param posterior_probability: posterior_probability
  */
  void add_event(uint64_t& position, string& kmer, float &descaled_event_mean, float &posterior_probability){
    throw_assert(position < num_positions,
                 "Position is out of range of initialized values: query (pos): "
                     + to_string(position) + " num_positions: "+ to_string(num_positions));
    positions[position].soft_add_kmer_event(kmer, descaled_event_mean, posterior_probability);
  }

  void add_event(uint64_t& position, string& kmer, Event& event){
    throw_assert(position < num_positions,
                 "Position is out of range of initialized values: query (pos): "
                     + to_string(position) + " num_positions: "+ to_string(num_positions));
    positions[position].soft_add_kmer_event(kmer, event);
  }
};



#endif //EMBED_FAST5_SRC_BASEKMER_HPP_
