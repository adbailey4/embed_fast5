//
// Created by Andrew Bailey on 5/20/20.
//

#ifndef EMBED_FAST5_SRC_BASEKMER_HPP_
#define EMBED_FAST5_SRC_BASEKMER_HPP_

// embed lib
#include "EmbedUtils.hpp"
// stdlib
#include <unordered_map>
#include <set>

using namespace std;
using namespace embed_utils;

class PosKmerIndex {
 public:
  /// Attributes ///
  string name;
  uint64_t name_length;
  uint64_t sequence_byte_index;
  uint64_t sequence_length;
};

class KmerIndex {
 public:
  /// Attributes ///
  string kmer;
  vector<shared_ptr<PosKmerIndex>> kmer_index_ptrs;
//  unordered_map<uint64_t, shared_ptr<PosKmer>> kmers;
  vector<uint64_t> positions;
  vector<string> contig_strands;
  KmerIndex() = default;
  ~KmerIndex() = default;

  void add_pos_kmer_index(string contig_strand, uint64_t position, shared_ptr<PosKmerIndex> kmer_index_ptr){
    contig_strands.push_back(contig_strand);
    positions.push_back(position);
    kmer_index_ptrs.push_back(kmer_index_ptr);
  }

};


class ByKmerIndex {
 public:
  /// Attributes ///
  ByKmerIndex() = default;
  ~ByKmerIndex() = default;
  void add_kmer_index_ptr(const string& contig_strand, const uint64_t& position, shared_ptr<PosKmerIndex> kmer){
    kmer_index_map[kmer->name].add_pos_kmer_index(contig_strand, position, kmer);
  }

  KmerIndex& get_kmer_index(const string& kmer){
    auto found = kmer_index_map.find(kmer);
    if (found != kmer_index_map.end()) {
      // Found it
      return found->second;
    } else {
      throw runtime_error(kmer + " was not found in kmer index map");
    }
  }

 private:
  unordered_map<string, KmerIndex> kmer_index_map;

};


class PositionIndex {
 public:
  /// Attributes ///
  uint64_t position;
  unordered_map<string, shared_ptr<PosKmerIndex>> kmer_indexes;
  uint64_t num_kmers;

  shared_ptr<PosKmerIndex> get_kmer_index(const string& kmer){
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
    for (auto &i : kmer_indexes) {
      v.push_back(i.first);
    }
    return v;
  }

};


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
struct PosKmer {
  string kmer;
  vector<Event> events;
  uint64_t max_events;
  bool max_events_set = false;
  uint64_t num_ignored_kmers;

  PosKmerIndex index;

  PosKmer(string kmer) :
      kmer(move(kmer)) {}
  PosKmer(string kmer, uint64_t max_events) :
      kmer(move(kmer)), max_events(move(max_events)), max_events_set(true) {
//    add one for extra when adding extras
    events.reserve(max_events + 1);
  }
  PosKmer() {}
  ~PosKmer() = default;
//  Kmer(const Kmer& mE)            = default;
//  Kmer(Kmer&& mE)                 = default;
  PosKmer(const PosKmer& other) :
      kmer{other.kmer},
      events{other.events},
      max_events{other.max_events},
      max_events_set{other.max_events_set},
      num_ignored_kmers{other.num_ignored_kmers}        {
//    cout << kmer + ": KMER COPY CONSTRUCTOR" << '\n';
  }

//  Kmer(Kmer&& mE) : a{move(mE.a)}, b{move(mE.b)} { }
  PosKmer& operator=(const PosKmer& other) {
    if(this != &other) {
      kmer = other.kmer;
      events = other.events;
      max_events = other.max_events;
      max_events_set = other.max_events_set;
      num_ignored_kmers = other.num_ignored_kmers;
//      cout << kmer + ": KMER COPY ASSIGNMENT" << '\n';
    }
    return *this;
  }

  PosKmer(PosKmer&& other) :
      kmer{move(other.kmer)},
      events{move(other.events)},
      max_events{move(other.max_events)},
      max_events_set{move(other.max_events_set)},
      num_ignored_kmers{move(other.num_ignored_kmers)} {
//    cout << kmer + ": KMER MOVE CONSTRUCTOR" << '\n';
  }

  PosKmer& operator=(PosKmer&& other) {
    if(this != &other) {
      kmer = move(other.kmer);
      events = move(other.events);
      max_events = move(other.max_events);
      max_events_set = move(other.max_events_set);
      num_ignored_kmers = move(other.num_ignored_kmers);
//      cout << kmer + ": KMER MOVE ASSIGNMENT" << '\n';
    }
    return *this;
  }


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

struct Kmer {
 public:
  string kmer;
  vector<shared_ptr<PosKmer>> kmers;
//  unordered_map<uint64_t, shared_ptr<PosKmer>> kmers;
  vector<uint64_t> positions;
  vector<string> contig_strands;
  Kmer(string kmer) :
    kmer(move(kmer)) {}
  ~Kmer() = default;

  void add_pos_kmer(const string& contig_strand, const uint64_t& position, shared_ptr<PosKmer> kmer_index_ptr){
    contig_strands.push_back(contig_strand);
    positions.push_back(position);
    kmers.push_back(kmer_index_ptr);
  }

  void soft_add_pos_kmer(const string& contig_strand, const uint64_t& position, shared_ptr<PosKmer> kmer_index_ptr){
    if (find_index(contig_strand, position) == -1){
      add_pos_kmer(contig_strand, position, kmer_index_ptr);
    }
  }

  int64_t find_index(const string& contig_strand, const uint64_t& pos){
    int64_t counter = 0;
    for (auto &cs: contig_strands) {
      if (cs != contig_strand) {
        counter += 1;
      } else {
        break;
      }
    }
    if (counter == contig_strands.size()){
      return -1;
    }
    int64_t counter2 = 0;
    for (auto &pos1: positions){
      if (pos1 != pos){
        counter2 += 1;
      } else {
        break;
      }
    }
    if (counter == counter2){
      return counter;
    } else {
      return -1;
    }

  }
};

class ByKmer {
 public:
  /// Attributes ///
  set<char> alphabet;
  uint64_t kmer_length;

  ByKmer(set<char> alphabet, uint64_t kmer_length){
    initialize_kmer_map(alphabet, kmer_length);
  }

  ByKmer() = default;
  ~ByKmer() = default;

  void initialize_kmer_map(set<char> alphabet1, uint64_t kmer_length1){
    alphabet = alphabet1;
    kmer_length = kmer_length1;
    string str_alphabet = char_set_to_string(this->alphabet);
    for (auto &k: all_string_permutations(str_alphabet, this->kmer_length)){
      this->add_kmer(k);
    }
  }

  void add_kmer_ptr(const string& contig_strand, const uint64_t& position, shared_ptr<PosKmer> kmer){
    this->get_kmer(kmer->kmer).soft_add_pos_kmer(contig_strand, position, kmer);
  }

  void add_kmer(string kmer){
    kmer_map.emplace(kmer, kmer);
  }

  Kmer& get_kmer(const string& kmer){
    throw_assert(this->has_kmer(kmer),
        "Kmer: " + kmer + " is not in kmer map with alphabet=" + char_set_to_string(alphabet) + " and kmer length="+ to_string(kmer_length))
    return kmer_map.at(kmer);
  }

  bool has_kmer(const string& kmer){
    return kmer_map.find(kmer) != kmer_map.end();
  }

 private:
  unordered_map<string, Kmer> kmer_map;

};


/**
Data structure for keeping track of kmers aligned to a position

@param position: reference position
*/
struct Position
{
 private:
  unordered_map<string, shared_ptr<PosKmer>> kmers;
 public:
  uint64_t position;
  bool has_data = false;
  bool populated = false;
  Position(uint64_t position) :
      position(move(position))
  {}
  Position() {}
  ~Position() = default;

  uint64_t num_kmers(){
    return kmers.size();
  }

  /**
  Add Kmer data structure to internal map via a move

  @param kmer: Kmer structure
  */
  void add_kmer(shared_ptr<PosKmer> kmer){
    throw_assert(kmers.find(kmer->kmer) == kmers.end(), "Kmer: " + kmer->kmer + " is already in Position.")
    string kmer_str(kmer->kmer);
    kmers.emplace(kmer_str, move(kmer));
//    string kmer_str(kmer.kmer);
//    kmers.insert(std::make_pair<string, Kmer>(static_cast<basic_string<char> &&>(kmer.kmer), forward<Kmer>(kmer)));
    has_data = true;
  }
  /**
  Add event to kmer or create and add kmer data structure

  @param kmer: kmer string
  @param event: Event data structure
  */
  void soft_add_kmer_event(string kmer, Event& event){
    soft_add_kmer_event(kmer, event.descaled_event_mean, event.posterior_probability);
  }
  /**
  Add event to kmer or create and add kmer data structure

  @param kmer: kmer string
  @param event: Event data structure
  */
  void soft_add_kmer_event(string kmer, const float &descaled_event_mean, const float &posterior_probability){
    auto search = kmers.find(kmer);
    if (search != kmers.end()) {
      kmers[kmer]->add_event(descaled_event_mean, posterior_probability);
    } else {
      shared_ptr<PosKmer> kmer1 = make_shared<PosKmer>(kmer);
      kmer1->add_event(descaled_event_mean, posterior_probability);
      kmers.insert({move(kmer), move(kmer1)});
    }
    has_data = true;
  }

  /**
  Return shared pointer to PosKmer pointer
  @param kmer: kmer string
  */
  shared_ptr<PosKmer> get_kmer(const string& kmer){
    throw_assert(kmers.find(kmer) != kmers.end(), "Kmer: " + kmer + " is not in Position.")
    return kmers[kmer];
  }

  /**
  Return true if kmer is in map, false if not
  @param kmer: kmer string
  */
  bool has_kmer(const string& kmer){
    return kmers.find(kmer) != kmers.end();
  }

  vector<string> get_kmer_strings(){
    vector<string> v;
    for (auto i : kmers) {
      v.push_back(i.first);
    }
    return v;
  }

  vector<shared_ptr<PosKmer>> get_kmer_pointers(){
    vector<shared_ptr<PosKmer>> v;
    for (auto i : kmers) {
      v.push_back(i.second);
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
//  PosKmerIndex& get_kmer_index(const uint64_t& position, const string& kmer) {
//    return get_position_index(position).get_kmer_index(kmer);
//  }
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
  void add_kmer(const uint64_t &position, shared_ptr<PosKmer> kmer){
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
