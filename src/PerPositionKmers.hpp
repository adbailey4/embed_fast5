//
// Created by Andrew Bailey on 2019-07-24.
//

#ifndef EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_
#define EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_

// embed source
#include "ReferenceHandler.hpp"
#include "AlignmentFile.hpp"
#include "EmbedUtils.hpp"
// boost lib
#include <boost/filesystem.hpp>
#include <boost/heap/priority_queue.hpp>
// std lib
#include <map>
#include <fstream>
#include <mutex>
#include <unordered_map>

using namespace std;
using namespace boost::filesystem;

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

struct Kmer {
  string kmer;
  boost::heap::priority_queue<Event> events;
  uint64_t max_kmers;
  bool max_kmer_set = false;

  Kmer(string kmer) :
      kmer(move(kmer)) {}
  Kmer(string kmer, uint64_t max_kmers) :
      kmer(move(kmer)), max_kmers(move(max_kmers)), max_kmer_set(true) {
//    add one for extra when adding extras
    events.reserve(max_kmers + 1);
  }
  Kmer() {}
  ~Kmer() = default;
  void add_event(Event event) {
    if (!max_kmer_set) {
      events.push(event);
    } else {
      if (events.size() < this->max_kmer_set) {
//    add to queue if not at capacity
        events.push(event);
//    if at capacity check to see if prob is greater than min
      } else if (events.top().posterior_probability < event.posterior_probability) {
        events.push(event);
        while (events.size() > this->max_kmer_set) {
          events.pop();
        }
      }
    }
  }
};

struct Position
{
  uint64_t position;
  unordered_map<string, Kmer> kmers;
  Position(uint64_t position) :
      position(move(position))
  {}
  Position() {}
  ~Position() = default;
  void add_kmer(Kmer kmer){
    auto search = kmers.find(kmer.kmer);
    throw_assert(search != kmers.end(), "Kmer: " + kmer.kmer + " is already in Position.")
    kmers.insert({kmer.kmer, kmer});
  }
  void add_event(string& kmer, Event event){
    auto search = kmers.find(kmer);
    if (search != kmers.end()) {
      kmers[kmer].add_event(event);
    } else {
      Kmer kmer1(kmer);
      kmer1.add_event(event);
      kmers.insert({kmer, kmer1});
    }
  }
  void soft_add_kmer_event(Kmer kmer, Event event){
    auto search = kmers.find(kmer.kmer);
    if (search != kmers.end()) {
      kmers[kmer.kmer].add_event(event);
    } else {
      kmer.add_event(event);
      kmers[kmer.kmer]= kmer;
    }
  }
  void add_kmer_event(Kmer kmer, Event event){
    auto search = kmers.find(kmer.kmer);
    throw_assert(search != kmers.end(), "Kmer: " + kmer.kmer + " is already in Position.")
    kmer.add_event(event);
    kmers.insert({kmer.kmer, kmer});
  }

  Kmer& get_kmer(string& kmer){
    auto search = kmers.find(kmer);
    throw_assert(search == kmers.end(), "Kmer: " + kmer + " is not in Position.")
    return kmers[kmer];
  }
};

struct ContigStrand
{
  string contig_strand;
  vector<Position> positions;
  uint64_t num_positions;
  uint64_t pos_offset;
  uint64_t name_length;
  ContigStrand() {}

  ContigStrand(string contig_strand, uint64_t num_positions, uint64_t pos_offset) :
      contig_strand(move(contig_strand)), num_positions(move(num_positions)), pos_offset(move(pos_offset))
  {
    name_length = contig_strand.length();
    this->initialize_positions();
  }
  ~ContigStrand() = default;

  void initialize_positions(){
    positions.resize(num_positions);
  }

  Position& get_position(uint64_t position){
    return positions[position+pos_offset];
  }

  string get_contig(){
    return contig_strand.substr (0,name_length-1);
  }

  string get_strand(){
    return contig_strand.substr(name_length-1, 1);
  }


};

class PerPositonKmers {
 public:
  PerPositonKmers() {}

  PerPositonKmers(uint64_t num_locks, path output_file, ReferenceHandler &reference) :
      output_file(std::move(output_file)), num_locks(num_locks), initialized_locks(true)
  {
    this->initialize_locks();
    this->initialize_reference(reference);
  }

  ~PerPositonKmers() = default;
  path output_file;
  //  deal with locks for multi threading
  uint64_t num_locks;
  bool initialized_locks;
  //  initialize internal dataset
  unordered_map<string, ContigStrand> data;

  void initialize_locks() {
    locks = std::vector<std::mutex>(this->num_locks);
  }
  void initialize_locks(uint64_t number_of_locks) {
    this->num_locks = number_of_locks;
    locks = std::vector<std::mutex>(this->num_locks);
  }

  void initialize_reference(ReferenceHandler &reference) {
    vector<string> chr_names = reference.get_chromosome_names();
    vector<string> strands{"+", "-"};
    uint64_t length;
    string contig_strand;
    for (auto &contig: chr_names){
      for (auto &strand: strands){
        contig_strand = contig+strand;
        length = reference.get_chromosome_sequence_length(contig);
        data[contig_strand] = ContigStrand(contig_strand, length, 0);
      }
    }
  }

  //  read in alignment file data
  void process_alignment(AlignmentFile &af) {
    char read_strand = af.strand.c_str()[0];
    string contig_strand = af.contig + af.strand;
    for (auto &event: af.iterate()){
      //    lock position
      std::unique_lock<std::mutex> lk(this->locks[event.reference_index % this->num_locks]);
      data[contig_strand].get_position(event.reference_index).add_event(event.path_kmer,
                                                                        Event(event.descaled_event_mean, event.posterior_probability));
//    unlock position
      lk.unlock();
    }
  }
 private:
//  mutexes
  std::vector<mutex> locks;

};

#endif //EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_
