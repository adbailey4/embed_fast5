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

struct Kmer
{
  string kmer;
  vector<Event> events;
  Kmer(string kmer) :
      kmer(move(kmer))
  {}
  Kmer() {}
  ~Kmer() = default;
  void add_event(Event event){
    events.push_back(event);
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
    kmers[kmer.kmer]= kmer;
  }
  void add_event(string& kmer, Event event){
    auto search = kmers.find(kmer);
    if (search != kmers.end()) {
      kmers[kmer].add_event(event);
    } else {
      Kmer kmer1(kmer);
      kmer1.add_event(event);
      kmers[kmer1.kmer]= kmer1;
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
  void add_event(Kmer kmer, Event event){
    auto search = kmers.find(kmer.kmer);
    if (search != kmers.end()) {
      kmers[kmer.kmer].add_event(event);
    } else {
      kmer.add_event(event);
      kmers[kmer.kmer]= kmer;
    }
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
  PerPositonKmers(uint64_t num_locks, path output_file, ReferenceHandler &reference);
  ~PerPositonKmers();
  path output_file;
  //  deal with locks for multi threading
  uint64_t num_locks;
  bool initialized_locks;
  void initialize_locks();
  void initialize_locks(uint64_t number_of_locks);
//  initialize internal dataset
  unordered_map<string, ContigStrand> data;
  void initialize_reference(ReferenceHandler &reference);
//  read in alignment file data
  void process_alignment(AlignmentFile &af);
 private:
//  mutexes
  std::vector<mutex> locks;

};

#endif //EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_
