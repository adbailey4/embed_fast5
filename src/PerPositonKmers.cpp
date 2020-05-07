//
// Created by Andrew Bailey on 2019-07-24.
//

#include "PerPositonKmers.hpp"
#include "EmbedUtils.hpp"
#include <boost/heap/priority_queue.hpp>

using namespace std;
using namespace embed_utils;

/**
Destroy locks at destruction of class
*/
PerPositonKmers::~PerPositonKmers() {
}

/**
Initialize locks and heaps and other important data structures
*/
PerPositonKmers::PerPositonKmers(uint64_t num_locks, path output_file, ReferenceHandler &reference) :
    output_file(std::move(output_file)), num_locks(num_locks), initialized_locks(true)
{
  this->initialize_locks();
  this->initialize_reference(reference);
}

/**
Initialize vector of locks
*/
void PerPositonKmers::initialize_locks() {
  locks = std::vector<std::mutex>(this->num_locks);
}

/**
Initialize vector of locks
*/
void PerPositonKmers::initialize_locks(uint64_t number_of_locks) {
  this->num_locks = number_of_locks;
  locks = std::vector<std::mutex>(this->num_locks);
}
void PerPositonKmers::process_alignment(AlignmentFile &af) {
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
void PerPositonKmers::initialize_reference(ReferenceHandler &reference) {
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

