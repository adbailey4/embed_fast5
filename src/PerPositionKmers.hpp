//
// Created by Andrew Bailey on 2019-07-24.
//

#ifndef EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_
#define EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_

// embed source
#include "BaseKmer.hpp"
#include "EventDataHandler.hpp"
#include "AlignmentFile.hpp"
#include "BinaryEventWriter.hpp"
// boost lib
#include <boost/filesystem.hpp>
// std lib
#include <map>
#include <fstream>
#include <mutex>

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
  PerPositionKmers(uint64_t num_locks, ReferenceHandler &reference, bool two_d = false) :
    num_locks(move(num_locks)), initialized_locks(true)
  {
    this->initialize_locks();
    data.initialize_reference(reference, two_d);
  }
  ~PerPositionKmers() = default;
//  data
  uint64_t num_locks;
  bool initialized_locks;
  EventDataHandler data;
// methods
// initalize locks
  void initialize_locks() {
    locks = std::vector<std::mutex>(this->num_locks);
  }
  //  read in alignment file data
  void process_alignment(AlignmentFile &af) {
    for (auto &event: af.iterate()){
      //    lock position
      std::unique_lock<std::mutex> lk(this->locks[event.reference_index % this->num_locks]);
      data.add_kmer_event(event.contig, af.strand, event.strand, event.reference_index, event.path_kmer, event.descaled_event_mean, event.posterior_probability);
//    unlock position
      lk.unlock();
    }
  }
//  write data to binary file
  void write_to_file(path& output_file) {
    data.write_to_file(output_file);
  }
 private:
//  mutexes
  std::vector<mutex> locks;
};

#endif //EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_
