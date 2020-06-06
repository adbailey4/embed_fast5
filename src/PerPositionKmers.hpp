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
  PerPositionKmers(ReferenceHandler &reference,
                   set<char> alphabet = {'A', 'C', 'G', 'T'},
                   uint64_t kmer_length = 6,
                   uint64_t num_locks = 1000,
                   bool two_d = false) :
      alphabet(move(alphabet)), kmer_length(move(kmer_length)), num_locks(move(num_locks)), initialized_locks(true)
  {
    this->initialize_locks();
    data.initialize(reference, this->alphabet, this->kmer_length, two_d, false);
  }

  ~PerPositionKmers() = default;
//  data
  uint64_t num_locks;
  uint64_t kmer_length;
  set<char> alphabet;
  bool initialized_locks;
  EventDataHandler data;
// methods
// initalize locks
  void initialize_locks() {
    locks = std::vector<std::mutex>(this->num_locks);
  }

  //  read in alignment file data
  void process_alignment(AlignmentFile &af) {
    uint64_t hash_value;
//    set<char> alphabet;
    for (auto &event: af.iterate()){
      //    lock position
      hash_value = compute_string_hash(event.path_kmer);
//      alphabet = add_string_to_set(alphabet, event.path_kmer);
//      cout << event.path_kmer << ' ' << hash_value << endl;
      std::unique_lock<std::mutex> lk(this->locks[hash_value % this->num_locks]);
      data.add_kmer_event(event.contig, af.strand, event.strand, event.reference_index, event.path_kmer, event.descaled_event_mean, event.posterior_probability);
//    unlock position
      lk.unlock();
    }
//    std::unique_lock<std::mutex> lk(alphabet_lock);
//    data.add_set_to_alphabet(alphabet);
//    lk.unlock();
  }
//  write data to binary file
  void write_to_file(path& output_file) {
    data.write_to_file(output_file);
  }
 private:
//  mutexes
  std::vector<mutex> locks;
//  mutex alphabet_lock;
};

#endif //EMBED_FAST5_SRC_PERPOSITIONKMERS_HPP_
