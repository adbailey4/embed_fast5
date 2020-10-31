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
PerPositonKmers::PerPositonKmers(uint64_t num_locks, path output_dir, const ReferenceHandler& reference) :
    num_locks(num_locks), output_dir(std::move(output_dir))
{
  this->initialize_map(reference);

}

PerPositonKmers::PerPositonKmers(uint64_t num_locks, path output_dir) :
    num_locks(num_locks), output_dir(std::move(output_dir))
{
}

void PerPositonKmers::initialize_map(ReferenceHandler reference) {
  vector<string> contigs = reference.get_chromosome_names();
  for (auto &i: contigs){
    uint64_t contig_length = reference.get_chromosome_sequence_length(i);

  }
}


path PerPositonKmers::create_file_path(const string& contig, uint64_t position, const string& kmer){
  path file_name = contig+"_"+to_string(position)+"_"+kmer+".tsv";
  path file_path = this->output_dir / file_name;
  return file_path;
}


//    if (this->initialized_locks){
//      omp_set_lock(&(this->locks[i.reference_index % this->num_locks]));
//    }
