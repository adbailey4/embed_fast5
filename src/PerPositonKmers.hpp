//
// Created by Andrew Bailey on 2019-07-24.
//

#ifndef EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_
#define EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_

#include "ReferenceHandler.hpp"
#include <boost/filesystem.hpp>
#include <omp.h>
#include <map>
#include <fstream>

using namespace std;
using namespace boost::filesystem;

class PerPositonKmers {
  PerPositonKmers(uint64_t num_locks, path output_dir, const ReferenceHandler& reference);
  PerPositonKmers(uint64_t num_locks, path output_dir);
  ~PerPositonKmers();
  uint64_t num_locks;
  path output_dir;
  bool initialized_locks = false;
  void initialize_map(ReferenceHandler reference);
  path create_file_path(const string& contig, uint64_t position, const string& kmer);
 private:
  map<string, std::fstream> file_handles;
  std::vector<omp_lock_t> locks;

};

#endif //EMBED_FAST5_SRC_PERPOSITONKMERS_HPP_
