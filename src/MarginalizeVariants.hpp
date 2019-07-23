//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_MARGINALIZEVARIANTS_HPP_
#define EMBED_FAST5_SRC_MARGINALIZEVARIANTS_HPP_

#include "AlignmentFile.hpp"
#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include <map>
#include <omp.h>



using namespace std;
using namespace boost::filesystem;

struct bed_line{
  uint64_t start;
  uint64_t stop;
  uint64_t coverage;
  vector<uint64_t> hits;
  string bases;
  bed_line(){
    start=-1;
    stop=-1;
    coverage=0;
  }
};


class MarginalizeVariants {

 public:
  MarginalizeVariants();
  explicit MarginalizeVariants(uint64_t num_locks);
  ~MarginalizeVariants();
  void load_variants(vector<VariantCall>* vector_of_calls);
  void write_to_file(path& path_to_bed);
  void initialize_locks();
  void initialize_locks(uint64_t number_of_locks);
//  per_genomic_position[contig][strand] = map of positions
  uint64_t num_locks{};
  bool initialized_locks{};
  std::vector<omp_lock_t> locks;
  map<std::string, pair<map<uint64_t, bed_line>, map<uint64_t, bed_line>>> per_genomic_position;
};

#endif //EMBED_FAST5_SRC_MARGINALIZEVARIANTS_HPP_
