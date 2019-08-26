//
// Created by Andrew Bailey on 2019-08-22.
//

#ifndef EMBED_FAST5_SRC_LOADVARIANTPATHS_HPP_
#define EMBED_FAST5_SRC_LOADVARIANTPATHS_HPP_

// embed libs
#include "VariantPath.hpp"

using namespace std;

class LoadVariantPaths {
 public:
  VariantPath vp;
  map<string, map<string, uint64_t>> read_id_path_id_map;
  map<string, map<uint64_t, uint64_t>> counts;

  explicit LoadVariantPaths(const string &positions_file_path, const string &full_sa_dir);
  ~LoadVariantPaths();
  void read_in_sa_files(const string &full_sa_dir);
  void write_per_read_calls(const string &output_path);
  void write_per_path_counts(const string &output_path);
};

#endif //EMBED_FAST5_SRC_LOADVARIANTPATHS_HPP_
