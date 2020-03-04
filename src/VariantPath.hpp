//
// Created by Andrew Bailey on 2019-08-21.
//

#ifndef EMBED_FAST5_SRC_VARIANTPATH_HPP_
#define EMBED_FAST5_SRC_VARIANTPATH_HPP_

// embed lib
#include "PositionsFile.hpp"
#include "VariantCall.hpp"
// std lib
#include <vector>
#include <string>


using namespace std;

class VariantPath {
 public:
  map<string, uint64_t> num_positions;
  map<string, uint64_t> num_ids;
  map<string, vector<vector<string>>> index_to_variant;
  map<string, map<uint64_t, uint64_t>> position_to_path_index;
  map<string, vector<uint64_t>> path_multiplier;
  string all_variant_chars;
  PositionsFile pf;

  explicit VariantPath(const string &positions_file_path);
  VariantPath();
  ~VariantPath();

  void load_positions_file(const string &path);
  uint64_t path_to_id(const string &contig_strand, vector<uint64_t> path);
  vector<uint64_t> id_to_path(const string &contig_strand, uint64_t path_id);
  vector<uint64_t> variant_call_to_path(const string &contig_strand, vector<VariantCall> &variants);
  uint64_t variant_call_to_id(const string &contig_strand, vector<VariantCall> &variants);
  string path_to_bases(const string &contig_strand, vector<uint64_t> &path);

//  vector<uint64_t> all_ids_through_node(const string &contig_strand, uint64_t position, uint64_t variant_index);
};

#endif //EMBED_FAST5_SRC_VARIANTPATH_HPP_
