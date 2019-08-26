//
// Created by Andrew Bailey on 2019-08-21.
//

// embed library
#include "VariantPath.hpp"

// std lib
#include <iostream>

using namespace std;

/**
 * Constructor for Variant Path
 * @param positions_file_path: path to positions file
 */
VariantPath::VariantPath(const string &positions_file_path) {
  this->load_positions_file(positions_file_path);
}

/**
 * Default constructor for Variant Path
 */
VariantPath::VariantPath()
= default;

/**
 * Default de-constructor for Variant Path
 */
VariantPath::~VariantPath()
= default;

/**
 * Load variants from positions file
 * @param positions_file_path: path to positions file
 */
void VariantPath::load_positions_file(const string &positions_file_path) {
  pf.load_positions_map(positions_file_path);
  uint64_t variant_count;
  std::set<string> variants;
  for (std::pair<std::string, map<uint64_t, tuple<string, string>>> contig_map : pf.positions_map) {
    uint64_t multiplier = 1;
    string change_to;
    uint64_t start;
    for (std::pair<uint64_t, tuple<string, string>> position_map: contig_map.second){
//      contig_map.first is the contig_strand
      path_multiplier[contig_map.first].push_back(multiplier);
      position_to_path_index[contig_map.first][position_map.first] = num_positions[contig_map.first];
      num_positions[contig_map.first] += 1;
//      get characters to change into
      change_to = get<1>(position_map.second);
      variant_count = change_to.length() + 1;
      multiplier *= variant_count;
      sort(change_to.begin(), change_to.end());
//      keep track of variant ambiguous bases
      variants.insert(pf.ambig_bases[change_to]);
  //      keep track index to base
      index_to_variant[contig_map.first] = {{0, "_"}};
      start = 1;
      for (auto &i: change_to){
        index_to_variant[contig_map.first][start] = i;
      }
    }
    num_ids[contig_map.first] = multiplier;
  }
  string variants_string;
  for (const string &v: variants){
    variants_string += v;
  }
  sort(variants_string.begin(), variants_string.end());
  all_variant_chars = variants_string;
}

/**
 * Convert encoded path into path ID
 * @param contig_strand: contig+strand string
 * @param path: vector of indices representing path through variant positions
 * @return: path id
 */
uint64_t VariantPath::path_to_id(const string &contig_strand, vector<uint64_t> path){
  uint64_t id = 0;
  uint64_t number_of_variants = num_positions[contig_strand];
  for (uint64_t i = 0; i < number_of_variants; ++i) {
    id += path[i] * path_multiplier[contig_strand][i];
  }
  return id;

}

/**
 * Convert path id into encoded path through variants
 @param contig_strand: contig+strand string
 @param: path_id: id of path through variants
 @return: vector of indices representing path through variant positions
 */
vector<uint64_t> VariantPath::id_to_path(const string &contig_strand, uint64_t path_id){
  uint64_t number_of_variants = num_positions.at(contig_strand);
  vector<uint64_t> path;
  path.resize(number_of_variants);
  uint64_t temp;

  for (int64_t i = number_of_variants; i > 0; --i) {
    temp = path_id / path_multiplier[contig_strand][i-1];
    path[i-1] = temp;
    path_id = path_id - (path_multiplier[contig_strand][i-1] * temp);
  }
  return path;

}

/**
 * Convert vector of variants to path
 @param variants: vector of variant calls
 */
vector<uint64_t> VariantPath::variant_call_to_path(vector<VariantCall> &variants){
  vector<uint64_t> path;
  string contig_strand = variants[0].contig+variants[0].strand;
  path.resize(num_positions[contig_strand], 0);
  for (auto &i: variants){
    uint64_t max_index = std::max_element(i.normalized_probs.begin(), i.normalized_probs.end())-i.normalized_probs.begin();
    path[position_to_path_index[contig_strand][i.reference_index]] = max_index + 1;
  }
  return path;
}


/**
 * Convert vector of variants to path
 @param variants: vector of variant calls
 */
uint64_t VariantPath::variant_call_to_id(vector<VariantCall> &variants){
  return this->path_to_id(variants[0].contig+variants[0].strand, this->variant_call_to_path(variants));
}

/**
 * Convert vector of indexes to characters
 * @param contig_strand: contig+strand string
 * @param path: vector of indices representing path through variant positions
 */
vector<string> VariantPath::path_to_bases(const string &contig_strand, vector<uint64_t> &path){
  vector<string> bases;
  bases.reserve(path.size());
  for (auto &index: path){
    bases.push_back(this->index_to_variant[contig_strand][index]);
  }
  return bases;
}



///**
// * Convert path id into encoded path through variants
// * @param contig_strand: contig+strand string
// * @param: path id
// * @return: vector of indices representing path through variant positions
// */
//vector<uint64_t> VariantPath::all_ids_through_node(const string &contig_strand, uint64_t position, uint64_t variant_index){
//  uint64_t n_variants = num_variants[contig_strand];
//  uint64_t step_size = path_multiplier[contig_strand][position];
//  uint64_t total_steps = step_size* ;
//
//  vector<uint64_t> path_ids;
//  for (uint64_t i = 0; i < n_variants; ++i){
//    if (i % )
//  }
//}
