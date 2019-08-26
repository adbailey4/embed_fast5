//
// Created by Andrew Bailey on 2019-08-22.
//

// embed libs
#include "LoadVariantPaths.hpp"
#include "AlignmentFile.hpp"
#include "EmbedUtils.hpp"

using namespace embed_utils;

/**
 * Constructor for LoadVariantPaths
 * @param positions_file_path: path to positions file
 * @param full_sa_dir: path to directory with "full" sa output tsv files

 */
LoadVariantPaths::LoadVariantPaths(const string &positions_file_path, const string &full_sa_dir) {
  this->vp.load_positions_file(positions_file_path);
  this->read_in_sa_files(full_sa_dir);
}

/**
 * Default de-constructor for LoadVariantPaths
 */
LoadVariantPaths::~LoadVariantPaths()
= default;

/** Read in "full" signal align files, generate variant calls and
 * @param full_sa_dir: path to directory with "full" sa output tsv files
 */
void LoadVariantPaths::read_in_sa_files(const string &full_sa_dir){
  path sa_dir = full_sa_dir;
  string contig_strand;
  uint64_t path_id;
  string ext = ".tsv";
  for (path &i: list_files_in_dir(sa_dir, ext)){
//    load file
    AlignmentFile af(i.string());
//    get variants
    vector<VariantCall> calls = af.get_variant_calls(vp.all_variant_chars, &vp.pf.ambig_bases);
//    get path corresponding to variants
    path_id = vp.variant_call_to_id(calls);
//  keep track fo variant path and counts per path
    contig_strand = calls[0].contig + calls[0].strand;
    read_id_path_id_map[contig_strand][af.read_id] = path_id;
    counts[contig_strand][path_id] += 1;
  }
}

void LoadVariantPaths::write_per_read_calls(const string &output_path){
  std::ofstream out_file;
  out_file.open(output_path);
  // loop through all contigs
  out_file << "contig" << '\t' << "read_id" << '\t' << "path_id" << '\t' << "path" << '\n';

  for (auto &contig_map: read_id_path_id_map){
//    loop through contig map
    for (auto &read_id: contig_map.second){
      out_file << contig_map.first << '\t' << read_id.first << '\t' << read_id.second <<  '\t';
      vector<uint64_t> path = vp.id_to_path(contig_map.first, read_id.second);
      for (auto &node: path){
        out_file << node;
      }
      out_file << '\t';
      for (auto &base: vp.path_to_bases(contig_map.first, path)){
        out_file << base;
      }
      out_file << "\n";
    }
  }
  out_file.close();
}

void LoadVariantPaths::write_per_path_counts(const string &output_path){
  std::ofstream out_file;
  out_file.open(output_path);
  // loop through all contigs
  out_file << "contig" << '\t' << "path_id" << '\t' << "path" << '\t' << "counts" << '\n';

  for (auto &contig_map: counts){
//    loop through contig map
    for (auto &id_map: contig_map.second){
      out_file << contig_map.first << '\t' << id_map.first << '\t' ;
      vector<uint64_t> path = vp.id_to_path(contig_map.first, id_map.first);
      for (auto &node: path){
        out_file << node;
      }
      out_file << '\t';

      for (auto &base: vp.path_to_bases(contig_map.first, path)){
        out_file << base;
      }
      out_file << '\t' << id_map.second << '\n';
    }
  }
  out_file.close();
}