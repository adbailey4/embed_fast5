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
LoadVariantPaths::LoadVariantPaths(const string &positions_file_path, const string &full_sa_dir, bool rna, uint64_t num_locks) {
  this->vp.load_positions_file(positions_file_path);
  this->num_locks = num_locks;
  this->initialize_locks();
  this->read_in_sa_files(full_sa_dir, rna);
}

/**
 * Default de-constructor for LoadVariantPaths
 */
LoadVariantPaths::~LoadVariantPaths()
= default;

/**
 * Default de-constructor for LoadVariantPaths
 */
void LoadVariantPaths::initialize_locks() {
  for (uint64_t i=0; i<this->num_locks; i++){
    omp_lock_t new_lock;
    omp_init_lock(&(new_lock));
    locks.push_back(new_lock);
  }
}


/** Read in "full" signal align files, generate variant calls and
 * @param full_sa_dir: path to directory with "full" sa output tsv files
 */
void LoadVariantPaths::read_in_sa_files(const string &full_sa_dir, bool rna){
  path sa_dir = full_sa_dir;

  string ext = ".tsv";
  vector<string> all_tsvs;
  for (path &i: list_files_in_dir(sa_dir, ext)) {
    all_tsvs.emplace_back(i.string());
  }
  int64_t number_of_files = all_tsvs.size();
  string* array_of_files = &all_tsvs[0];
  read_id_path_id_map.resize(number_of_files);
// looping through the files
#pragma omp parallel for shared(array_of_files, number_of_files, rna) default(none)
  for(int64_t i=0; i < number_of_files; i++) {
    //    load file
    AlignmentFile af(array_of_files[i], rna);
//    get variants
    vector<VariantCall> calls = af.get_variant_calls(this->vp.all_variant_chars, &this->vp.pf.ambig_bases);
//    get path corresponding to variants
    string contig_strand = af.contig+af.strand;
    uint64_t path_id = this->vp.variant_call_to_id(contig_strand, calls);
//  keep track fo variant path and counts per path
    read_id_path_id_map[i] = make_tuple(contig_strand, af.read_id, path_id);
//  set lock for creating and addition of path id if needed
    omp_set_lock(&(locks[path_id % num_locks]));
    counts[contig_strand][path_id] += 1;
    omp_unset_lock(&(locks[path_id % num_locks]));
  }
}

/** Write read_id, path_id, path, and nucleotide path for each read
 * @param output_path: path to output file
 */
void LoadVariantPaths::write_per_read_calls(const string &output_path){
  std::ofstream out_file;
  out_file.open(output_path);
  // loop through all contigs
  out_file << "contig_strand" << '\t' << "read_id" << '\t' << "path_id" << '\t' << "path" << '\t' << "nucleotide_path" << '\n';

  for (auto &id_tuple: read_id_path_id_map){
//    loop through contig map
    out_file << get<0>(id_tuple) << '\t' << get<1>(id_tuple) << '\t' << get<2>(id_tuple) <<  '\t';
    vector<uint64_t> path = vp.id_to_path(get<0>(id_tuple), get<2>(id_tuple));
    for (auto &node: path){
      out_file << node;
    }
    out_file << '\t';
    for (auto &base: vp.path_to_bases(get<0>(id_tuple), path)){
      out_file << base;
    }
    out_file << "\n";
  }
  out_file.close();
}

/** For each path_id/ path/ nucleotide path write the number of reads falling into that path
 * @param output_path: path to output file
 */
void LoadVariantPaths::write_per_path_counts(const string &output_path){
  std::ofstream out_file;
  out_file.open(output_path);
  // loop through all contigs
  out_file << "contig_strand" << '\t' << "path_id" << '\t' << "path" << '\t' << "nucleotide_path" << '\t' << "counts" << '\n';

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