//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_MAXKMERS_HPP_
#define EMBED_FAST5_SRC_MAXKMERS_HPP_

#include "AssignmentFile.hpp"
#include <boost/filesystem.hpp>
#include <boost/heap/priority_queue.hpp>
#include "omp.h"

using namespace std;
using namespace boost::filesystem;

class MaxKmers{
 public:
  explicit MaxKmers(size_t heap_size, string alphabet, int kmer_length);
  ~MaxKmers();

  vector<string> create_kmers(string& alphabet1, int kmer_length1);
  size_t get_kmer_index(string& kmer);
  string get_index_kmer(size_t kmer_index);
  void initialize_heap();
  void initialize_locks();
  void destroy_locks();
  void add_to_heap(eventkmer kmer_struct);
  void write_to_file(path& output_path);
  void write_to_file(path& output_path, path& log_path);

  string alphabet;
  int alphabet_size;
  int kmer_length;
  int n_kmers;
  size_t max_heap;
  std::vector<boost::heap::priority_queue<eventkmer>> kmer_queues;
  std::vector<omp_lock_t> locks;
};

#endif //EMBED_FAST5_SRC_MAXKMERS_HPP_
