//
// Created by Andrew Bailey on 2019-06-17.
//

#ifndef EMBED_FAST5_TOP_KMERS_H
#define EMBED_FAST5_TOP_KMERS_H

#include <stdlib.h>
#include <string>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/heap/priority_queue.hpp>
#include <boost/coroutine2/all.hpp>
#include <boost/filesystem.hpp>
#include "omp.h"

using namespace std;
using namespace boost::heap;

struct eventkmer
{
  string kmer;
  float mean;
  string strand;
  float prob;
  eventkmer(string kmer, float mean, string strand, float prob) :
      kmer(move(kmer)), mean(mean), strand(move(strand)), prob(prob)
  {
  }
  ~eventkmer()
  = default;
};


typedef boost::coroutines2::coroutine<eventkmer> event_kmer_coro;


class AssignmentFile
{
public:
  explicit AssignmentFile(const string& input_reads_filename);
  ~AssignmentFile();

  void assignment_coroutine(event_kmer_coro::push_type& yield);
  event_kmer_coro::pull_type iterate();
  int64_t get_k();

  string file_path;
  int64_t k;

};


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
  void write_to_file(boost::filesystem::path& output_path);


  string alphabet;
  int alphabet_size;
  int kmer_length;
  int n_kmers;
  size_t max_heap;
  std::vector<boost::heap::priority_queue<eventkmer>> kmer_queues;
  std::vector<omp_lock_t> locks;

};

void generate_master_assignment_table(string assignment_dir, string output_dir, int heap_size, string& alphabet);
int top_kmers_main(int argc, char** argv);


#endif //EMBED_FAST5_TOP_KMERS_H
