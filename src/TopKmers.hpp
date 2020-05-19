//
// Created by Andrew Bailey on 2019-06-17.
//

#ifndef EMBED_FAST5_TOP_KMERS_H
#define EMBED_FAST5_TOP_KMERS_H

#include "AssignmentFile.hpp"
#include "AlignmentFile.hpp"
#include "MaxKmers.hpp"
#include <unordered_set>
#include <thread>
#include <atomic>

using namespace std;
static std::exception_ptr globalExceptionPtr = nullptr;

void generate_master_kmer_table_wrapper(vector<string> event_table_files,
                                        string &output_file,
                                        string &log_file,
                                        uint64_t heap_size,
                                        string &alphabet,
                                        double min_prob,
                                        uint64_t n_threads,
                                        bool verbose,
                                        bool write_full);


/**
 * Worker for generate_master_kmer_table. Add kmers to heap from alignment file
 *
 * @tparam T1: Event table file parsing class
 * @tparam T2: Event table data type
 * @param signalalign_output_files: reference to vector of signalalign files
 * @param max_kmers: templated reference to thread safe queue
 * @param job_index: atomic index for selecting output files to process
 * @param n_files: max number of files to process
 * @param verbose: option for printing files processed
 */
template<class T1, class T2>
void bin_max_kmer_worker(vector<path>& signalalign_output_files, T2& max_kmers, atomic<uint64_t>& job_index, uint64_t n_files, bool& verbose) {
  try {
    while (job_index < n_files and !globalExceptionPtr) {
      // Fetch add
      uint64_t thread_job_index = job_index.fetch_add(1);
      path current_file = signalalign_output_files[thread_job_index];
      if (verbose) {
//      cout << current_file << "\n";
        // Print status update to stdout
        cerr << "\33[2K\rParsed: " << current_file << flush;

      }
      T1 af(current_file.string());
      for (auto &event: af.iterate()) {
        max_kmers.add_to_heap(event);
      }
    }
  } catch(...){
    globalExceptionPtr = std::current_exception();
  }
}

/**
 * Generate the master table by parsing event table files and outputting the top n kmers to a files
 *
 * @tparam T1: Event table file parsing class
 * @tparam T2: Event table data type
 * @param output_file: path to output file
 * @param log_file: path to output log file
 * @param alphabet: alphabet used to generate kmers
 * @param heap_size: number of max kmers to keep for each kmer
 * @param min_prob: minimum probability
 * @param n_threads: set number of threads to use: default 2
 * @param verbose: boolean verbose option
 */
template<class T1, class T2>
void generate_master_kmer_table(vector<string> &sa_output_paths,
                                string &output_file,
                                string &log_file,
                                string &alphabet,
                                int heap_size,
                                double min_prob = 0.0,
                                unsigned int n_threads = 1,
                                bool verbose = false,
                                bool write_full = false) {

//  filter out empty files and check if there are any left
  vector<path> all_tsvs = filter_emtpy_files(sa_output_paths, ".tsv");
  throw_assert(!all_tsvs.empty(), "There are no valid .tsv files")
  int64_t number_of_files = all_tsvs.size();
//  get kmer length
  T1 af(all_tsvs[0].string());
  int64_t kmer_length = af.get_k();
//  initialize heap, job index, threads, and count number of files
  MaxKmers<T2> mk(heap_size, alphabet, kmer_length, min_prob);
  atomic<uint64_t> job_index(0);
  vector<thread> threads;
  globalExceptionPtr = nullptr;
  // Launch threads
  for (uint64_t i=0; i<n_threads; i++){
      threads.emplace_back(thread(bin_max_kmer_worker<T1, MaxKmers<T2>>,
                                  ref(all_tsvs),
                                  ref(mk),
                                  ref(job_index),
                                  ref(number_of_files),
                                  ref(verbose)));
  }
  // Wait for threads to finish
  for (auto& t: threads){
    t.join();
  }
  if (globalExceptionPtr){
    std::rethrow_exception(globalExceptionPtr);
  }
  if (verbose){
    cerr << "\n" << flush;
  }
  path output_path(output_file);
  path log_path(log_file);

  mk.write_to_file(output_path, log_path, write_full);
}

auto top_kmers_main(int argc, char** argv) -> int;


#endif //EMBED_FAST5_TOP_KMERS_H
