//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_MAXKMERS_HPP_
#define EMBED_FAST5_SRC_MAXKMERS_HPP_

#include "AssignmentFile.hpp"
#include "AlignmentFile.hpp"

#include "EmbedUtils.hpp"
#include <boost/filesystem.hpp>
#include <boost/heap/priority_queue.hpp>
#include <mutex>

using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;


template <class T>
class MaxKmers{
 public:

  /**
  Initialize locks and heaps and other important data structures
  */
  MaxKmers(size_t heap_size, string alphabet, uint64_t kmer_length, double min_prob= 0.0) :
      alphabet(sort_string(alphabet)), alphabet_size(alphabet.length()),
      kmer_length(kmer_length), n_kmers(pow(this->alphabet_size, this->kmer_length)), max_heap(heap_size),
      min_prob(min_prob)
  {
    throw_assert(this->kmer_length > 0, "Kmer length must be greater than 0.")
    this->initialize_heap();
    this->initialize_locks();
  }
  ~MaxKmers() = default;
//    this->destroy_locks();


  string alphabet;
  int alphabet_size;
  uint64_t kmer_length;
  int n_kmers;
  size_t max_heap;
  double min_prob;
  std::vector<boost::heap::priority_queue<T>> kmer_queues;
  std::vector<mutex> locks;

  /**
  Create all kmers given an alphabet and kmer length

  eg. create_kmers(BA, 1) = {"A", "B"}

  @param alphabet: string with each character is apart of the alphabet
  @param kmer_length: length of the string
  @return vector of alphabetically ordered strings
  */
  vector<string> create_kmers(string& alphabet1, uint64_t kmer_length1) {
    return all_string_permutations(alphabet1, kmer_length1);
  }

  /**
  Get kmer index based on alphabet and kmer length

  eg. get_kmer_index(AAAAA) = 0

  @param kmer: input kmer
  @return index
  */
  size_t get_kmer_index(string& kmer){
    int64_t id = 0;
    int64_t step = 1;
    throw_assert(kmer.length() == this->kmer_length, "Kmer length is different than expected: " << kmer.length() << " != " << this->kmer_length)
    uint64_t index;

    for (int64_t i = this->kmer_length - 1; i >= 0; i--) {
      index = this->alphabet.find(kmer[i]);
      throw_assert( index != string::npos, "Kmer (" + kmer + ") has character not in (" + this->alphabet + ") alphabet")
      id += step * index;
      step *= this->alphabet_size;
    }
    return id;
  }

  /**
  Get kmer based on kmer index.

  eg. get_index_kmer(0) = AAAAA

  @param index: kmer index
  @return kmer
  */
  string get_index_kmer(size_t &kmer_index) {
    string kmer;
    size_t id_remainder = kmer_index;

    for (int64_t i = this->kmer_length - 1; i >= 0; i--) {
      kmer += this->alphabet[id_remainder % this->alphabet_size];
      id_remainder /= this->alphabet_size;
    }
    reverse(kmer.begin(), kmer.end());
    return kmer;
  }

  /**
  Initialize vector of heaps
  */
  void initialize_heap() {
    kmer_queues.resize(n_kmers);
  }

  /**
  Initialize vector of locks
  */
  void initialize_locks() {
    locks = std::vector<std::mutex>(this->n_kmers);;
  }

//  /**
//Destroy vector of locks
//*/
//  void destroy_locks() {
//    for (int i=0; i<this->n_kmers; i++){
//      omp_destroy_lock(&(locks[i]));
//    }
//  }

  /**
   * Write all kmers in the kmer_queues to output path
   * @param output_path
   */
  void write_to_file(boost::filesystem::path &output_path, bool write_full) {
    std::ofstream out_file;
    out_file.open(output_path.string());
    for (auto &pq: this->kmer_queues){
      for (auto &event: pq){
        out_file << event.format_line(write_full);
      }
    }
    out_file.close();
  }

  /**
   * Write all kmers in the kmer_queues to output path and write info about the kmers to log_path
   * @param output_path: path to output file
   * @param log_path: path to output log file
   * @param write_full: boolean option to write full output
   */
  void write_to_file(boost::filesystem::path &output_path, boost::filesystem::path &log_path, bool write_full) {
    std::ofstream out_file;
    out_file.open(output_path.string());

    std::ofstream out_log;
    out_log.open(log_path.string());
    out_log << "kmers" << '\t' << "num_events" << '\t' << "min_prob" << '\n';

    size_t kmer_index = 0;
    for (auto &pq: this->kmer_queues){
// loop through all queues
      for (auto &event: pq){
        out_file << event.format_line(write_full);
//      keep track of number of events
      }
      //    log info about queue
      int n_events = pq.size();
      float min_p;
      if (n_events > 0){
        auto top = pq.top();
        min_p = top.posterior_probability;
      } else{
        min_p = 0.0;
      }
      out_log << this->get_index_kmer(kmer_index) << '\t' << n_events << '\t' << min_p << '\n';
      kmer_index += 1;

    }
    out_file.close();
    out_log.close();
  }

  /**
  Add kmer data to the heap data structure for said kmer if probability is greater than the smallest probability

  @param kmer: kmer to add
  @param b: second eventkmer

  @return a.posterior_probability > b.posterior_probability
  */
  void add_to_heap(T& kmer_struct){
    size_t index = this->get_kmer_index(kmer_struct.path_kmer);
    if (kmer_struct.posterior_probability >= min_prob){
      std::unique_lock<std::mutex> lock(this->locks[index]);
      boost::heap::priority_queue<T>& queue = this->kmer_queues[index];
      if (queue.empty()) {
//    add to queue if not at capacity
        queue.push(kmer_struct);
//    if at capacity check to see if prob is greater than min
      } else if (queue.top().posterior_probability < kmer_struct.posterior_probability ||
          queue.size() < this->max_heap) {
        queue.push(kmer_struct);
        while (queue.size() > this->max_heap){
          queue.pop();
        }
      }
      lock.unlock();
    }
  }
};

bool operator<(const eventkmer& a, const eventkmer& b);
bool operator<(const FullSaEvent& a, const FullSaEvent& b);


#endif //EMBED_FAST5_SRC_MAXKMERS_HPP_
