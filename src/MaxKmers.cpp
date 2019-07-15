//
// Created by Andrew Bailey on 2019-07-15.
//

#include "MaxKmers.hpp"
#include "EmbedUtils.hpp"

using namespace embed_utils;

/**
Destroy locks at destruction of class
*/
MaxKmers::~MaxKmers() {
  this->destroy_locks();
}

/**
Initialize locks and heaps and other important data structures
*/
MaxKmers::MaxKmers(size_t heap_size, string alphabet, int kmer_length) :
    alphabet(sort_string(alphabet)), alphabet_size(alphabet.length()),
    kmer_length(kmer_length), n_kmers(pow(this->alphabet_size, this->kmer_length)), max_heap(heap_size)
{
  this->initialize_heap();
  this->initialize_locks();
}

/**
Create all kmers given an alphabet and kmer length

eg. create_kmers(BA, 1) = {"A", "B"}

@param alphabet: string with each character is apart of the alphabet
@param kmer_length: length of the string
@return vector of alphabetically ordered strings
*/

vector<string> MaxKmers::create_kmers(string& alphabet1, int kmer_length1) {
  return all_string_permutations(alphabet1, kmer_length1);
}

/**
Get kmer index based on alphabet and kmer length

eg. get_kmer_index(AAAAA) = 0

@param kmer: input kmer
@return index
*/

size_t MaxKmers::get_kmer_index(string& kmer){
  assert(this->alphabet_size > 0);
  int64_t id = 0;
  int64_t step = 1;
  int64_t kmer_len = kmer.length();
  int64_t index;

  for (int64_t i = kmer_len - 1; i >= 0; i--) {
    index = this->alphabet.find(kmer[i]);
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

string MaxKmers::get_index_kmer(size_t kmer_index) {
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
Override < operator in order to make the priority queue a min heap

@param a: first eventkmer
@param b: second eventkmer

@return a.prob > b.prob
*/

bool operator<(const eventkmer& a, const eventkmer& b) {
  return a.prob > b.prob;
}

/**
Initialize vector of heaps
*/
void MaxKmers::initialize_heap() {
  for (int i=0; i < this->n_kmers; i++){
    boost::heap::priority_queue<eventkmer> kmer_queue;
    this->kmer_queues.push_back(kmer_queue);
  }
}

/**
Initialize vector of locks
*/
void MaxKmers::initialize_locks() {
  for (int i=0; i<this->n_kmers; i++){
    omp_lock_t new_lock;
    omp_init_lock(&(new_lock));
    locks.push_back(new_lock);
  }
}

/**
Destroy vector of locks
*/
void MaxKmers::destroy_locks() {
  for (int i=0; i<this->n_kmers; i++){
    omp_destroy_lock(&(locks[i]));
  }
}

/**
Add kmer data to the heap data structure for said kmer if probability is greater than the smallest probability

@param kmer: kmer to add
@param b: second eventkmer

@return a.prob > b.prob
*/
void MaxKmers::add_to_heap(eventkmer kmer_struct){
  size_t index = this->get_kmer_index(kmer_struct.kmer);
  omp_set_lock(&(this->locks[index]));
  if (this->kmer_queues[index].empty()) {
//    add to queue if not at capacity
    this->kmer_queues[index].push(kmer_struct);
//    if at capacity check to see if prob is greater than min
  } else if (this->kmer_queues[index].top().prob < kmer_struct.prob || this->kmer_queues[index].size() < this->max_heap)
  {
    this->kmer_queues[index].push(kmer_struct);
    while (this->kmer_queues[index].size() > this->max_heap){
      this->kmer_queues[index].pop();
    }
  }
  omp_unset_lock(&(this->locks[index]));

}
/**
 * Write all kmers in the kmer_queues to output path
 * @param output_path
 */
void MaxKmers::write_to_file(boost::filesystem::path& output_path){
  std::ofstream out_file;
  out_file.open(output_path.string());
  for (auto &pq: this->kmer_queues){
    for (auto &event: pq){
      out_file << event.kmer << '\t' << event.strand << '\t' << event.mean << '\t' <<  event.prob << '\n';
    }
  }
  out_file.close();
}

/**
 * Write all kmers in the kmer_queues to output path and write info about the kmers to log_path
 * @param output_path
 * @param log_path
 */
void MaxKmers::write_to_file(boost::filesystem::path& output_path, boost::filesystem::path& log_path){
  std::ofstream out_file;
  out_file.open(output_path.string());

  std::ofstream out_log;
  out_log.open(log_path.string());
  out_log << "kmers" << '\t' << "num_events" << '\t' << "min_prob" << '\n';

  size_t kmer_index = 0;
  for (auto &pq: this->kmer_queues){
// loop through all queues
    for (auto &event: pq){
      out_file << event.kmer << '\t' << event.strand << '\t' << event.mean << '\t' <<  event.prob << '\n';
//      keep track of number of events
    }
    //    log info about queue
    int n_events = pq.size();
    float min_p;
    if (n_events > 0){
      eventkmer top = pq.top();
      min_p = top.prob;
    } else{
      min_p = 0.0;
    }
    out_log << this->get_index_kmer(kmer_index) << '\t' << n_events << '\t' << min_p << '\n';
    kmer_index += 1;

  }
  out_file.close();
  out_log.close();
}