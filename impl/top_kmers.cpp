//
// Created by Andrew Bailey on 2019-06-17.
//

#include <math.h>
#include <getopt.h>
#include <sstream>
#include <boost/coroutine2/all.hpp>
#include <boost/filesystem.hpp>
#include "omp.h"
#include "top_kmers.hpp"
#include "embed_utils.hpp"

using namespace std;
using namespace embed_utils;


AssignmentFile::~AssignmentFile()
= default;

AssignmentFile::AssignmentFile(const string& input_reads_filename){
    this->file_path = input_reads_filename;
}

/**
Create a push type coroutine for parsing an assignment file

 file format:
GCCTTA	t	83.709275	1.000000

*/
void AssignmentFile::assignment_coroutine(event_kmer_coro::push_type& yield){
    std::ifstream in_file(this->file_path.c_str());
    if(in_file.good()) {
        // read the file
        std::string line;
        while(getline(in_file, line)) {
            std::vector<std::string> fields = embed_utils::split(line, '\t');
            string kmer = fields[0];
            string strand = fields[1];
            float mean = convert_to_float(fields[2]);
            float prob = convert_to_float(fields[3]);
            yield(eventkmer(kmer, mean, strand, prob));
        }
    } else  {
        cout << "Error loading file: " << this->file_path << "\n";
    }
    in_file.close();
}

/**
Call push type coroutine to create a pull type coroutine
*/
event_kmer_coro::pull_type AssignmentFile::iterate(){
  event_kmer_coro::pull_type seq {bind(&AssignmentFile::assignment_coroutine, this, std::placeholders::_1)};
  return seq;
}

/**
Get the kmer size for a given file
*/
int64_t AssignmentFile::get_k(){
  std::ifstream in_file(this->file_path.c_str());
  if (in_file.good()) {
    // read the file
    std::string line;
    getline(in_file, line);
    std::vector<std::string> fields = split(line, '\t');
    this->k = fields[0].length();
  }
  in_file.close();
  return this->k;
}

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

  for (auto &pq: this->kmer_queues){
    int counter = 0;
    float min_p = 2.0;
    string kmer;
    for (auto &event: pq){
      out_file << event.kmer << '\t' << event.strand << '\t' << event.mean << '\t' <<  event.prob << '\n';
      counter += 1;
      if (min_p < event.prob) {
        min_p = event.prob;
        kmer = event.kmer;
      }
    }
    out_log << kmer << '\t' << counter << '\t' << min_p << '\n';
  }
  out_file.close();
  out_log.close();
}

/**
 * Generate the master assignment table by parsing assignment files and outputting the top n kmers to a files
 *
 * @param assignment_dir: path to assignment files directory
 * @param output_dir: path to output directory where new builtAssignment.tsv will be written
 * @param heap_size: number of max kmers to keep for each kmer
 * @param alphabet: alphabet used to generate kmers
 */

void generate_master_assignment_table(string assignment_dir, string output_dir, int heap_size, string& alphabet){

  path p(assignment_dir);
  path output_path = make_dir(output_dir);

  directory_iterator end_itr;
//    Get all tsvs to process
  vector<path> all_tsvs;
  int kmer_length = 0;
  int counter = 0;
  for (directory_iterator itr(p); itr != end_itr; ++itr) {
//        filter for files that are regular, end with tsv and are not empty
    if (is_regular_file(itr->path()) and itr->path().extension().string() == ".tsv" and getFilesize(itr->path().string()) > 0) {
      all_tsvs.push_back(itr->path());
      if (counter == 0){
        AssignmentFile af = AssignmentFile(itr->path().string());
        kmer_length = af.get_k();
        counter += 1;
      }
    }
  }

  MaxKmers mk = MaxKmers(heap_size, alphabet, kmer_length);

  int64_t number_of_files = all_tsvs.size();
  path* array_of_files = &all_tsvs[0];

// looping through the files
#pragma omp parallel for shared(array_of_files, mk, number_of_files)
  for(int64_t i=0; i < number_of_files; i++) {

    path current_file = array_of_files[i];
    cout << current_file << "\n";
    AssignmentFile af = AssignmentFile(current_file.string());
    for (auto &event: af.iterate()){
      mk.add_to_heap(event);
    }
//        if (current_file.filename().string() == "0a4e473d-4713-4c7f-9e18-c465ea6d5b8c.sm.forward.tsv"){
//        }
  }
  path output_file = output_path / "builtAssignment.tsv";
  path log_file = output_path / "built_log.tsv";

  mk.write_to_file(output_file, log_file);
}




// Getopt
//
#define SUBPROGRAM "top_kmers"
#define TOP_KMER_VERSION "0.0.1"
#define THIS_NAME "top_kmers"
#define PACKAGE_BUGREPORT2 "None"

static const char *TOP_KMER_VERSION_MESSAGE =
    SUBPROGRAM " Version " TOP_KMER_VERSION "\n";

static const char *TOP_KMER_USAGE_MESSAGE =
    "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --assignment_dir=DIR --output_dir=DIR --threads=NUMBER --heap_size=NUMBER --alphabet=STRING\n"
    "Filters all assignment files into a master assignment table with the top n kmers.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -d, --assignment_dir=DIR             directory of signalalign assignment files\n"
    "  -o, --output_dir=DIR                 path to directory to output master assignment table\n"
    "  -t, --threads=NUMBER                 number of threads\n"
    "  -s, --heap_size=NUMBER               size of heap for each kmer\n"
    "  -a, --alphabet=STRING                alphabet for kmers\n"


    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static std::string assignment_dir;
static std::string output_dir;
static unsigned int threads;
static unsigned int heap_size;
static string alphabet;
static int num_threads = 1;
}

static const char* shortopts = "a:d:s:t:o:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "assignment_dir",   required_argument, nullptr, 'd' },
    { "output_dir",       required_argument, nullptr, 'o' },
    { "heap_size",        required_argument, nullptr, 's' },
    { "alphabet",         required_argument, nullptr, 'a' },
    { "threads",          optional_argument, nullptr, 't' },
    { "help",             no_argument,       nullptr, OPT_HELP },
    { "version",          no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

void parse_top_kmers_main_options(int argc, char** argv)
{
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
    std::istringstream arg(optarg != nullptr ? optarg : "");
    switch (c) {
      case 'd': arg >> opt::assignment_dir; break;
      case 'o': arg >> opt::output_dir; break;
      case 't': arg >> opt::threads; break;
      case 'a': arg >> opt::alphabet; break;
      case 's': arg >> opt::heap_size; break;
      case 'v': opt::verbose++; break;
      case OPT_HELP:
        std::cout << TOP_KMER_USAGE_MESSAGE;
        exit(EXIT_SUCCESS);
      case OPT_VERSION:
        std::cout << TOP_KMER_VERSION_MESSAGE;
        exit(EXIT_SUCCESS);
    }
  }

  if (argc - (optind+1) > 0) {
    cout << argc << " " << optind << "\n";
    std::cerr << SUBPROGRAM ": too many arguments\n";
    die = true;
  }

  if(opt::assignment_dir.empty()) {
    std::cerr << SUBPROGRAM ": a --assignment_dir directory must be provided\n";
    die = true;
  }
  if(opt::alphabet.empty()) {
    std::cerr << SUBPROGRAM ": a --alphabet must be provided\n";
    die = true;
  }
  if(opt::heap_size <= 0) {
    std::cerr << SUBPROGRAM ": a --heap_size must be provided\n";
    die = true;
  }
  if(opt::output_dir.empty()) {
    std::cerr << SUBPROGRAM ": a --output_dir file must be provided\n";
    die = true;
  }
  if (die)
  {
    std::cout << "\n" << TOP_KMER_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}


int top_kmers_main(int argc, char** argv)
{
  parse_top_kmers_main_options(argc, argv);

#ifndef H5_HAVE_THREADSAFE
  if(opt::num_threads > 1) {
    fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
    fprintf(stderr, "Please recompile built-in libhdf5 or run with -t 1\n");
    exit(1);
  }
#endif

  omp_set_num_threads(opt::threads); // Use 4 threads for all consecutive parallel regions
  generate_master_assignment_table(opt::assignment_dir, opt::output_dir, opt::heap_size, opt::alphabet);

  return EXIT_SUCCESS;
}




