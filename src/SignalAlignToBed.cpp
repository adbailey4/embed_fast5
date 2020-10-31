//
// Created by Andrew Bailey on 2019-07-14.
//

#include "SignalAlignToBed.hpp"
#include "EmbedUtils.hpp"
#include "MarginalizeVariants.hpp"
#include "ConcurrentQueue.hpp"
#include <getopt.h>
#include <iostream>
#include <boost/filesystem.hpp>
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;
static std::exception_ptr globalExceptionPtr = nullptr;


/**
 * worker which parses "full" signalalign file and passes a vector of variant calls to a ConcurrentQueue object
 *
 * @param signalalign_output_files: reference to vector of signalalign files
 * @param mv: MarginalizeVariants class object
 * @param variant_queue: templated reference to thread safe queue
 * @param job_index: atomic index for selecting output files to process
 * @param n_files: max number of files to process
 * @param verbose: option for printing files processed
 * @param rna: boolean option if reads are rna
 * @param ambig_bases: string of all ambiguous characters
 * @param ambig_bases_map: map of ambig bases to expected chars
 */
void get_variants_worker(
    vector<path>& signalalign_output_files,
    MarginalizeVariants& mv,
    ConcurrentQueue<tuple<string, vector<VariantCall>>>& variant_queue,
    atomic<uint64_t>& job_index,
    int64_t& n_files,
    bool& verbose,
    bool& rna,
    string &ambig_bases,
    std::map<string, string>& ambig_bases_map) {
  try {
    tuple<string, vector<VariantCall>> read_id_and_variants;
    while (job_index < n_files and !globalExceptionPtr) {
      // Fetch add
      uint64_t thread_job_index = job_index.fetch_add(1);
      if (thread_job_index < n_files){
        path current_file = signalalign_output_files[thread_job_index];
        AlignmentFile af(current_file.string(), rna);
        vector<VariantCall> vc_calls = af.get_variant_calls(ambig_bases, &ambig_bases_map);
        mv.load_variants(&vc_calls);
        read_id_and_variants = make_tuple(af.read_id, vc_calls);
        variant_queue.push(read_id_and_variants);
        if (verbose) {
//      cout << current_file << "\n";
          // Print status update to stdout
          cerr << "\33[2K\rParsed: " << current_file << flush;
        }
      }
    }
  } catch(...){
    globalExceptionPtr = std::current_exception();
  }
  variant_queue.stop();
}

/**
 * Get the max number of variants given in a map<string, string>
 *
 * @param ambig_bases_map: ambig bases map
 */
uint64_t get_max_n_variants_in_model_file(std::map<string, string>& ambig_bases_map){
  uint64_t max_n_variants = 0;
  uint64_t curr_var_len;
  for (auto const& x : ambig_bases_map)
  {
    curr_var_len = x.second.length();
    if (max_n_variants < curr_var_len){
      max_n_variants = curr_var_len;
    }
  }
  return max_n_variants;
}

/**
 * Worker which reads and writes variants from a queue.
 *
 * @param variant_queue: instance of ConcurrentQueue object made up of vector of variant calls
 * @param max_n_variants: the max number of variants and thus the max number of columns for each row
 * @param output_file: path to output file
 */
void write_tsv_file_worker(
    ConcurrentQueue<tuple<string, vector<VariantCall>>>& variant_queue,
    uint64_t& max_n_variants,
    path& output_file){
  std::ofstream my_file(output_file.string());
  if (my_file.is_open())
  {
    my_file << "read_id,contig,reference_index,strand,variants";
    for (uint64_t i=0; i < max_n_variants; i++){
      my_file << ",prob" << to_string(i+1);
    }
    my_file << endl;

    tuple<string, vector<VariantCall>> nvc;
    vector<VariantCall> vc;
    string read_id;
    uint64_t delta = 0;
    while (variant_queue.wait_and_pop(nvc)){
      read_id = get<0>(nvc);
      vc = get<1>(nvc);
      for (const auto& variant: vc){
        my_file << read_id << "," << variant.contig << "," << variant.reference_index << "," << variant.strand << "," << variant.bases;
        for (auto &prob: variant.normalized_probs){
          my_file << "," << prob;
        }
        delta = max_n_variants - variant.normalized_probs.size();
        for (uint64_t j=0; j < delta; j++){
          my_file << ",";
        }
        my_file << endl;
      }
    }
  } else{
    cout << "Unable to open file: " << output_file.string() << "\n";
  }
}


/**
 Dump signalalign "full" variant calls into a csv

 @param sa_output_paths: vector of paths to sa files
 @param output_file_path: path to output bed file
 @param ambig_bases: possible ambiguous bases to search for
 @param n_locks: number of locks for writing to common data structure
 @param n_threads: number of threads to process files
 @param ambig_model: path to ambig model if not using default
 @param verbose: boolean option to output file names as they are being processed (not helpful)
*/
void dump_signalalign_variant_calls(vector<string> &sa_output_paths,
                                    string &output_file_path,
                                    string ambig_bases,
                                    uint64_t n_threads=2,
                                    uint64_t n_locks=2,
                                    bool rna=false,
                                    string ambig_model = "",
                                    bool verbose=true,
                                    bool overwrite=false) {
  path output_file(output_file_path);
  path output_tsv_file = change_extension(output_file_path, "csv");
  if (!overwrite){
    throw_assert(!exists(output_file),
        output_file_path+" already exists: overwrite to true")
    throw_assert(!exists(output_tsv_file),
        output_tsv_file.string()+" already exists: overwrite to true")
  }
// create ambig model
  throw_assert(exists(ambig_model), ambig_model+" does not exist")
  map<string, string> ambig_bases_map = create_ambig_bases2(ambig_model);
  uint64_t max_n_variants = get_max_n_variants_in_model_file(ambig_bases_map);
// check tsv files
  vector<path> all_tsvs = filter_emtpy_files(sa_output_paths, ".tsv");
  throw_assert(!all_tsvs.empty(), "There are no valid .tsv files")
  auto number_of_files = (int64_t) all_tsvs.size();
// create thread safe queue
  ConcurrentQueue<tuple<string, vector<VariantCall>>> variant_queue;
//  create marginalize variants
  MarginalizeVariants mv(n_locks);
//  get kmer length
  atomic<uint64_t> job_index(0);
  vector<thread> threads;
  globalExceptionPtr = nullptr;
  // Launch threads
  if (n_threads == 1){
    thread t1(get_variants_worker,
           ref(all_tsvs),
           ref(mv),
           ref(variant_queue),
           ref(job_index),
           ref(number_of_files),
           ref(verbose),
           ref(rna),
           ref(ambig_bases),
           ref(ambig_bases_map));
    t1.join();
  } else {
    for (uint64_t i=0; i<n_threads; i++){
      threads.emplace_back(thread(get_variants_worker,
                                  ref(all_tsvs),
                                  ref(mv),
                                  ref(variant_queue),
                                  ref(job_index),
                                  ref(number_of_files),
                                  ref(verbose),
                                  ref(rna),
                                  ref(ambig_bases),
                                  ref(ambig_bases_map)));
    }
//    threads.emplace_back(thread(write_tsv_file_worker,
//                                ref(variant_queue),
//                                ref(max_n_variants),
//                                ref(output_tsv_file)));
    // Wait for threads to finish
    for (auto& t: threads){
      t.join();
    }
  }
  if (globalExceptionPtr){
    std::rethrow_exception(globalExceptionPtr);
  }
  write_tsv_file_worker(variant_queue, max_n_variants, output_tsv_file);
  mv.write_to_file(output_file);
  if (verbose){
    cerr << "\n" << flush;
  }
}

// Getopt
//
#define SUBPROGRAM "sa2bed"
#define SA_TO_BED_VERSION "0.0.1"
#define THIS_NAME "embed"
#define PACKAGE_BUGREPORT2 "None"

static const char *SA2BED_ALIGNMENT_VERSION_MESSAGE =
    SUBPROGRAM " Version " SA_TO_BED_VERSION "\n";

static const char *SA2BED_ALIGNMENT_USAGE_MESSAGE =
    "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --alignment_files ALIGNMENT_FILES --positions_file POSITIONS_FILE --output_dir OUTPUT_DIR\n"
    "Collapses Alignment files into a bed file of called variants or modifications.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -d, --alignment_dir=DIR              directory of signalalign alignment files\n"
    "  -a, --ambig_model=PATH               path to ambig model for ambig char encoding\n"
    "  -o, --output=PATH                    path and name of output bed file\n"
    "  -c, --ambig_chars=NUCLEOTIDES        a string containing all ambiguous characters (default is P which corresponds to cytosine and 5methylcytosine \n"
    "  -t, --threads=NUMBER                 number of threads\n"
    "  -l, --locks=NUMBER                   number of locks for multithreading\n"
    "  -r, --rna                            set if rna reads\n"

    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static uint64_t num_locks = 10000;
static std::string alignment_files;
static std::string output;
static std::string ambig_chars;
static unsigned int threads;
static std::string ambig_model;
static bool rna=false;
static bool overwrite=false;
}

static const char* shortopts = "a:t:c:o:d:r:l:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "alignment_files",  required_argument, nullptr, 'd' },
    { "output",           required_argument, nullptr, 'o' },
    { "ambig_chars",      optional_argument, nullptr, 'c' },
    { "ambig_model",      optional_argument, nullptr, 'a' },
    { "locks",            optional_argument, nullptr, 'l' },
    { "threads",          optional_argument, nullptr, 't' },
    { "help",             no_argument,       nullptr, OPT_HELP },
    { "rna",              no_argument,       nullptr, 'r' },
    { "overwrite",        no_argument,       nullptr, 'b'},
    { "version",          no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

void parse_sa2bed_main_options(int argc, char** argv)
{
  bool die = false;
  char c;
  for (int i; (i = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
    c = (char) i;
    std::istringstream arg(optarg != nullptr ? optarg : "");
    switch (c) {
      case 'd': arg >> opt::alignment_files; break;
      case 'c': arg >> opt::ambig_chars; break;
      case 'o': arg >> opt::output; break;
      case 't': arg >> opt::threads; break;
      case 'l': arg >> opt::num_locks; break;
      case 'a': arg >> opt::ambig_model; break;
      case 'r': opt::rna = true; break;
      case 'b': opt::overwrite = true; break;
      case 'v': opt::verbose++; break;
      case OPT_HELP:
        std::cout << SA2BED_ALIGNMENT_USAGE_MESSAGE;
        exit(EXIT_SUCCESS);
      case OPT_VERSION:
        std::cout << SA2BED_ALIGNMENT_VERSION_MESSAGE;
        exit(EXIT_SUCCESS);
      default:
        string error = ": unreconized argument -";
        error += c;
        error += " \n";
        std::cerr << SUBPROGRAM + error;
        exit(EXIT_FAILURE);
    }
  }

  if (argc - optind > 0) {
    std::cerr << SUBPROGRAM ": too many arguments\n";
    die = true;
  }

  if(opt::alignment_files.empty()) {
    std::cerr << SUBPROGRAM ": --alignment_files file must be provided\n";
    die = true;
  }
  if(opt::output.empty()) {
    std::cerr << SUBPROGRAM ": --output file must be provided\n";
    die = true;
  }
  if(opt::output.empty()) {
    std::cerr << SUBPROGRAM ": --ambig_chars must be provided\n";
    die = true;
  }

  if (die)
  {
    std::cout << "\n" << SA2BED_ALIGNMENT_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}


int sa2bed_main(int argc, char** argv)
{
  parse_sa2bed_main_options(argc, argv);
  path input_dir(opt::alignment_files);
  vector<string> all_tsvs;
  string extension = ".tsv";
  for (auto &i: list_files_in_dir(input_dir, extension)) {
    all_tsvs.push_back(i.string());
  }

  auto bound_funct = bind(dump_signalalign_variant_calls,
      all_tsvs,
      opt::output,
      opt::ambig_chars,
      opt::threads,
      opt::num_locks,
      opt::rna,
      opt::ambig_model,
      false,
      opt::overwrite);
  string funct_time = get_time_string(bound_funct);
  cout << funct_time;
  return EXIT_SUCCESS;
}
