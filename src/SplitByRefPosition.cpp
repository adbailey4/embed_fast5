//
// Created by Andrew Bailey on 2019-07-24.
//

// embed lib
#include "SplitByRefPosition.hpp"
#include "PerPositionKmers.hpp"
#include "EmbedUtils.hpp"
// boost lib
#include <boost/filesystem.hpp>
// std lib
#include <getopt.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <functional>

using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;
static std::exception_ptr globalExceptionPtr = nullptr;


/**
 * worker which parses "full" signalalign file by position
 *
 * @param signalalign_output_files: reference to vector of signalalign files
 * @param ppk: PerPositonKmers class object
 * @param job_index: atomic index for selecting output files to process
 * @param n_files: max number of files to process
 * @param verbose: option for printing files processed
 * @param rna: boolean option if reads are rna
 */
void per_position_worker(
    vector<path>& signalalign_output_files,
    PerPositionKmers& ppk,
    atomic<uint64_t>& job_index,
    uint64_t& n_files,
    bool& verbose,
    bool& rna,
    ProgressBar& pb) {
  uint64_t step = floor(n_files / 100) + 1;
  try {
    tuple<string, vector<VariantCall>> read_id_and_variants;
    while (job_index < n_files and !globalExceptionPtr) {
      // Fetch add
      uint64_t thread_job_index = job_index.fetch_add(1);
      if (thread_job_index < n_files){
        path current_file = signalalign_output_files[thread_job_index];
        AlignmentFile af(current_file.string(), rna);
        ppk.process_alignment(af);
        if (verbose) {
          // Print status update to stdout
          cerr << "\33[2K\rParsed: " << current_file << flush;
        }
        if (thread_job_index % step == 0){
          pb.write((double)thread_job_index / (double)n_files);
        }
      }
    }
  } catch(...){
    globalExceptionPtr = std::current_exception();
  }
}

/**
 Split full signalalign output by reference position and write a file with top n most probable kmers

 @param sa_input_dir: path to sa files
 @param output_file_path: path to output bed file
 @param ambig_bases: possible ambiguous bases to search for
 @param num_locks: number of locks for writing to common data structure
 @param n_threads: number of threads to process
 @return tuple of uint64_t's [hours, minutes, seconds, microseconds]
*/
void split_signal_align_by_ref_position(const vector<string> &sa_input_dir,
                                        string &output_file_path,
                                        string reference,
                                        uint64_t num_locks,
                                        uint64_t n_threads,
                                        bool verbose,
                                        bool rna,
                                        bool two_d,
                                        set<char> alphabet) {
  //  check output file does not exist
  path output_file(output_file_path);
  throw_assert(!exists(output_file), output_file_path+" already exists")

//  process all tsvs from directories
  vector<path> all_tsvs;
  string extension = ".tsv";
  cout << "Reading in Files..\n";
  for (auto &in_dir: sa_input_dir){
    path input_dir(in_dir);
    for (auto &i: list_files_in_dir(input_dir, extension)) {
      all_tsvs.push_back(i);
    }
  }
  uint64_t number_of_files = all_tsvs.size();
//  initialize per-position dataset using info from reference
  ReferenceHandler rh(reference);
  PerPositionKmers ppk(rh, alphabet, AlignmentFile(all_tsvs[0].string()).get_k(), num_locks, two_d);
  //  creat job index, threads and reset exception pointer
  atomic<uint64_t> job_index(0);
  vector<thread> threads;
  globalExceptionPtr = nullptr;
  cout << "\33[2K\rStarting threads..\n ";
  // Launch threads
  ProgressBar progress{std::cout, 70u, "Working"};
  for (uint64_t i=0; i<n_threads; i++){
    threads.emplace_back(thread(per_position_worker,
                                ref(all_tsvs),
                                ref(ppk),
                                ref(job_index),
                                ref(number_of_files),
                                ref(verbose),
                                ref(rna),
                                ref(progress)));
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
  progress.~ProgressBar();
  cout << "\33[2K\rWriting to file.. \n ";
  ppk.write_to_file(output_file);
}

// Getopt
//
#define SUBPROGRAM "split_by_position"
#define SPLIT_BY_REF_VERSION "0.0.1"
#define THIS_NAME "embed"
#define PACKAGE_BUGREPORT2 "None"

static const char *SPLIT_BY_REF_VERSION_MESSAGE =
    SUBPROGRAM " Version " SPLIT_BY_REF_VERSION "\n";

static const char *SPLIT_BY_REF_USAGE_MESSAGE =
    "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --alignment_files ALIGNMENT_FILES --output OUTPUT_PATH\n"
    "Filters alignment files into per reference position kmer tables.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -a, --alignment_files=DIR            directory of signalalign alignment files\n"
    "  -o, --output=PATH                    path and name of output bed file\n"
    "  -t, --threads=NUMBER                 number of threads\n"
    "  -l, --locks=NUMBER                   number of locks for multithreading\n"
    "  -r, --reference=PATH                 reference sequence (fa format)\n"
    "  -c, --alphabet=PATH                  characters that make up alphabet\n"
    "  --rna                                boolean option if reads are rna\n"
    "  --two_d                              boolean option if reads are 2d\n"
    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static uint64_t num_locks = 10000;
vector<string> alignment_files;
static std::string output;
static unsigned int threads = 1;
static string reference;
static bool rna=false;
static bool two_d=false;
static string alphabet;

}

static const char* shortopts = "a:t:o:r:l:d:b:c:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "alignment_files",  required_argument, nullptr, 'a' },
    { "output",           required_argument, nullptr, 'o' },
    { "reference",        required_argument, nullptr, 'r' },
    { "alphabet",         required_argument, nullptr, 'c' },
    { "locks",            optional_argument, nullptr, 'l' },
    { "threads",          optional_argument, nullptr, 't' },
    { "rna",              no_argument,       nullptr, 'b' },
    { "two_d",            no_argument,       nullptr, 'd' },
    { "help",             no_argument,       nullptr, OPT_HELP },
    { "version",          no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

void parse_split_by_ref_main_options(int argc, char** argv)
{
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
    std::istringstream arg(optarg != nullptr ? optarg : "");
    switch (c) {
      case 'a': opt::alignment_files.push_back(arg.str()); break;
      case 'o': arg >> opt::output; break;
      case 't': arg >> opt::threads; break;
      case 'l': arg >> opt::num_locks; break;
      case 'r': arg >> opt::reference; break;
      case 'c': arg >> opt::alphabet; break;
      case 'b': opt::rna = true; break;
      case 'd': opt::two_d = true; break;
      case 'v': opt::verbose++; break;
      case OPT_HELP:
        std::cout << SPLIT_BY_REF_USAGE_MESSAGE;
        exit(EXIT_SUCCESS);
      case OPT_VERSION:
        std::cout << SPLIT_BY_REF_VERSION_MESSAGE;
        exit(EXIT_SUCCESS);
      default:
        string error = ": unreconized argument -";
        error += c;
        error += " \n";
        std::cerr << SUBPROGRAM + error;
        exit(EXIT_FAILURE);
    }
  }

  if (argc - (optind+1) > 0) {
    std::cerr << SUBPROGRAM ": too many arguments\n";
    die = true;
  }

  if(opt::alignment_files.empty()) {
    std::cerr << SUBPROGRAM ": a --alignment_files file must be provided\n";
    die = true;
  }
  if(opt::output.empty()) {
    std::cerr << SUBPROGRAM ": a --output file must be provided\n";
    die = true;
  }
  if(opt::alphabet.empty()) {
    std::cerr << SUBPROGRAM ": a --alphabet string must be provided\n";
    die = true;
  }
  if (die)
  {
    std::cout << "\n" << SPLIT_BY_REF_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}


int split_by_ref_main(int argc, char** argv)
{
  parse_split_by_ref_main_options(argc, argv);
  set<char> alphabet = string_to_char_set(opt::alphabet);
  auto bound_funct = bind(split_signal_align_by_ref_position,
                          opt::alignment_files,
                          opt::output,
                          opt::reference,
                          opt::num_locks,
                          opt::threads,
                          opt::verbose,
                          opt::rna,
                          opt::two_d,
                          alphabet);
  string funct_time = get_time_string(bound_funct);
  cout << funct_time;

  return EXIT_SUCCESS;
}

