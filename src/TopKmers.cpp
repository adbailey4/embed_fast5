//
// Created by Andrew Bailey on 2019-06-17.
//

#include "TopKmers.hpp"
#include "AssignmentFile.hpp"
#include "AlignmentFile.hpp"
#include "MaxKmers.hpp"
#include <getopt.h>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;

/**
 * Wrapper over generate_master_kmer_table so we can infer file type by looking at the number of columns.
 * Generate the master table by parsing event table files and outputting the top n kmers to a files
 *
 * @param event_table_files: vector of strings of files
 * @param output_file: path to output file
 * @param log_file: path to log file
 * @param heap_size: number of max kmers to keep for each kmer
 * @param alphabet: alphabet used to generate kmers
 * @param n_threads: set number of threads to use: default 2
 * @param verbose: print out files as they are being processed
 */
void generate_master_kmer_table_wrapper(vector<string> event_table_files,
                                        string &output_file,
                                        string &log_file,
                                        uint64_t heap_size,
                                        string &alphabet,
                                        double min_prob,
                                        uint64_t n_threads,
                                        bool verbose,
                                        bool write_full) {
  uint64_t n_col = number_of_columns(event_table_files[0]);
  throw_assert(n_col == 16 or n_col == 4,
               "Incorrect number of columns in tsv: " + event_table_files[0])
  if (n_col == 4) {
    generate_master_kmer_table<AssignmentFile, eventkmer>(event_table_files, output_file, log_file,
                                                          alphabet, heap_size, min_prob, n_threads,
                                                          verbose, write_full);

  } else if (n_col == 16) {
    generate_master_kmer_table<AlignmentFile, FullSaEvent>(event_table_files, output_file, log_file,
                                                           alphabet, heap_size, min_prob, n_threads,
                                                           verbose, write_full);
  }
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
    "  -d, --event_directory=DIR            directory of signalalign event table files\n"
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
static double min_prob = 0.0;
}

static const char* shortopts = "a:d:s:t:o:m:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "event_directory",  required_argument, nullptr, 'd' },
    { "output_dir",       required_argument, nullptr, 'o' },
    { "heap_size",        required_argument, nullptr, 's' },
    { "alphabet",         required_argument, nullptr, 'a' },
    { "min_prob",         required_argument, nullptr, 'm' },
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
      case 'm': arg >> opt::min_prob; break;
      case 's': arg >> opt::heap_size; break;
      case 'v': opt::verbose++; break;
      case OPT_HELP:
        std::cout << TOP_KMER_USAGE_MESSAGE;
        exit(EXIT_SUCCESS);
      case OPT_VERSION:
        std::cout << TOP_KMER_VERSION_MESSAGE;
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


auto top_kmers_main(int argc, char** argv) -> int
{
  parse_top_kmers_main_options(argc, argv);

#ifndef H5_HAVE_THREADSAFE
  if(opt::num_threads > 1) {
    fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
    fprintf(stderr, "Please recompile built-in libhdf5 or run with -t 1\n");
    exit(1);
  }
#endif
  vector<string> all_files;
  path p(opt::assignment_dir);
  string ext = ".tsv";
  for (auto& file: list_files_in_dir(p, ext)){
    all_files.push_back(file.string());
  }
  path out_dir_path(opt::output_dir);
  string output_file = (out_dir_path / "buildAlignment.tsv").string();
  string log_file = (out_dir_path / "bA_log.tsv").string();

  generate_master_kmer_table_wrapper(all_files,
                                     output_file,
                                     log_file,
                                     opt::heap_size,
                                     opt::alphabet,
                                     opt::min_prob,
                                     opt::threads,
                                     opt::verbose, true);

  return EXIT_SUCCESS;
}




