//
// Created by Andrew Bailey on 2019-06-17.
//

#include "TopKmers.hpp"
#include "AssignmentFile.hpp"
#include "MaxKmers.hpp"
#include "EmbedUtils.hpp"
#include <getopt.h>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;


/**
 * Generate the master assignment table by parsing assignment files and outputting the top n kmers to a files
 *
 * @param assignment_dir: path to assignment files directory
 * @param output_dir: path to output directory where new builtAssignment.tsv will be written
 * @param heap_size: number of max kmers to keep for each kmer
 * @param alphabet: alphabet used to generate kmers
 */

void generate_master_assignment_table(string assignment_dir, string& output_dir, int heap_size, string& alphabet){

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
#pragma omp parallel for shared(array_of_files, mk, number_of_files, cout)
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




