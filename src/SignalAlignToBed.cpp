//
// Created by Andrew Bailey on 2019-07-14.
//

#include "SignalAlignToBed.hpp"
#include "EmbedUtils.hpp"
#include "MarginalizeVariants.hpp"
#include <getopt.h>
#include <iostream>
#include <boost/filesystem.hpp>
#include "omp.h"
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;

/**
 Convert signalalign "full" outputs into a bed file

 @param sa_input_dir: path to sa files
 @param output_file_path: path to output bed file
 @param ambig_bases: possible ambiguous bases to search for
 @param num_locks: number of locks for writing to common data structure
 @return tuple of uint64_t's [hours, minutes, seconds, microseconds]
*/
void signalalign_to_bed(string& sa_input_dir, string& output_file_path, string ambig_bases, uint64_t num_locks){
  path input_dir(sa_input_dir);
  path output_file(output_file_path);
  throw_assert(!exists(output_file), output_file_path+" already exists")
  std::map<string, string> ambig_bases_map = create_ambig_bases();
  MarginalizeVariants mv(num_locks);
  vector<path> all_tsvs;
  string extension = ".tsv";
  for (auto &i: list_files_in_dir(input_dir, extension)) {
    all_tsvs.push_back(i);
  }
  uint64_t number_of_files = all_tsvs.size();
//  uint64_t number_of_files = 4;
#pragma omp parallel for shared(all_tsvs, mv, number_of_files, ambig_bases, ambig_bases_map, cout)
  for(uint64_t i=0; i < number_of_files; i++) {
    cout << all_tsvs[i] << "\n";
    path current_file = all_tsvs[i];
//    std::fstream file(current_file.string());
    AlignmentFile af(current_file.string());
    vector<VariantCall> vc_calls = af.get_variant_calls(ambig_bases, &ambig_bases_map);
    mv.load_variants(&vc_calls);
  }
  mv.write_to_file(output_file);
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
    "Converts alignment files into assignment files with kmers covering certain positions.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -a, --alignment_files=DIR            directory of signalalign alignment files\n"
    "  -o, --output=PATH                    path and name of output bed file\n"
    "  -c, --ambig_chars=NUCLEOTIDES        a string containing all ambiguous characters (default is P which corresponds to cytosine and 5methylcytosine \n"
    "  -t, --threads=NUMBER                 number of threads\n"
    "  -l, --locks=NUMBER                   number of locks for multithreading\n"

    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static uint64_t num_locks = 10000;
static std::string alignment_files;
static std::string output;
static std::string ambig_chars = "P";
static unsigned int threads;
}

static const char* shortopts = "a:t:c:o:l:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "alignment_files",  required_argument, nullptr, 'a' },
    { "output",           required_argument, nullptr, 'o' },
    { "ambig_chars",      optional_argument, nullptr, 'c' },
    { "locks",            optional_argument, nullptr, 'l' },
    { "threads",          optional_argument, nullptr, 't' },
    { "help",             no_argument,       nullptr, OPT_HELP },
    { "version",          no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

void parse_sa2bed_main_options(int argc, char** argv)
{
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
    std::istringstream arg(optarg != nullptr ? optarg : "");
    switch (c) {
      case 'a': arg >> opt::alignment_files; break;
      case 'c': arg >> opt::ambig_chars; break;
      case 'o': arg >> opt::output; break;
      case 't': arg >> opt::threads; break;
      case 'l': arg >> opt::num_locks; break;

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
    std::cerr << SUBPROGRAM ": a --alignment_files file must be provided\n";
    die = true;
  }
  if(opt::output.empty()) {
    std::cerr << SUBPROGRAM ": a --output file must be provided\n";
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
  omp_set_num_threads(opt::threads); // Use 4 threads for all consecutive parallel regions
  auto bound_funct = bind(signalalign_to_bed, opt::alignment_files, opt::output, opt::ambig_chars, opt::num_locks);
  string funct_time = get_time_string(bound_funct);
  cout << funct_time;

  return EXIT_SUCCESS;
}
