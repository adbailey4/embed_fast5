//
// Created by Andrew Bailey on 2019-07-24.
//

#ifndef EMBED_FAST5_SRC_SCRIPTS_SPLITBYREFPOSITION_HPP_
#define EMBED_FAST5_SRC_SCRIPTS_SPLITBYREFPOSITION_HPP_

#include "SplitByRefPosition.hpp"
#include "AlignmentFile.hpp"
#include "PerPositonKmers.hpp"
#include "ReferenceHandler.hpp"
#include "EmbedUtils.hpp"
#include <getopt.h>
#include <iostream>
#include <boost/filesystem.hpp>
#include "omp.h"
#include <chrono>

using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;

/**
 Split full signalalign output by reference position and write a file with top n most probable kmers

 @param sa_input_dir: path to sa files
 @param output_file_path: path to output bed file
 @param ambig_bases: possible ambiguous bases to search for
 @param num_locks: number of locks for writing to common data structure
 @return tuple of uint64_t's [hours, minutes, seconds, microseconds]
*/
void split_signal_align_by_ref_position(string& sa_input_dir, string& output_file_path, string reference, uint64_t num_locks){
  path input_dir(sa_input_dir);
  path output_file(output_file_path);
  throw_assert(!exists(output_file), output_file_path+" already exists")
  std::map<string, string> ambig_bases_map = create_ambig_bases();
  vector<path> all_tsvs;
  string extension = ".tsv";
  for (auto &i: list_files_in_dir(input_dir, extension)) {
    all_tsvs.push_back(i);
  }
  uint64_t number_of_files = all_tsvs.size();
  PerPositonKmers ppk();
//  uint64_t number_of_files = 4;
//#pragma omp parallel for shared(all_tsvs, mv, number_of_files, ambig_bases, ambig_bases_map)
  for(uint64_t i=0; i < number_of_files; i++) {
    cout << all_tsvs[i] << "\n";
    path current_file = all_tsvs[i];
//    std::fstream file(current_file.string());
    AlignmentFile af(current_file.string());
//    ppk.split_alignment_file(af);
  }
}

// Getopt
//
#define SUBPROGRAM "split"
#define SPLIT_BY_REF_VERSION "0.0.1"
#define THIS_NAME "embed"
#define PACKAGE_BUGREPORT2 "None"

static const char *SPLIT_BY_REF_VERSION_MESSAGE =
    SUBPROGRAM " Version " SPLIT_BY_REF_VERSION "\n";

static const char *SPLIT_BY_REF_USAGE_MESSAGE =
    "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --alignment_files ALIGNMENT_FILES --output_dir OUTPUT_DIR\n"
    "Filters alignment files into per reference position kmer tables.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -a, --alignment_files=DIR            directory of signalalign alignment files\n"
    "  -o, --output=PATH                    path and name of output bed file\n"
    "  -t, --threads=NUMBER                 number of threads\n"
    "  -l, --locks=NUMBER                   number of locks for multithreading\n"
    "  -r, --reference=NUMBER               reference sequence (fa format)\n"

    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static uint64_t num_locks = 10000;
static std::string alignment_files;
static std::string output;
static unsigned int threads = 1;
static string reference;
}

static const char* shortopts = "a:t:o:r:l:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "alignment_files",  required_argument, nullptr, 'a' },
    { "output",           required_argument, nullptr, 'o' },
    { "reference",        required_argument, nullptr, 'r' },
    { "locks",            optional_argument, nullptr, 'l' },
    { "threads",          optional_argument, nullptr, 't' },
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
      case 'a': arg >> opt::alignment_files; break;
      case 'o': arg >> opt::output; break;
      case 't': arg >> opt::threads; break;
      case 'l': arg >> opt::num_locks; break;
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
  if (die)
  {
    std::cout << "\n" << SPLIT_BY_REF_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}


int split_by_ref_main(int argc, char** argv)
{
  parse_split_by_ref_main_options(argc, argv);
  omp_set_num_threads(opt::threads);
  auto bound_funct = bind(split_signal_align_by_ref_position, opt::alignment_files, opt::output, opt::reference, opt::num_locks);
  string funct_time = get_time_string(bound_funct);
  cout << funct_time;

  return EXIT_SUCCESS;
}

#endif //EMBED_FAST5_SRC_SCRIPTS_SPLITBYREFPOSITION_HPP_
