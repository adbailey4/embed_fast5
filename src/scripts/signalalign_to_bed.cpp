//
// Created by Andrew Bailey on 2019-07-14.
//

#include "signalalign_to_bed.hpp"
#include "EmbedUtils.hpp"
#include "BedFile.hpp"
#include <getopt.h>
#include <iostream>
#include <boost/filesystem.hpp>
#include "omp.h"

using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;

void signalalign_to_bed(string& sa_input_dir, string& output_bed_path){
  path input_dir(sa_input_dir);
  path output_bed(output_bed_path);
  throw_assert(!exists(output_bed), output_bed_path+" already exists");
  directory_iterator end_itr;
  vector<path> all_tsvs;
  int counter = 0;
  for (directory_iterator itr(input_dir); itr != end_itr; ++itr) {
//        filter for files that are regular, end with tsv and are not empty
    if (is_regular_file(itr->path()) and itr->path().extension().string() == ".tsv" and getFilesize(itr->path().string()) > 0) {
      all_tsvs.push_back(itr->path());
      if (counter == 0){
//        AlignmentFile af = AlignmentFile(itr->path().string());
//        k = af.k;
        counter += 1;
      }
    }
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
    "Converts alignment files into assignment files with kmers covering certain positions.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -a, --alignment_files=DIR            directory of signalalign alignment files\n"
    "  -o, --output=PATH                path and name of output bed file\n"

    "  -t, --threads=NUMBER                 number of threads\n"

    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static std::string alignment_files;
static std::string output;
static unsigned int threads;
}

static const char* shortopts = "a:t:o:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "alignment_files",  required_argument, nullptr, 'a' },
    { "output",           required_argument, nullptr, 'o' },
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
      case 'o': arg >> opt::output; break;
      case 't': arg >> opt::threads; break;
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
  signalalign_to_bed(opt::alignment_files, opt::output);

  return EXIT_SUCCESS;
}
