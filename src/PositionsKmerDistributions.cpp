//
// Created by Andrew Bailey on 6/2/20.
//

// embed lib
#include "PositionsKmerDistributions.hpp"
#include "PerPositionKmers.hpp"
#include "AmbigModel.hpp"
// boost lib
#include <boost/filesystem.hpp>
// std lib
#include <getopt.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <functional>
#include<algorithm>


using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;
static std::exception_ptr globalExceptionPtr = nullptr;

void write_kmer_distribution_file(const vector<uint64_t>& data, const path& output_path, const float& min, const float& max,
                                  const uint64_t& steps, const float& threshold){
  boost::filesystem::ofstream myfile(output_path);
  throw_assert(data.size() == steps,
               "Number of steps:" + to_string(steps) + " does not equal length of the data: " + to_string(steps))
  if (myfile.is_open())
  {
    myfile << min << ',' << max << ',' << steps << ',' << threshold << '\n';

    for (uint64_t i=0; i < steps-1; ++i){
      myfile << data[i] << ',';
    }
    myfile << data[steps-1] << '\n';
    myfile.close();
  } else{
    cout << "Unable to open file: " << output_path.string() << "\n";
  }
}

void write_plot_kmer_dist_file(const path& output_path, const float& min, const float& max,
                               const uint64_t& steps, const float& threshold,
                               const vector<pair<string, vector<uint64_t>>>& data){
  boost::filesystem::ofstream myfile(output_path);
  if (myfile.is_open())
  {
    myfile << min << ',' << max << ',' << steps << ',' << threshold << ',' << data.size() <<'\n';
    for (auto &my_pair: data){
      throw_assert(my_pair.second.size() == steps,
                   "Number of steps:" + to_string(steps) + " does not equal length of the data: " + to_string(steps))
      myfile << my_pair.first << '\n';
      for (uint64_t i=0; i < steps-1; ++i){
        myfile << my_pair.second[i] << ',';
      }
      myfile << my_pair.second[steps-1] << '\n';
    }
    myfile.close();
  } else{
    cout << "Unable to open file: " << output_path.string() << "\n";
  }
}


/**
 Get kmer positional distributions based on the positions file

 @param positions_file_path: path to sa files
 @param output_dir: path to output directory
 @param reference: path to reference file
 @param event_file: path to event file
 @param nanopore_strand: "t" or "c" nanopore strand to process
 @param min: low bound on histogram
 @param max: high bound on histogram
 @param size: number of bars in histogram
 @param min_prob_threshold: minimum probability to process
*/
void get_kmer_distributions_by_position(const string &positions_file_path,
                                        const string &output_dir,
                                        const string &reference,
                                        const string &event_file,
                                        const uint64_t& n_threads,
                                        const string &ambig_model,
                                        const string &nanopore_strand,
                                        const float& min,
                                        const float& max,
                                        const uint64_t& size,
                                        const float& min_prob_threshold) {
//  check output dirs
  path output_dir_path(output_dir);
  throw_assert(exists(output_dir_path), output_dir+" does not exist")
  path position_output_dir = output_dir_path / "positions";
  create_directory(position_output_dir);
//  initialize per-position dataset using info from reference
  AmbigModel am(ambig_model);
  ReferenceHandler rh(reference);
  PositionsFile pf(positions_file_path);
  EventDataHandler edh(rh, event_file);
  uint64_t kmer_length = edh.get_kmer_length();
//  kmer histograms of correct size
  vector<uint64_t> kmer_hist(size, 0);
  vector<uint64_t> pos_hist(size, 0);
  string ref_pos;
  path pos_specific_dir;
  path pos_file;
  path kmer_file;
//  string prev_contig = "";
//  uint64_t prev_pos = -1;
//  string prev_strand = "";

  for (auto &line: pf.iterate()){
    cout << line.contig << " " << line.position << '\n';
    ref_pos = rh.get_reference_sequence(line.contig, line.position, line.position+1);
//    if ((prev_contig == line.contig) & (prev_strand == line.strand) & (line.position - prev_pos) >= kmer_length) {
//    }
    pos_specific_dir = position_output_dir / path(line.contig+line.strand+to_string(line.position));
    create_directory(pos_specific_dir);
    throw_assert( ref_pos == line.change_from,
                  "Reference position does not match positions file. Check both.\ncontig:" + line.contig +
                  "\n" + "strand:" + line.strand + "\n" + "position:" + to_string(line.position) + "\n" +
                  "change_from:" + line.change_from + "\n" + "ref base:" + ref_pos + "\n");
    try{
      for (uint64_t i=0; i < kmer_length; ++i){
        vector<pair<string, vector<uint64_t>>> data;
        Position& pos = edh.get_position(line.contig, line.strand, nanopore_strand, line.position-i);
        set<string> kmers = pos.get_kmer_strings();
        set<string> canonical_kmers = am.get_canonical_kmers(kmers);
        kmers.insert(canonical_kmers.begin(), canonical_kmers.end());

        for (auto &k: kmers){
          if (edh.has_kmer(k)) {
            Kmer& kmer = edh.get_kmer(k);
            kmer_hist = kmer.get_hist(min, max, size, min_prob_threshold);
            data.push_back(make_pair(k, kmer_hist));
            for (auto &pos_k: kmer.pos_kmer_map){
              ContigStrandPosition csp = kmer.split_pos_kmer_map_key(pos_k.first);
              data.push_back(make_pair(csp.contig+"_"+csp.strand+"_"+to_string(csp.position)+"_"+k, pos_k.second->get_hist(min, max, size, min_prob_threshold)));
            }
          }
        }
        pos_file = pos_specific_dir / path(line.contig+"_"+line.strand+"_"+to_string(line.position-i)+".csv");
        write_plot_kmer_dist_file(pos_file, min, max, size, min_prob_threshold, data);
      }
//    prev_contig = line.contig;
//    prev_pos = line.position;
//    prev_strand = line.strand;
    } catch(std::runtime_error& e){
      cout << e.what() << "\n";
    }
  }
}


// Getopt
//
#define SUBPROGRAM "get_kmer_distributions"
#define GET_KMER_DISTRIBUTIONS_VERSION "0.0.1"
#define THIS_NAME "embed"
#define PACKAGE_BUGREPORT2 "None"

static const char *GET_KMER_DISTRIBUTIONS_VERSION_MESSAGE =
    SUBPROGRAM " Version " GET_KMER_DISTRIBUTIONS_VERSION "\n";

static const char *GET_KMER_DISTRIBUTIONS_USAGE_MESSAGE =
    "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --positions_file POSITIONS_FILE "
                                       "--output OUTPUT_DIR --event_file EVENT_FILE --reference REFERENCE "
                                       "--ambig_model AMBIG_MODEL_PATH\n"
    "Output kernal density estimates for kmers covering positions from a positions file.\n"
    "\n"
    "  -v, --verbose                        display verbose output\n"
    "      --version                        display version\n"
    "      --help                           display this help and exit\n"
    "  -p, --positions_file=DIR             path to positions file\n"
    "  -o, --output=PATH                    path to output directory\n"
    "  -e, --event_file=PATH                path to .event file generated from SplitByRefPosition\n"
    "  -r, --reference=PATH                 reference sequence (fa format)\n"
    "  -a, --ambig_model=PATH               ambig_model path \n"
    "  -m, --min_prob=PROB                  Minimum probability threshold \n"

    "  -t, --threads=NUMBER                 number of threads\n"
    "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
static unsigned int verbose;
static std::string positions_file;
static std::string output;
static uint64_t threads = 1;
static string reference;
static string event_file;
static string ambig_model;
static float min_prob = 0.0;
}

static const char* shortopts = "p:t:m:o:r:a:e:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",          no_argument,       nullptr, 'v' },
    { "positions_file",   required_argument, nullptr, 'p' },
    { "output",           required_argument, nullptr, 'o' },
    { "reference",        required_argument, nullptr, 'r' },
    { "ambig_model",      required_argument, nullptr, 'a' },
    { "event_file",       required_argument, nullptr, 'e' },
    { "threads",          optional_argument, nullptr, 't' },
    { "min_prob",         optional_argument, nullptr, 'm' },
    { "help",             no_argument,       nullptr, OPT_HELP },
    { "version",          no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

void parse_get_kmer_distributions_main_options(int argc, char** argv)
{
  bool die = false;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
    std::istringstream arg(optarg != nullptr ? optarg : "");
    switch (c) {
      case 'p': arg >> opt::positions_file; break;
      case 'o': arg >> opt::output; break;
      case 't': arg >> opt::threads; break;
      case 'e': arg >> opt::event_file; break;
      case 'r': arg >> opt::reference; break;
      case 'a': arg >> opt::ambig_model; break;
      case 'm': arg >> opt::min_prob; break;
      case 'v': opt::verbose++; break;
      case OPT_HELP:
        std::cout << GET_KMER_DISTRIBUTIONS_USAGE_MESSAGE;
        exit(EXIT_SUCCESS);
      case OPT_VERSION:
        std::cout << GET_KMER_DISTRIBUTIONS_VERSION_MESSAGE;
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

  if(opt::positions_file.empty()) {
    std::cerr << SUBPROGRAM ": a --positions_file file must be provided\n";
    die = true;
  }
  if(opt::output.empty()) {
    std::cerr << SUBPROGRAM ": a --output file must be provided\n";
    die = true;
  }
  if(opt::event_file.empty()) {
    std::cerr << SUBPROGRAM ": a --event_file file must be provided\n";
    die = true;
  }
  if(opt::ambig_model.empty()) {
    std::cerr << SUBPROGRAM ": a --ambig_model file must be provided\n";
    die = true;
  }
  if(opt::reference.empty()) {
    std::cerr << SUBPROGRAM ": a --reference string must be provided\n";
    die = true;
  }
  if (die)
  {
    std::cout << "\n" << GET_KMER_DISTRIBUTIONS_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}


int get_kmer_distributions_main(int argc, char** argv)
{
  parse_get_kmer_distributions_main_options(argc, argv);
  auto bound_funct = bind(get_kmer_distributions_by_position,
                          opt::positions_file,
                          opt::output,
                          opt::reference,
                          opt::event_file,
                          opt::threads,
                          opt::ambig_model,
                          "t",
                          0.0,
                          200.0,
                          2000,
                          opt::min_prob);
  cout << get_time_string(bound_funct);

  return EXIT_SUCCESS;
}

