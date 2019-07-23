//
// Created by Andrew Bailey on 2019-06-07.
//
#include "FilterAlignments.hpp"
#include "AlignmentFile.hpp"
#include "EmbedUtils.hpp"
#include "omp.h"
#include <getopt.h>



using namespace std;
using namespace boost::filesystem;
using namespace embed_utils;

void filter_alignment_files(string input_reads_dir, const string& positions_file, string output_dir, string& bases){

    path p(input_reads_dir);
    path output_path = make_dir(output_dir);

    directory_iterator end_itr;
//    Get all tsvs to process
    vector<path> all_tsvs;
    int64_t k;
    int counter = 0;
    for (directory_iterator itr(p); itr != end_itr; ++itr) {
//        filter for files that are regular, end with tsv and are not empty
        if (is_regular_file(itr->path()) and itr->path().extension().string() == ".tsv" and getFilesize(itr->path().string()) > 0) {
            all_tsvs.push_back(itr->path());
            if (counter == 0){
                AlignmentFile af(itr->path().string());
                k = af.k;
                counter += 1;
            }
        }
    }

    PositionsFile pf = PositionsFile(positions_file, k);

    int64_t number_of_files = all_tsvs.size();
    path* array_of_files = &all_tsvs[0];
// looping through the files
    #pragma omp parallel for shared(array_of_files, pf, number_of_files)
    for(int64_t i=0; i < number_of_files; i++) {

        path current_file = array_of_files[i];
        cout << current_file << "\n";
        AlignmentFile af(current_file.string());
        path output_file = output_path / current_file.filename();
//        if (current_file.filename().string() == "0a4e473d-4713-4c7f-9e18-c465ea6d5b8c.sm.forward.tsv"){
        af.filter_by_positions(&pf, output_file, bases);
//        }
    }

}


// Getopt
//
#define SUBPROGRAM "filter_alignments"
#define FILTER_VERSION "0.0.1"
#define THIS_NAME "filter"
#define PACKAGE_BUGREPORT2 "None"

static const char *FILTER_ALIGNMENT_VERSION_MESSAGE =
        SUBPROGRAM " Version " FILTER_VERSION "\n";

static const char *FILTER_ALIGNMENT_USAGE_MESSAGE =
        "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --alignment_files ALIGNMENT_FILES --positions_file POSITIONS_FILE --output_dir OUTPUT_DIR\n"
        "Converts alignment files into assignment files with kmers covering certain positions.\n"
        "\n"
        "  -v, --verbose                        display verbose output\n"
        "      --version                        display version\n"
        "      --help                           display this help and exit\n"
        "  -a, --alignment_files=DIR            directory of signalalign alignment files\n"
        "  -p, --positions_file=FILE            path to positions file\n"
        "  -o, --output_dir=DIR                 path to directory to output filtered alignment files\n"
        "  -t, --threads=NUMBER                 number of threads\n"
        "  -n, --bases=BASES                    nucleotides to filter out unless in positions file\n"

        "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string alignment_files;
    static std::string positions_file;
    static std::string output_dir;
    static unsigned int threads;
    static std::string bases;
}

static const char* shortopts = "a:p:t:o:n:vh";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
        { "verbose",          no_argument,       nullptr, 'v' },
        { "alignment_files",  required_argument, nullptr, 'a' },
        { "positions_file",   required_argument, nullptr, 'p' },
        { "output_dir",       required_argument, nullptr, 'o' },
        { "threads",          optional_argument, nullptr, 't' },
        { "bases",            optional_argument, nullptr, 'n' },
        { "help",             no_argument,       nullptr, OPT_HELP },
        { "version",          no_argument,       nullptr, OPT_VERSION },
        { nullptr, 0, nullptr, 0 }
};

void parse_filter_main_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
            case 'a': arg >> opt::alignment_files; break;
            case 'p': arg >> opt::positions_file; break;
            case 'o': arg >> opt::output_dir; break;
            case 't': arg >> opt::threads; break;
            case 'n': arg >> opt::bases; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FILTER_ALIGNMENT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FILTER_ALIGNMENT_VERSION_MESSAGE;
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

    if(opt::positions_file.empty()) {
        std::cerr << SUBPROGRAM ": a --positions_file file must be provided\n";
        die = true;
    }
    if(opt::alignment_files.empty()) {
        std::cerr << SUBPROGRAM ": a --alignment_files file must be provided\n";
        die = true;
    }
    if(opt::output_dir.empty()) {
        std::cerr << SUBPROGRAM ": a --output_dir file must be provided\n";
        die = true;
    }
    if (die)
    {
        std::cout << "\n" << FILTER_ALIGNMENT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}


int filter_alignments_main(int argc, char** argv)
{
    parse_filter_main_options(argc, argv);

#ifndef H5_HAVE_THREADSAFE
    if(opt::threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    omp_set_num_threads(opt::threads); // Use 4 threads for all consecutive parallel regions
    filter_alignment_files(opt::alignment_files, opt::positions_file, opt::output_dir, opt::bases);

    return EXIT_SUCCESS;
}
