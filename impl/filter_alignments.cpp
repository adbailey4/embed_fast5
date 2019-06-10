//
// Created by Andrew Bailey on 2019-06-07.
//

#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <string>
#include "omp.h"
#include <boost/icl/discrete_interval.hpp>
#include <boost/filesystem.hpp>
#include "filter_alignments.hpp"


using namespace std;
using namespace boost;
using namespace boost::icl;
using namespace boost::filesystem;

PositionsFile::~PositionsFile()
= default;

PositionsFile::PositionsFile()
= default;


PositionsFile::PositionsFile(const std::string& input_reads_filename, int64_t k){
    this->load(input_reads_filename, k);
}


// Split a string into parts based on the delimiter
std::vector<std::string> split(std::string &in, char delimiter)
{
    std::vector<std::string> out;
    size_t lastPos = 0;
    size_t pos = in.find_first_of(delimiter);

    while(pos != std::string::npos)
    {
        out.push_back(in.substr(lastPos, pos - lastPos));
        lastPos = pos + 1;
        pos = in.find_first_of(delimiter, lastPos);
    }
    out.push_back(in.substr(lastPos));
    return out;
}

//https://www.systutorials.com/131/convert-string-to-int-and-reverse/
int64_t convert_to_int(std::string& str_int){
    int64_t number;
    std::istringstream iss (str_int);
    iss >> number;
    if (!iss.good ()) {
        return number;
    }
}


void PositionsFile::load(const std::string& input_reads_filename, int64_t k)
{
    // generate input filenames
    std::ifstream in_file(input_reads_filename.c_str());
    if(in_file.good()) {
        // read the file
        std::string line;
        while(getline(in_file, line)) {
            std::vector<std::string> fields = split(line, '\t');
            string contig = fields[0];
            string strand = fields[2];
            string contig_strand = contig+strand;
            int64_t start_position = convert_to_int(fields[1]);
            int64_t end_position = start_position - (k - 1);

            if ( m_data.find(contig_strand) == m_data.end() ) {
                // not found
                interval_set<int64_t> intervalSet;
                m_data[contig_strand] = intervalSet;
            }
//            create interval
            discrete_interval<int64_t>::type positions(end_position, start_position, interval_bounds::closed());
            m_data[contig_strand].insert(positions);
        }
    }
    in_file.close();
}

//check to see if position is in interval tree
bool PositionsFile::is_in(string& contig, int64_t position){
    if ( m_data.find(contig) == m_data.end() ) {
        // not found
        return false;
    }
    return contains(m_data[contig], position) ;
}

AlignmentFile::~AlignmentFile()
= default;

AlignmentFile::AlignmentFile(const string& input_reads_filename){
    this->file_path = input_reads_filename;
    this->get_strand();
    this->get_k();
}

int64_t AlignmentFile::get_k(){
    std::ifstream in_file(this->file_path.c_str());
    if (in_file.good()) {
        // read the file
        std::string line;
        getline(in_file, line);
        std::vector<std::string> fields = split(line, '\t');
        this->k = fields[2].length();
    }
    in_file.close();
    return this->k;
}

void AlignmentFile::filter(PositionsFile* pf, boost::filesystem::path& output_file) {
    std::ofstream out_file;
    out_file.open(output_file.string());
    std::ifstream in_file(this->file_path.c_str());

    if (in_file.good()) {
        // read the file
        std::string line;
        while (getline(in_file, line)) {
            std::vector<std::string> fields = split(line, '\t');
            string contig = fields[0];
            int64_t reference_index = convert_to_int(fields[1]);
//            string reference_kmer = fields[2];
//            string read_file = fields[3];
            string read_strand = fields[4];
//            string event_index = fields[5];
//            string event_mean = fields[6];
//            string event_noise = fields[7];
//            string event_duration = fields[8];
//            string aligned_kmer = fields[9];
//            string scaled_mean_current = fields[10];
//            string scaled_noise = fields[11];
            string posterior_probability = fields[12];
            string descaled_event_mean = fields[13];
//            string ont_model_mean = fields[14];
            string path_kmer = fields[15];
            string contig_strand = contig+this->strand;
            if (pf->is_in(contig_strand, reference_index))
                out_file << path_kmer << '\t' << read_strand << '\t' << descaled_event_mean << '\t' <<  posterior_probability << '\n';
        }
    }
    out_file.close();
    in_file.close();
}

string AlignmentFile::get_strand(){
    std::vector<std::string> fields = split(this->file_path, '.');
    if (fields.end()[-2] == "backward"){
        this->strand = "-";
    } else if (fields.end()[-2] == "forward") {
        this->strand = "+";
    } else {
        fprintf(stderr, "error: could not infer strand from  %s\n", this->file_path.c_str());
        fprintf(stderr, "Please check input file is full alignment file from signalalign\n");
        exit(EXIT_FAILURE);

    }
    return this->strand;
}

path make_dir(path output_path){
    if (!exists(output_path)){
        create_directory(output_path);
    }
    return output_path;
}

void filter_alignment_files(string input_reads_dir, const string& positions_file, string output_dir){

    path p(input_reads_dir);
    path output_path = make_dir(output_dir);

    directory_iterator end_itr;
//    Get all tsvs to process
    vector<path> all_tsvs;
    int64_t k;
    int counter = 0;
    for (directory_iterator itr(p); itr != end_itr; ++itr) {
        if (is_regular_file(itr->path()) and itr->path().extension().string() == ".tsv") {
            all_tsvs.push_back(itr->path());
        }
        if (counter == 0){
            AlignmentFile af = AlignmentFile(itr->path().string());
            k = af.k;
            counter += 1;
        }
    }

    PositionsFile pf = PositionsFile(positions_file, k);

    int64_t number_of_files = all_tsvs.size();
    path* array_of_files = &all_tsvs[0];
// looping through the files
    #pragma omp parallel for shared(array_of_files, pf)
    for(int64_t i=0; i < number_of_files; i++) {

        path current_file = array_of_files[i];
        AlignmentFile af = AlignmentFile(current_file.string());

        path output_file = output_path / current_file.filename();
        af.filter(&pf, output_file);
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
        "  -p, --positions_file=FILE            path to positons file\n"
        "  -o, --output_dir=DIR                 path to directory to output filteredalignment files\n"
        "  -t, --threads=NUMBER                 number of threads\n"

        "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string alignment_files;
    static std::string positions_file;
    static std::string output_dir;
    static unsigned int threads;
    static int num_threads = 1;
}

static const char* shortopts = "r:i:t:f:vn";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
        { "verbose",          no_argument,       nullptr, 'v' },
        { "alignment_files",  required_argument, nullptr, 'a' },
        { "positions_file",   required_argument, nullptr, 'p' },
        { "output_dir",       required_argument, nullptr, 'o' },
        { "threads",          optional_argument, nullptr, 't' },
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
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FILTER_ALIGNMENT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FILTER_ALIGNMENT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
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
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    omp_set_num_threads(opt::threads); // Use 4 threads for all consecutive parallel regions
    filter_alignment_files(opt::alignment_files, opt::positions_file, opt::output_dir);

    return EXIT_SUCCESS;
}
