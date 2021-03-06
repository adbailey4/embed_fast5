//
// Created by Andrew Bailey on 03/15/19.
//


#include <cassert>
#include <iostream>
#include <string>
#include <fast5.hpp>
#include <boost/filesystem.hpp>
#include "EmbedFast5.hpp"
#include "TopKmers.hpp"
#include "FilterAlignments.hpp"
#include "SignalAlignToBed.hpp"
#include "SplitByRefPosition.hpp"
#include "PositionsKmerDistributions.hpp"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_index.h"
#include "nanopolish_extract.h"
#include "nanopolish_eventalign.h"

using namespace std;
using namespace boost::filesystem;



int print_usage(int argc, char **argv);

static std::map< std::string, std::function<int(int, char**)> > programs = {
    {"help",                     print_usage},
    {"--help",                   print_usage},
    {"embed",                    embed_fast5_main},
    {"index",                    index_main},
    {"extract",                  extract_main},
    {"eventalign",               eventalign_main},
    {"filter_by_positions",      filter_alignments_main},
    {"top_kmers",                top_kmers_main},
    {"split_by_position",        split_by_ref_main},
    {"kmer_distributions",       get_kmer_distributions_main},
    {"sa2bed",                   sa2bed_main}
};

int print_usage(int, char **)
{
    std::cout << "usage: embed [options]" << std::endl;
    std::cout << "  valid commands: " << std::endl;
    for (const auto &item : programs){
        std::cout << "    " << item.first << std::endl;
    }
    std::cout << "  for help on given command, type main command --help" << std::endl;
    return 0;
}


int main(int argc, char** argv) {
    // Turn off HDF's exception printing, which is generally unhelpful for users
//    H5Eset_auto(0, NULL, NULL);

    int ret = 0;
    if(argc <= 1) {
        printf("error: no command provided\n");
        print_usage(argc - 1 , argv + 1);
        return 0;
    } else {
        std::string command(argv[1]);
        auto iter = programs.find(command);
        if (iter != programs.end())
            ret = iter->second( argc - 1, argv + 1);
        else
            ret = print_usage( argc - 1, argv + 1);
    }
}
