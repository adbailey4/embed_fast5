//
// Created by Andrew Bailey on 2019-06-07.
//

#ifndef EMBED_FAST5_FILTER_ALIGNMENTS_H
#define EMBED_FAST5_FILTER_ALIGNMENTS_H


#include "PositionsFile.hpp"

using namespace std;

int filter_alignments_main(int argc, char** argv);
void filter_alignment_files(string &input_reads_dir, const string& positions_file, string &output_dir, string& bases);

#endif //EMBED_FAST5_FILTER_ALIGNMENTS_H
