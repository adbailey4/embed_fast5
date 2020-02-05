//
// Created by Andrew Bailey on 2019-06-17.
//

#ifndef EMBED_FAST5_TOP_KMERS_H
#define EMBED_FAST5_TOP_KMERS_H

#include <string>

using namespace std;

string generate_master_assignment_table(string assignment_dir, string& output_dir, int heap_size, string& alphabet, unsigned int n_threads);
auto top_kmers_main(int argc, char** argv) -> int;


#endif //EMBED_FAST5_TOP_KMERS_H
