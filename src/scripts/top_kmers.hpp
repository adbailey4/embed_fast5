//
// Created by Andrew Bailey on 2019-06-17.
//

#ifndef EMBED_FAST5_TOP_KMERS_H
#define EMBED_FAST5_TOP_KMERS_H

#include <string>

using namespace std;

void generate_master_assignment_table(string assignment_dir, string& output_dir, int heap_size, string& alphabet);
int top_kmers_main(int argc, char** argv);


#endif //EMBED_FAST5_TOP_KMERS_H
