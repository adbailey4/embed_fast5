//
// Created by Andrew Bailey on 2019-07-24.
//

#ifndef EMBED_FAST5_SRC_SCRIPTS_SPLITBYREFPOSITION_HPP_
#define EMBED_FAST5_SRC_SCRIPTS_SPLITBYREFPOSITION_HPP_

#include <string>

using namespace std;

int split_by_ref_main(int argc, char** argv);
void split_signal_align_by_ref_position(
    string& sa_input_dir,
    string& output_file_path,
    string reference,
    uint64_t num_locks,
    uint64_t n_threads,
    bool verbose,
    bool rna,
    bool two_d);

#endif //EMBED_FAST5_SRC_SCRIPTS_SPLITBYREFPOSITION_HPP_
