//
// Created by Andrew Bailey on 6/2/20.
//

#ifndef EMBED_FAST5_SRC_POSITIONSKMERDISTRIBUTIONS_HPP_
#define EMBED_FAST5_SRC_POSITIONSKMERDISTRIBUTIONS_HPP_

#include <string>

using namespace std;

int get_kmer_distributions_main(int argc, char** argv);

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
                                        const float& min_prob_threshold);

#endif //EMBED_FAST5_SRC_POSITIONSKMERDISTRIBUTIONS_HPP_
