//
// Created by Andrew Bailey on 2019-07-19.
//

#ifndef EMBED_FAST5_SRC_VARIANTCALL_HPP_
#define EMBED_FAST5_SRC_VARIANTCALL_HPP_

#include <string>
#include <vector>

using namespace std;

class VariantCall {
  public:
    explicit VariantCall(string contig1, string strand1, uint64_t reference_index1, string bases1);
    VariantCall();
    ~VariantCall();
    string contig;
    string strand;
    uint64_t reference_index;
    string bases;
    vector<double> positional_probs;
    vector<double> positional_probs2;
    vector<double> normalized_probs;
};

#endif //EMBED_FAST5_SRC_VARIANTCALL_HPP_
