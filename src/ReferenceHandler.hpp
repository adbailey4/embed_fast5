//
// Created by Kishwar Shafin on 6/14/18.
// Edited by Andrew Bailey on 2019-07-25.
//

#ifndef EMBED_FAST5_SRC_REFERENCEHANDLER_HPP_
#define EMBED_FAST5_SRC_REFERENCEHANDLER_HPP_

#include "htslib/faidx.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;


class ReferenceHandler {
 public:
  ReferenceHandler(const string& path);
  // this start and stop is zero-based. Should we change everything to one-based?
  string get_reference_sequence(const string& region, long long start, long long stop);
  uint64_t get_chromosome_sequence_length(const string& chromosome_name);
  vector<string> get_chromosome_names();
  ~ReferenceHandler();
 private:
  faidx_t* fasta;

};


#endif //EMBED_FAST5_SRC_REFERENCEHANDLER_HPP_
