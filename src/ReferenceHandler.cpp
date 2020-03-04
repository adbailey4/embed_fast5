//
// Created by Kishwar Shafin on 6/14/18.
// Edited by Andrew Bailey on 2019-07-25.
//

#include "ReferenceHandler.hpp"

ReferenceHandler::ReferenceHandler(const string& path) {
  this->fasta = fai_load(path.c_str());
  if (fasta == NULL) {
    cerr<<"INVALID FASTA FILE. PLEASE CHECK IF PATH IS CORRECT AND FILE IS INDEXED: "<<path<<endl;
    exit (EXIT_FAILURE);
  }
}

ReferenceHandler::~ReferenceHandler() {
  fai_destroy(this->fasta);
}

vector<string> ReferenceHandler::get_chromosome_names() {
  vector<string> chromosome_names;
  int number_of_chromosome_sequences = faidx_nseq(this->fasta);

  for(int i=0; i < number_of_chromosome_sequences; i++) {
    string chromosome_name = faidx_iseq(this->fasta, i);
    chromosome_names.push_back(chromosome_name);
  }

  return chromosome_names;
}

string ReferenceHandler::get_reference_sequence(const string& region, long long start, long long stop) {
  // the fetch is zero-based, in one-based system the previous base would be the start.
  // this is done to make the reference sequence compatible with the BAM file's read sequences.
  // start += 1;
  // stop += 1;

  int len = 0;
  string sequence;

  sequence = faidx_fetch_seq(fasta, region.c_str(), start, stop - 1, &len);
  //-2 if c_name not present, -1 general error
  if(len == -2){
    cerr<<"CHROMOSOME NAME NOT PRESENT IN REFERENCE FASTA FILE: "<<region<<" "<<start<<" "<<stop<<endl;
    return nullptr;
  } else if(len == -1){
    cerr<<"ENCOUNTERED ERROR IN FETCHING REFERENCE FASTA FILE: "<<region<<" "<<start<<" "<<stop<<endl;
    return nullptr;
  }
  return sequence;
}

uint64_t ReferenceHandler::get_chromosome_sequence_length(const string& chromosome_name) {
  uint64_t len = faidx_seq_len(this->fasta, chromosome_name.c_str());
  return len;
}
