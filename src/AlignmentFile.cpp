#include <utility>

//
// Created by Andrew Bailey on 2019-07-15.
//

#include "AlignmentFile.hpp"
#include "EmbedUtils.hpp"
#include <numeric>
#include <iostream>
#include <fstream>
#include <stdexcept>


using namespace embed_utils;
using namespace std;


/**
 * Deconstructor for AlignmentFile. Need to close file
 */
AlignmentFile::~AlignmentFile(){
  this->in_file.close();
}

/**
 * Constructor for AlignmentFile.
 * Sets file_path, opens file, reports error to screen if file is not good then gets strand and kmer len
 */
AlignmentFile::AlignmentFile(string input_reads_filename) :
    file_path(std::move(input_reads_filename)){
  this->in_file.open(file_path);
  this->good_file = in_file.good();
  if (!this->good_file){
    cout << "Error loading file: " << this->file_path << "\n";
  } else{
    strand = this->get_strand();
    k = this->get_k();
  }
}

/**
 * Constructor2 for AlignmentFile.
 * Sets file_path, opens file, reports error to screen if file is not good then gets strand and kmer len
 */
AlignmentFile::AlignmentFile(string input_reads_filename, bool is_rna) :
    file_path(std::move(input_reads_filename)){
  this->in_file.open(file_path);
  this->good_file = in_file.good();
  if (!this->good_file){
    cout << "Error loading file: " << this->file_path << "\n";
  } else{
    strand = this->get_strand();
    k = this->get_k();
  }
  rna = is_rna;
}


/**
 * Get the strand of the read based on the file naming
*/
string AlignmentFile::get_strand(){
  std::vector<std::string> fields = split_string(this->file_path, '.');
  string this_strand;
  if (fields.end()[-2] == "backward"){
    this_strand = "-";
  } else if (fields.end()[-2] == "forward") {
    this_strand = "+";
  } else {
    fprintf(stderr, "error: could not infer strand from  %s\n", this->file_path.c_str());
    fprintf(stderr, "Please check input file is full alignment file from signalalign\n");
    exit(EXIT_FAILURE);
  }
  return this_strand;
}

/**
 * Reads in a line of the file to get kmer length from the reference kmer field.
 *
 * @return: -1 if file is not good otherwise the kmer length
 */
int64_t AlignmentFile::get_k(){
  int64_t kmer_len = -1;
  if (this->good_file) {
    in_file.clear();
    in_file.seekg(0, ios::beg);
    // read the file
    std::string line;
    getline(in_file, line);
    std::vector<std::string> fields = split_string(line, '\t');
    kmer_len = fields[2].length();
    read_id = fields[3];
    contig = fields[0];
  }
  return kmer_len;
}

/**
 * Filters out events which do not cover positions in the positions file and writes to a new file.
 *
 * @param pf: PositionsFile
 * @param output_file: path to filtered output file
 * @param bases: characters to exclude if kmers contain the character unless the position is in the PositionsFile
 */
void AlignmentFile::filter_by_positions(PositionsFile *pf, boost::filesystem::path &output_file, string bases) {
  std::ofstream out_file;
  for (auto &event: this->iterate()) {
      string contig_strand = event.contig+this->strand;
      if (pf->is_in(contig_strand, event.reference_index) || !are_characters_in_string(bases, event.path_kmer)) {
        out_file << event.path_kmer << '\t' << event.strand << '\t' << event.descaled_event_mean << '\t' << event.posterior_probability << '\n';
    }
  }
  out_file.close();
}


/**
 * Create a push type coroutine for parsing an alignment file
*/
void AlignmentFile::push_iterate(full_sa_coro::push_type& yield){
  if(this->good_file) {
    in_file.clear();
    in_file.seekg(0, ios::beg);
    string line;
    while(getline(this->in_file, line)) {
      vector<string> fields = split_string(line, '\t');
      if (fields.size() != 16) {
        continue;
      } else {
        FullSaEvent sa(fields[0],
                       string_to_int(fields[1]), fields[2], fields[3], fields[4],
                       string_to_int(fields[5]), string_to_float(fields[6]),
                       string_to_float(fields[7]), string_to_float(fields[8]), fields[9], string_to_float(fields[10]),
                       string_to_float(fields[11]),
                       string_to_float(fields[12]),
                       string_to_float(fields[13]),
                       string_to_float(fields[14]), fields[15]);
        yield(sa);
      }
    }
  }
}


/**
 * Iterate over all rows of the file
 * Call push type coroutine to create a pull type coroutine
*/
full_sa_coro::pull_type AlignmentFile::iterate() {
  full_sa_coro::pull_type event{bind(&AlignmentFile::push_iterate, this, std::placeholders::_1)};
  return event;
}


/**
 * Filters out events which do not have specific characters within the reference kmer.
 *
 * @param yield: full_sa_coro push type
 * @param bases: characters to search for within reference kmer
 * @return: -1 if file is not good otherwise the kmer length
 */
void AlignmentFile::push_filter_by_ref_bases(full_sa_coro::push_type& yield, string bases) {
  for (auto &event: this->iterate()) {
    if (are_characters_in_string(bases, event.reference_kmer)){
      yield(event);
    }
  }
}

/**
 * Iterate over all rows of the file
 * Call push type coroutine to create a pull type coroutine
 * @param bases: characters to search for within reference kmer
*/
full_sa_coro::pull_type AlignmentFile::filter_by_ref_bases(string& bases) {
  full_sa_coro::pull_type event{bind(&AlignmentFile::push_filter_by_ref_bases, this, std::placeholders::_1, bases)};
  return event;
}

/**
 * Get variant calls from full alignment file
 *
 * @param ambig_bases: the list of possible ambiguous calls
 * @param ambig_bases_map: map of ambig bases to represented nucleotides
 * @return: vector of variant calls
*/
vector<VariantCall> AlignmentFile::get_variant_calls(string &ambig_bases, std::map<string, string> *ambig_bases_map) {

  std::map<uint64_t, VariantCall> variant_calls;
  std::map<uint64_t, VariantCall>::iterator it;
  std::size_t location;
  uint64_t position;
  uint64_t path_kmer_pos;
  uint64_t index;
  string possible_bases;
  VariantCall position_call;
  vector<VariantCall> final_output;
  string s;

//  loop through file
  try {
    for (auto &event: this->iterate()) {
      for (char &c : ambig_bases) {
        path_kmer_pos = event.aligned_kmer.find(c);
        if (path_kmer_pos != std::string::npos) {
//        get position of ambiguous base
          if (rna){
            position = event.reference_index - path_kmer_pos + (k-1);
          } else {
            if (this->strand == "+") {
              position = event.reference_index + path_kmer_pos;
            } else {
              position = event.reference_index + (this->k - path_kmer_pos - 1);
            }
          }
          s = string(1, c);
          possible_bases = (*ambig_bases_map).at(s);
//        positon_call = VariantCalls[position];
          it = variant_calls.find(position);
          if (it == variant_calls.end()) {
//        initialize VariantCall if empty
//          position_call = VariantCall();
            position_call = VariantCall(event.contig, this->strand, position, possible_bases);
            for (int i = 0; i < possible_bases.length(); ++i) {
              position_call.normalized_probs.push_back(0.0);
              position_call.positional_probs.push_back(0.0);
              position_call.positional_probs2.push_back(0.0);
            }
          } else {
            position_call = variant_calls[position];
          }
//        get corresponding index for base call
          index = possible_bases.find(event.path_kmer[path_kmer_pos]);
          if (index != std::string::npos) {
            if (location == 0) {
              position_call.positional_probs[index] += event.posterior_probability;
            } else {
              position_call.positional_probs2[index] += event.posterior_probability;
            }
            variant_calls[position] = position_call;
          } else {
            throw runtime_error("Programmer Error: This should never happen yo.");
          }
        }
      }
    }

    // Iterate over the map using c++11 range based for loop
    for (std::pair<uint64_t, VariantCall> element : variant_calls) {
      double sum_of_elements = std::accumulate(element.second.positional_probs.begin(),
                                               element.second.positional_probs.end(), 0.0);

      if (sum_of_elements > 0.0) {
        for (size_t i = 0; i < element.second.positional_probs.size(); ++i)
          element.second.normalized_probs[i] =
              element.second.positional_probs[i] / sum_of_elements;
      } else {
        sum_of_elements = std::accumulate(element.second.positional_probs2.begin(),
                                                 element.second.positional_probs2.end(), 0.0);
        for (size_t i = 0; i < element.second.positional_probs2.size(); ++i)
          element.second.normalized_probs[i] =
              element.second.positional_probs2[i] / sum_of_elements;
      }
      final_output.push_back(element.second);
    }
  } catch (const std::out_of_range& oor) {
    throw runtime_error("Programmer Error: ambig_bases not in ambig_bases_map. base: " + s);
  }

  return final_output;
}

