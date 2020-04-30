//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_
#define EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_

#include "PositionsFile.hpp"
#include "VariantCall.hpp"
#include <boost/filesystem.hpp>
#include <boost/coroutine2/all.hpp>
#include <utility>
#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;
using namespace boost::filesystem;
using namespace boost::coroutines2;


class FullSaEvent {
 public:
  string contig;
  uint64_t reference_index;
  string reference_kmer;
  string read_file;
  string strand;
  uint64_t event_index;
  double event_mean;
  double event_noise;
  double event_duration;
  string aligned_kmer;
  double scaled_mean_current;
  double scaled_noise;
  double posterior_probability;
  double descaled_event_mean;
  double ont_model_mean;
  string path_kmer;
  FullSaEvent(string contig, uint64_t reference_index, string reference_kmer, string read_file, string read_strand,
              uint64_t event_index, double event_mean, double event_noise, double event_duration, string aligned_kmer,
              double scaled_mean_current, double scaled_noise, double posterior_probability, double descaled_event_mean,
              double ont_model_mean, string path_kmer) :
      contig(std::move(contig)), reference_index(reference_index), reference_kmer(std::move(reference_kmer)), read_file(std::move(read_file)),
      strand(std::move(read_strand)), event_index(event_index), event_mean(event_mean), event_noise(event_noise),
      event_duration(event_duration), aligned_kmer(std::move(aligned_kmer)), scaled_mean_current(scaled_mean_current),
      scaled_noise(scaled_noise), posterior_probability(posterior_probability),
      descaled_event_mean(descaled_event_mean), ont_model_mean(ont_model_mean), path_kmer(std::move(path_kmer))
  {
  }

  ~FullSaEvent() = default;
  string format_line(bool write_full=true) const{
    ostringstream person_info;
    if (write_full) {
      person_info << contig << '\t' << reference_index << '\t' << reference_kmer << '\t' << read_file << '\t'
                  << strand << '\t' << event_index << '\t' << event_mean << '\t' << event_noise << '\t'
                  << event_duration << '\t' << aligned_kmer << '\t' << scaled_mean_current << '\t'
                  << scaled_noise << '\t' << posterior_probability << '\t' << descaled_event_mean << '\t'
                  << ont_model_mean << '\t' << path_kmer << '\n';
    } else {
      person_info << path_kmer << '\t' << strand << '\t' << descaled_event_mean << '\t'
                  << posterior_probability << '\n';
    }
    return person_info.str();
  };

};

typedef coroutine<FullSaEvent> full_sa_coro;


class AlignmentFile
{
 public:
  explicit AlignmentFile(string input_reads_filename);
  AlignmentFile(string input_reads_filename, bool is_rna);
  ~AlignmentFile();
  string get_strand();
  int64_t get_k();
  void filter_by_positions(PositionsFile *pf, path &output_file, string bases);
  full_sa_coro::pull_type iterate();
  full_sa_coro::pull_type filter_by_ref_bases(string& bases);
  vector<VariantCall> get_variant_calls(string& ambig_bases, std::map<string, string> *ambig_bases_map);
    //
  string file_path;
  bool good_file;
  bool rna = false;
  string strand;
  string contig;
  string read_id;
  int64_t k;
  AlignmentFile(const AlignmentFile&) = delete;
  AlignmentFile& operator=(const AlignmentFile&) = delete;

 private:
  std::ifstream in_file;
  void push_iterate(full_sa_coro::push_type& yield);
  void push_filter_by_ref_bases(full_sa_coro::push_type& yield, string bases);

};

#endif //EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_
