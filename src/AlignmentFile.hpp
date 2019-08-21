#include <utility>

//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_
#define EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_

#include "PositionsFile.hpp"
#include "VariantCall.hpp"
#include "FullSaEvent.hpp"
#include <boost/filesystem.hpp>
#include <boost/coroutine2/all.hpp>
#include <utility>
#include <fstream>
#include <iostream>


using namespace std;
using namespace boost::filesystem;
using namespace boost::coroutines2;


class AlignmentFile
{
 public:
  explicit AlignmentFile(string input_reads_filename);
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
  string strand;
  int64_t k;
  AlignmentFile(const AlignmentFile&) = delete;
  AlignmentFile& operator=(const AlignmentFile&) = delete;

 private:
  std::ifstream in_file;
  void push_iterate(full_sa_coro::push_type& yield);
  void push_filter_by_ref_bases(full_sa_coro::push_type& yield, string bases);

};

#endif //EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_
