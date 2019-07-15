//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_
#define EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_

#include "PositionsFile.hpp"
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::filesystem;

class AlignmentFile
{
 public:
  explicit AlignmentFile(const string& input_reads_filename);
  ~AlignmentFile();
  string get_strand();
  int64_t get_k();
  void filter(PositionsFile* pf, path& output_file, string bases);
  //
  string strand;
  string file_path;
  int64_t k;
};

#endif //EMBED_FAST5_SRC_ALIGNMENTFILE_HPP_
