//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_POSITIONSFILE_HPP_
#define EMBED_FAST5_SRC_POSITIONSFILE_HPP_

#include <boost/icl/interval_set.hpp>
#include <boost/filesystem.hpp>
#include <map>

using namespace std;
using namespace boost::icl;

class PositionsFile
{
 public:
  PositionsFile();
  PositionsFile(const string& input_reads_filename, int64_t k);
  ~PositionsFile();

  void load(const string& input_reads_filename, int64_t k);
  bool is_in(string& contig, int64_t position);
  //
  map<string, interval_set<int64_t>> m_data;

};


#endif //EMBED_FAST5_SRC_POSITIONSFILE_HPP_
