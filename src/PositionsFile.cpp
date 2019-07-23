//
// Created by Andrew Bailey on 2019-07-15.
//

#include "PositionsFile.hpp"
#include "EmbedUtils.hpp"

using namespace std;
using namespace embed_utils;


PositionsFile::~PositionsFile()
= default;

PositionsFile::PositionsFile()
= default;


PositionsFile::PositionsFile(const std::string& input_reads_filename, int64_t k){
  throw_assert(k!=0, "k is zero. Check initialization of PositionsFile")
  this->load(input_reads_filename, k);
}


void PositionsFile::load(const std::string& input_reads_filename, int64_t k)
{
  // generate input filenames
  std::ifstream in_file(input_reads_filename.c_str());
  if(in_file.good()) {
    // read the file
    std::string line;
    while(getline(in_file, line)) {
      std::vector<std::string> fields = embed_utils::split_string(line, '\t');
      string contig = fields[0];
      string strand = fields[2];
      string contig_strand = contig+strand;
      int64_t start_position = embed_utils::convert_to_int(fields[1]);
      int64_t end_position = start_position - (k - 1);

      if ( m_data.find(contig_strand) == m_data.end() ) {
        // not found
        interval_set<int64_t> intervalSet;
        m_data[contig_strand] = intervalSet;
      }
//            create interval
      discrete_interval<int64_t>::type positions(end_position, start_position, interval_bounds::closed());
      m_data[contig_strand].insert(positions);
    }
  }
  in_file.close();
}

//check to see if position is in interval tree
bool PositionsFile::is_in(string& contig, int64_t position){
  const auto& iter = m_data.find(contig);
  if(iter == m_data.end()) {
    return false;
  } else {
    return contains(m_data[contig], position) ;
  }
}