//
// Created by Andrew Bailey on 2019-07-15.
//


// Embed libs
#include "PositionsFile.hpp"
#include "EmbedUtils.hpp"

using namespace std;
using namespace embed_utils;

/**
 * Deconstructor for PositionsFile.
 */
PositionsFile::~PositionsFile()
= default;

/**
 * Default constructor for PositionsFile.
 */
PositionsFile::PositionsFile(){
  this->ambig_bases = create_ambig_bases();
}
/**
 * Constructor which only initializes file path and does not "load" data.
 * @param input_reads_filename: path to positions file
 */
PositionsFile::PositionsFile(const std::string& input_reads_filename)
{
  this->file_path = input_reads_filename;
  this->ambig_bases = create_ambig_bases();
}
/**
 * Constructor for PositionsFile.
 * loads data and throws error if k is zero
 * @param input_reads_filename: path to positions file
 * @param k: kmer size
 */
PositionsFile::PositionsFile(const std::string& input_reads_filename, int64_t k){
  throw_assert(k!=0, "k is zero. Check initialization of PositionsFile")
  this->file_path = input_reads_filename;
  this->ambig_bases = create_ambig_bases();
  this->load_interval_map(k);
}

/**
 * Load positions file into interval map
  @param input_reads_filename: path to positions file
  @param k: kmer length
 */
void PositionsFile::load_interval_map(int64_t k)
{
  for (auto &position_line: this->iterate()){
    int64_t end_position = position_line.position - (k - 1);
    string contig_strand = position_line.contig+position_line.strand;

    if ( m_data.find(contig_strand) == m_data.end() ) {
      // not found
      interval_set<int64_t> intervalSet;
      m_data[contig_strand] = intervalSet;
    }
//            create interval
    discrete_interval<int64_t>::type positions(end_position, position_line.position, interval_bounds::closed());
    m_data[contig_strand].insert(positions);

  }
}

/**
 * Load positions into map
 *
 * [contig+strand][position] = {change_from, change_to}
 *
  @param input_reads_filename: path to positions file
  @param k: kmer length
 */
void PositionsFile::load_positions_map(const std::string& input_reads_filename){
  this->file_path = input_reads_filename;
  this->load_positions_map();
}


/**
 * Load positions into map
 *
 * [contig+strand][position] = {change_from, change_to}
 *
  @param k: kmer length
 */
void PositionsFile::load_positions_map()
{
  for (auto &position_line: this->iterate()){
    string contig_strand = position_line.contig+position_line.strand;

    if ( positions_map.find(contig_strand) == positions_map.end() ) {
      // create contig strand map and position map
      map<uint64_t, tuple<string, string>> contig_positions = {{position_line.position, make_tuple(position_line.change_from, position_line.change_to)}};
      positions_map[contig_strand] = contig_positions;
    } else {
//      add position to positions map
      positions_map[contig_strand][position_line.position] = make_tuple(position_line.change_from, position_line.change_to);
    }
  }
}

/**
 * check to see if position is in interval tree
 */
bool PositionsFile::is_in(string& contig, int64_t position){
  const auto& iter = m_data.find(contig);
  if(iter == m_data.end()) {
    return false;
  } else {
    return contains(m_data[contig], position) ;
  }
}

/**
Generate push type positions coroutine object
 file format:
contig  position    strand  change_from change_to
ecoli     20         +         A          T

*/
void PositionsFile::positions_coroutine(positions_coro::push_type& yield){
  std::ifstream in_file(this->file_path.c_str());
  if(in_file.good()) {
    // read the file
    std::string line;
    while(getline(in_file, line)) {
      std::vector<std::string> fields = embed_utils::split_string(line, '\t');
      string contig = fields[0];
      int64_t start_position = embed_utils::convert_to_int(fields[1]);
      string strand = fields[2];
      string change_from = fields[3];
      string change_to = fields[4];
      yield(PositionLine(contig, start_position, strand, change_from, change_to));
    }
  } else  {
    cout << "Error loading file: " << this->file_path << "\n";
  }
  in_file.close();
}

/**
Call push type coroutine to create a pull type coroutine
*/
positions_coro::pull_type PositionsFile::iterate(){
  positions_coro::pull_type seq {bind(&PositionsFile::positions_coroutine, this, std::placeholders::_1)};
  return seq;
}
