//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_POSITIONSFILE_HPP_
#define EMBED_FAST5_SRC_POSITIONSFILE_HPP_

// Boost
#include <boost/icl/interval_set.hpp>
#include <boost/filesystem.hpp>
#include <boost/coroutine2/all.hpp>


// std lib
#include <map>
#include <string>
#include <tuple>

using namespace std;
using namespace boost::icl;
using namespace boost::coroutines2;



/** Positions Line
 *
 * used to pass the
 * @param contig: contig
 * @param position: position on the contig
 * @param strand: strand of contig
 * @param change_from: reference nucleotide
 * @param change_to: expected nucleotide or nucleotides
 *
 */
struct PositionLine
{
  string contig;
  uint64_t position;
  string strand;
  string change_from;
  string change_to;
  PositionLine(string contig, uint64_t position, string strand, string change_from, string change_to) :
      contig(std::move(contig)), position(move(position)), strand(std::move(strand)),
      change_from(std::move(change_from)), change_to(std::move(change_to))
  {
  }
  PositionLine()
  = default;
  ~PositionLine()
  = default;
  bool operator==(const PositionLine &p) const{
    return (p.contig == contig) and (p.position == position) and (p.strand == strand) and
        (p.change_from == change_from) and (p.change_to == change_to);
  }
  bool operator!=(const PositionLine &other) const {
    return !(*this == other);
  }

};


typedef coroutine<PositionLine> positions_coro;


class PositionsFile
{
 public:
  PositionsFile();
  PositionsFile(const string& input_reads_filename, int64_t k);
  explicit PositionsFile(const string& input_reads_filename);
  ~PositionsFile();
  void positions_coroutine(positions_coro::push_type& yield);
  positions_coro::pull_type iterate();
  void load_interval_map(int64_t k);
  void load_positions_map();
  void load_positions_map(const std::string& input_reads_filename);
  bool is_in(string& contig, int64_t position);
  //
  string file_path;
  map<string, string> ambig_bases;

  map<string, interval_set<int64_t>> m_data;
  map<string, map<uint64_t, tuple<string, string>>> positions_map;

};


#endif //EMBED_FAST5_SRC_POSITIONSFILE_HPP_
