#include <climits>
//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_ASSIGNMENTFILE_HPP_
#define EMBED_FAST5_SRC_ASSIGNMENTFILE_HPP_

#include <boost/coroutine2/all.hpp>
#include <string>
#include <sstream>


using namespace std;
using namespace boost::coroutines2;

struct eventkmer
{
  string path_kmer;
  float descaled_event_mean;
  string strand;
  float posterior_probability;
  eventkmer(string kmer, float mean, string strand, float prob) :
      path_kmer(move(kmer)), descaled_event_mean(mean), strand(move(strand)), posterior_probability(prob)
  {
  }
  ~eventkmer() = default;
  string format_line(__unused bool trim=false) const{
    ostringstream person_info;
    person_info << path_kmer << '\t' << strand << '\t' << descaled_event_mean << '\t' << posterior_probability << '\n';
    return person_info.str();
  };

};


typedef coroutine<eventkmer> event_kmer_coro;



class AssignmentFile
{
 public:
  explicit AssignmentFile(const string& input_reads_filename);
  ~AssignmentFile();

  void assignment_coroutine(event_kmer_coro::push_type& yield);
  event_kmer_coro::pull_type iterate();
  int64_t get_k();

  string file_path;
  int64_t k;

};

#endif //EMBED_FAST5_SRC_ASSIGNMENTFILE_HPP_
