//
// Created by Andrew Bailey on 2019-07-15.
//

#ifndef EMBED_FAST5_SRC_ASSIGNMENTFILE_HPP_
#define EMBED_FAST5_SRC_ASSIGNMENTFILE_HPP_

#include <boost/coroutine2/all.hpp>
#include <string>


using namespace std;
using namespace boost::coroutines2;

struct eventkmer
{
  string kmer;
  float mean;
  string strand;
  float prob;
  eventkmer(string kmer, float mean, string strand, float prob) :
      kmer(move(kmer)), mean(mean), strand(move(strand)), prob(prob)
  {
  }
  ~eventkmer()
  = default;
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
