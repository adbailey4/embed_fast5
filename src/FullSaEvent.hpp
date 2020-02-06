//
// Created by Andrew Bailey on 2019-07-19.
//

#ifndef EMBED_FAST5_SRC_FULLSAEVENT_HPP_
#define EMBED_FAST5_SRC_FULLSAEVENT_HPP_


#include <boost/coroutine2/all.hpp>

using namespace std;
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
  ~FullSaEvent();
  FullSaEvent(string contig, uint64_t reference_index, string reference_kmer, string read_file, string strand,
              uint64_t event_index, double event_mean, double event_noise, double event_duration, string aligned_kmer,
              double scaled_mean_current, double scaled_noise, double posterior_probability, double descaled_event_mean,
              double ont_model_mean, string path_kmer);
};

typedef coroutine<FullSaEvent> full_sa_coro;



#endif //EMBED_FAST5_SRC_FULLSAEVENT_HPP_
