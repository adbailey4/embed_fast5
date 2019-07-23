//
// Created by Andrew Bailey on 2019-07-19.
//

#include "FullSaEvent.hpp"

FullSaEvent::FullSaEvent(string contig, uint64_t reference_index, string reference_kmer, string read_file, string read_strand,
                          uint64_t event_index, double event_mean, double event_noise, double event_duration, string aligned_kmer,
                          double scaled_mean_current, double scaled_noise, double posterior_probability, double descaled_event_mean,
                          double ont_model_mean, string path_kmer) :
  contig(std::move(contig)), reference_index(reference_index), reference_kmer(std::move(reference_kmer)), read_file(std::move(read_file)),
  read_strand(std::move(read_strand)), event_index(event_index), event_mean(event_mean), event_noise(event_noise),
  event_duration(event_duration), aligned_kmer(std::move(aligned_kmer)), scaled_mean_current(scaled_mean_current),
  scaled_noise(scaled_noise), posterior_probability(posterior_probability),
  descaled_event_mean(descaled_event_mean), ont_model_mean(ont_model_mean), path_kmer(std::move(path_kmer))
{
}

FullSaEvent::~FullSaEvent()
= default;