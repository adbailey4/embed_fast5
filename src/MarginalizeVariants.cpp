//
// Created by Andrew Bailey on 2019-07-15.
//

#include "MarginalizeVariants.hpp"
#include <iostream>
#include <fstream>



/**
 * Deconstructor for MarginalizeVariants.
 */
MarginalizeVariants::~MarginalizeVariants()
=default;

/**
 * Constructor for MarginalizeVariants.
 */
MarginalizeVariants::MarginalizeVariants() :
    num_locks(0), initialized_locks(false)
{
}

/**
 * Constructor for MarginalizeVariants.
 */
MarginalizeVariants::MarginalizeVariants(uint64_t num_locks) :
num_locks(num_locks), initialized_locks(true) {
  this->initialize_locks();

}

/**
Initialize vector of locks
*/
void MarginalizeVariants::initialize_locks() {
  locks = std::vector<std::mutex>(this->num_locks);
}

/**
Initialize vector of locks
*/
void MarginalizeVariants::initialize_locks(uint64_t number_of_locks) {
  this->num_locks = number_of_locks;
  locks = std::vector<std::mutex>(this->num_locks);
}


/**
 * Load a vector of variants into the per_genomic position data structure.
 *
 * @param vector_of_calls: vector of variant_calls
 */
void MarginalizeVariants::load_variants(vector<VariantCall>* vector_of_calls){
  for (auto& call: *vector_of_calls){
    if (this->initialized_locks){
      std::unique_lock<std::mutex> lk(this->locks[call.reference_index % this->num_locks]);
//      omp_set_lock(&(this->locks[i.reference_index % this->num_locks]));
      MarginalizeVariants::load_variant(call);
      lk.unlock();
//      omp_unset_lock(&(this->locks[i.reference_index % this->num_locks]));
    } else {
      MarginalizeVariants::load_variant(call);
    }
  }
}

void MarginalizeVariants::load_variant(VariantCall& call){
  bed_line bed_entry;
  int max_index;
  if (call.strand == "+"){
    bed_entry = this->per_genomic_position[call.contig].first[call.reference_index];
  } else {
    bed_entry = this->per_genomic_position[call.contig].second[call.reference_index];
  }

  if (bed_entry.start == -1) {
    bed_entry.start = call.reference_index;
    bed_entry.stop = call.reference_index + 1;
    bed_entry.bases = call.bases;
    size_t n_entries = call.normalized_probs.size();
    bed_entry.hits.resize(call.normalized_probs.size());
    for (size_t j = 0; j < n_entries; ++j) {
      bed_entry.hits.push_back(0.0);
    }
  }
  bed_entry.coverage += 1;
  max_index = max_element(call.normalized_probs.begin(), call.normalized_probs.end()) - call.normalized_probs.begin();
  bed_entry.hits[max_index] += 1;
  if (call.strand == "+"){
    this->per_genomic_position[call.contig].first[call.reference_index] = bed_entry;
  } else {
    this->per_genomic_position[call.contig].second[call.reference_index] = bed_entry;
  }
}

/**
 * Write per_genomic_position to file
 *
 * chrom chromStart chromEnd strand coverage characters hits
 * @param path_to_bed: path to bed file
 */
void MarginalizeVariants::write_to_file(path& path_to_bed) {
  boost::filesystem::ofstream myfile(path_to_bed);
  if (myfile.is_open())
  {
//    go through each contig
    for (pair<std::string, pair<map<uint64_t, bed_line>, map<uint64_t, bed_line>>> contig: this->per_genomic_position){
//      first of pair is plus strand
      for (pair<uint64_t, bed_line> strand: contig.second.first){
        myfile << contig.first << "\t" << strand.second.start << "\t" << strand.second.stop << "\t" << "+" << "\t" << strand.second.coverage << "\t" << strand.second.bases;
        for (uint64_t hit: strand.second.hits) {
          myfile << "\t";
          myfile << hit;
        }
        myfile << "\n";
      }
//      second of pair is minus strand
      for (pair<uint64_t, bed_line> strand2: contig.second.second){
        myfile << contig.first << "\t" << strand2.second.start << "\t" << strand2.second.stop << "\t" << "-" << "\t" << strand2.second.coverage << "\t" << strand2.second.bases;
        for (uint64_t hit: strand2.second.hits) {
          myfile << "\t";
          myfile << hit;
        }
        myfile << "\n";
      }
    }
    myfile.close();
  } else{
    cout << "Unable to open file: " << path_to_bed.string() << "\n";
  }
}