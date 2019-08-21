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
  for (uint64_t i=0; i<this->num_locks; i++){
    omp_lock_t new_lock;
    omp_init_lock(&(new_lock));
    locks.push_back(new_lock);
  }
}
/**
Initialize vector of locks
*/
void MarginalizeVariants::initialize_locks(uint64_t number_of_locks) {
  this->num_locks = number_of_locks;
  for (uint64_t i=0; i<this->num_locks; i++){
    omp_lock_t new_lock;
    omp_init_lock(&(new_lock));
    locks.push_back(new_lock);
  }
}


/**
 * Load a vector of variants into the per_genomic position data structure.
 *
 * @param vector_of_calls: vector of variant_calls
 */
void MarginalizeVariants::load_variants(vector<VariantCall>* vector_of_calls){
  bed_line bed_entry;
  int max_index;
  for (auto& i: *vector_of_calls){
    if (this->initialized_locks){
      omp_set_lock(&(this->locks[i.reference_index % this->num_locks]));
    }
    if (i.strand == "+"){
      bed_entry = this->per_genomic_position[i.contig].first[i.reference_index];
    } else {
      bed_entry = this->per_genomic_position[i.contig].second[i.reference_index];
    }

    if (bed_entry.start == -1) {
      bed_entry.start = i.reference_index;
      bed_entry.stop = i.reference_index + 1;
      bed_entry.bases = i.bases;
      size_t n_entries = i.normalized_probs.size();
      bed_entry.hits.resize(i.normalized_probs.size());
      for (size_t j = 0; j < n_entries; ++j) {
        bed_entry.hits.push_back(0.0);
      }
    }
    bed_entry.coverage += 1;
    max_index = max_element(i.normalized_probs.begin(), i.normalized_probs.end()) - i.normalized_probs.begin();
    bed_entry.hits[max_index] += 1;
    if (i.strand == "+"){
      this->per_genomic_position[i.contig].first[i.reference_index] = bed_entry;
    } else {
      this->per_genomic_position[i.contig].second[i.reference_index] = bed_entry;
    }
    if (this->initialized_locks){
      omp_unset_lock(&(this->locks[i.reference_index % this->num_locks]));
    }
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
        myfile << contig.first << "\t" << strand.second.start << "\t" << strand.second.stop << "\t" << "+" << "\t" << strand.second.coverage << "\t" << strand.second.bases << "\t";
        for (uint64_t hit: strand.second.hits) {
          myfile << hit << "\t" ;
        }
        myfile << "\n";
      }
//      second of pair is minus strand
      for (pair<uint64_t, bed_line> strand2: contig.second.second){
        myfile << contig.first << "\t" << strand2.second.start << "\t" << strand2.second.stop << "\t" << "-" << "\t" << strand2.second.coverage << "\t" << strand2.second.bases << "\t";
        for (uint64_t hit: strand2.second.hits) {
          myfile << hit << "\t" ;
        }
        myfile << "\n";
      }
    }
    myfile.close();
  } else{
    cout << "Unable to open file: " << path_to_bed.string() << "\n";
  }
}