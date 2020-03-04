//
// Created by Andrew Bailey on 2019-07-15.
//

#include "MaxKmers.hpp"


/**
Override < operator in order to make the priority queue a min heap

@param a: first eventkmer
@param b: second eventkmer

@return a.posterior_probability > b.posterior_probability
*/
bool operator<(const eventkmer& a, const eventkmer& b) {
  return a.posterior_probability > b.posterior_probability;
}

/**
Override < operator in order to make the priority queue a min heap

@param a: first eventkmer
@param b: second eventkmer

@return a.posterior_probability > b.posterior_probability
*/
bool operator<(const FullSaEvent& a, const FullSaEvent& b) {
  return a.posterior_probability > b.posterior_probability;
}
