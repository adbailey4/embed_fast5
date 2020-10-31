//
// Created by Andrew Bailey on 6/7/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_AMBIGMODELTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_AMBIGMODELTESTS_HPP_

// embed source
#include "TestFiles.hpp"
#include "AmbigModel.hpp"
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace std;
using namespace embed_utils;
using namespace test_files;


TEST (AmbigModelTests, test_get_canonical_kmers) {
//  Redirect a(true, true);
  AmbigModel am(AMBIG_MOD_FILE.string());
  std::set<string> test = am.get_canonical_kmers("ATGaAa");
  std::set<string> answer = {"ATGAAA"};
  ASSERT_THAT(answer, ElementsAreArray(test));
  test = am.get_canonical_kmers("ATGAAa");
  ASSERT_THAT(answer, ElementsAreArray(test));
  std::set<string> a{"ATGaAa", "ATGAAa"};
  test = am.get_canonical_kmers(a);
  ASSERT_THAT(answer, ElementsAreArray(test));
}



#endif //EMBED_FAST5_TESTS_SRC_AMBIGMODELTESTS_HPP_
