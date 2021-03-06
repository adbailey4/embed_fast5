//
// Created by Andrew Bailey on 03/15/19.
//

// embed tests
#include "BinaryIOTests.hpp"
#include "ConcurrentQueueTests.hpp"
#include "EmbedUtilsTests.hpp"
#include "Fast5Tests.hpp"
#include "AmbigModelTests.hpp"
#include "FileTests.hpp"
#include "FilterAlignmentsTests.hpp"
#include "MarginalizeVariantsTests.hpp"
#include "MaxKmersTests.hpp"
#include "TopKmersTests.hpp"
#include "VariantPathTests.hpp"
#include "PerPositionKmersTests.hpp"
#include "BaseKmerTests.hpp"
#include "BinaryEventTests.hpp"

// boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>
// Standard Libray
#include <iostream>

using namespace test_files;

int main(int argc, char **argv) {
  H5Eset_auto(0, nullptr, nullptr);
//  deal with paths to test files
  ::testing::InitGoogleMock(&argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  std::cout << argc - 1;
  if (argc == 2){
    path home = argv[1];
    prepend_home(home);
  } else{
    throw_assert(argc == 2, "There is an incorrect number of arguments to tests")
  }
//  google test stuff
//  specific tests to run
//  ::testing::GTEST_FLAG(filter) = "EmbedUtilsTests*";

//  ::testing::GTEST_FLAG(filter) = "BaseKmerTests*";
//  ::testing::GTEST_FLAG(filter) = "PerPositionKmersTests*";
//  ::testing::GTEST_FLAG(filter) = "TopKmersTests*";
//  ::testing::GTEST_FLAG(filter) = "EmbedUtilsTests*";
//  ::testing::GTEST_FLAG(filter) = "BinaryEventTests*";
//  ::testing::GTEST_FLAG(filter) = "AmbigModelTests*";

  return RUN_ALL_TESTS();
}
