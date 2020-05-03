//
// Created by Andrew Bailey on 5/2/20.
//

// embed source
#include "BinaryEventWriter.hpp"
#include "BinaryIO.hpp"
#include "EmbedUtils.hpp"

// boost
#include <boost/filesystem.hpp>

// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace boost::filesystem;
using namespace std;
using namespace embed_utils;


TEST (BinaryFileTest, check1) {
  EXPECT_EQ(1, 1);
}
