//
// Created by Andrew Bailey on 5/2/20.
//

// embed source
#include "BinaryIO.hpp"
#include "EmbedUtils.hpp"
#include "TestFiles.hpp"

// boost
#include <boost/filesystem.hpp>

// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace boost::filesystem;
using namespace std;
using namespace embed_utils;
using namespace test_files;


TEST (BinaryFileTests, test_read_write_string_to_binary) {
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  write
  std::ofstream file_handle = std::ofstream(test_file.string(), std::ofstream::binary);
  string write_me = "hello!";
  write_string_to_binary(file_handle, write_me);
  file_handle.close();
// read
  std::ifstream file_handle2 = std::ifstream(test_file.string(), std::ofstream::binary);

  string read_to_me;
  read_string_from_binary(file_handle2, read_to_me, 6);
  EXPECT_EQ(write_me, read_to_me);
  EXPECT_EQ("hello!", read_to_me);
  file_handle2.close();
//  pread
  int sequence_file_descriptor = ::open(test_file.c_str(), O_RDONLY);
  string read_to_me2;
  off_t byte_index = off_t(0);
  pread_string_from_binary(sequence_file_descriptor, read_to_me2, 6, byte_index);
  EXPECT_EQ("hello!", read_to_me2);
  ::close(sequence_file_descriptor);
}

TEST (BinaryFileTests, test_read_write_vector_to_binary) {
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  write
  std::ofstream file_handle = std::ofstream(test_file.string(), std::ofstream::binary);
  vector<int> write_me = {2, 3};
  write_vector_to_binary(file_handle, write_me);
  file_handle.close();
// read
  std::ifstream file_handle2 = std::ifstream(test_file.string(), std::ofstream::binary);
  vector<int> read_to_me;
  read_vector_from_binary(file_handle2, read_to_me, 2);
  EXPECT_EQ(write_me[0], read_to_me[0]);
  EXPECT_EQ(write_me[1], read_to_me[1]);
  file_handle2.close();
//  pread
  int sequence_file_descriptor = ::open(test_file.c_str(), O_RDONLY);
  vector<int> read_to_me2;
  off_t byte_index = off_t(0);
  pread_vector_from_binary(sequence_file_descriptor, read_to_me2, 2, byte_index);
  EXPECT_EQ(write_me[0], read_to_me2[0]);
  EXPECT_EQ(write_me[1], read_to_me2[1]);

}

TEST (BinaryFileTests, test_read_write_value_to_binary) {
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  write
  std::ofstream file_handle = std::ofstream(test_file.string(), std::ofstream::binary);
  int write_me = 2;
  write_value_to_binary(file_handle, write_me);
  file_handle.close();
// read
  std::ifstream file_handle2 = std::ifstream(test_file.string(), std::ofstream::binary);
  int read_to_me;
  read_value_from_binary(file_handle2, read_to_me);
  EXPECT_EQ(write_me, read_to_me);
  file_handle2.close();
  //  pread
  int sequence_file_descriptor = ::open(test_file.c_str(), O_RDONLY);
  int read_to_me2;
  off_t byte_index = off_t(0);
  pread_value_from_binary(sequence_file_descriptor, read_to_me2, byte_index);
  EXPECT_EQ(write_me, read_to_me2);
}
