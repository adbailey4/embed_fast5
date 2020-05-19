//
// Created by Andrew Bailey on 5/2/20.
//

// embed source
#include "BinaryIO.hpp"
#include "EmbedUtils.hpp"
#include "TestFiles.hpp"
#include "PerPositionKmers.hpp"

// boost
#include <boost/filesystem.hpp>
#include <boost/heap/priority_queue.hpp>

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
  ::close(sequence_file_descriptor);
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
  ::close(sequence_file_descriptor);
}

TEST (BinaryFileTests, test_read_write_vector_of_Events_to_binary) {
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  write
  std::ofstream file_handle = std::ofstream(test_file.string(), std::ofstream::binary);
  vector<Event> write_me;
  write_me.push_back(Event(1, 2));
  write_me.push_back(Event(1.1, 2.2));
  write_vector_to_binary(file_handle, write_me);
  file_handle.close();
// read
  std::ifstream file_handle2 = std::ifstream(test_file.string(), std::ofstream::binary);
  vector<Event> read_to_me;
  read_vector_from_binary(file_handle2, read_to_me, 2);
  file_handle2.close();
  EXPECT_EQ(write_me[0].posterior_probability, read_to_me[0].posterior_probability);
  EXPECT_EQ(write_me[1].posterior_probability, read_to_me[1].posterior_probability);
  EXPECT_EQ(write_me[0].descaled_event_mean, read_to_me[0].descaled_event_mean);
  EXPECT_EQ(write_me[1].descaled_event_mean, read_to_me[1].descaled_event_mean);
//  pread
  int sequence_file_descriptor = ::open(test_file.c_str(), O_RDONLY);
  vector<Event> read_to_me2;
  off_t byte_index = off_t(0);
  pread_vector_from_binary(sequence_file_descriptor, read_to_me2, 2, byte_index);
  EXPECT_EQ(write_me[0].posterior_probability, read_to_me2[0].posterior_probability);
  EXPECT_EQ(write_me[1].posterior_probability, read_to_me2[1].posterior_probability);
  EXPECT_EQ(write_me[0].descaled_event_mean, read_to_me2[0].descaled_event_mean);
  EXPECT_EQ(write_me[1].descaled_event_mean, read_to_me2[1].descaled_event_mean);
  ::close(sequence_file_descriptor);
}

//template <class T, class S, class C>
//S& Container(boost::heap::priority_queue<T, S, C>& q) {
//  struct HackedQueue : private boost::heap::<T, S, C> {
//    static S& Container(priority_queue<T, S, C>& q) {
//      return q.*&HackedQueue::q_;
//    }
//  };
//  return HackedQueue::Container(q);
//}


TEST (BinaryFileTests, test_read_write_queue_of_Events_to_binary) {
  path tempdir = temp_directory_path();
  path test_file = tempdir / "test.event";
  if (exists(test_file)){
    remove(test_file);
  }
//  write
  std::ofstream file_handle = std::ofstream(test_file.string(), std::ofstream::binary);
  Kmer write_me;
  write_me.add_event(1, 2);
  write_me.add_event(2, 3);
  write_vector_to_binary(file_handle, write_me.events);
  file_handle.close();
// read
  std::ifstream file_handle2 = std::ifstream(test_file.string(), std::ofstream::binary);
  vector<Event> read_to_me;
  read_vector_from_binary(file_handle2, read_to_me, 2);
  file_handle2.close();
  EXPECT_EQ(2, read_to_me[0].posterior_probability);
  EXPECT_EQ(3, read_to_me[1].posterior_probability);
  EXPECT_EQ(1, read_to_me[0].descaled_event_mean);
  EXPECT_EQ(2, read_to_me[1].descaled_event_mean);
//  pread
  int sequence_file_descriptor = ::open(test_file.c_str(), O_RDONLY);
  vector<Event> read_to_me2;
  off_t byte_index = off_t(0);
  pread_vector_from_binary(sequence_file_descriptor, read_to_me2, 2, byte_index);
  EXPECT_EQ(2, read_to_me[0].posterior_probability);
  EXPECT_EQ(3, read_to_me[1].posterior_probability);
  EXPECT_EQ(1, read_to_me[0].descaled_event_mean);
  EXPECT_EQ(2, read_to_me[1].descaled_event_mean);
  ::close(sequence_file_descriptor);
}