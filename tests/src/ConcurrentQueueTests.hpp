//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_CONCURRENTQUEUETESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_CONCURRENTQUEUETESTS_HPP_

// embed source
#include "ConcurrentQueue.hpp"
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>
// Standard Libray
#include <future>

using namespace std;

template<typename T>
void add_to_queue(ConcurrentQueue<T>& cq){
  cq.push("Test");
}

template<typename T>
void remove_from_queue(ConcurrentQueue<T>& cq, std::promise<T> && p){
  T value;
  cq.wait_and_pop(value);
  p.set_value(value);
}


TEST (ConcurrentQueueTests, test_ConcurrentQueue) {
  ConcurrentQueue<string> cq{};
  std::promise<string> p;
  auto f = p.get_future();
  std::thread t2(remove_from_queue<string>, ref(cq), std::move(p));
  std::thread t1(add_to_queue<string>, ref(cq));
  t1.join();
  t2.join();
  string i = f.get();
  EXPECT_EQ("Test", i);
}

#endif //EMBED_FAST5_TESTS_SRC_CONCURRENTQUEUETESTS_HPP_
