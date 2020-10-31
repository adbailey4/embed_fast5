//
// Created by Andrew Bailey on 4/10/20.
//

#ifndef EMBED_FAST5_SRC_CONCURRENTQUEUE_HPP_
#define EMBED_FAST5_SRC_CONCURRENTQUEUE_HPP_

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

/**
 * Create a concurrent queue which can be written to and read from by multiple threads
 * source: https://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html

 * @tparam Data: data class for queue
 */

template<typename Data>
class ConcurrentQueue {
 private:
  std::queue<Data> the_queue;
  mutable std::mutex the_mutex;
  std::condition_variable the_condition_variable;

 public:
  ConcurrentQueue() : no_additional_data(false) {}
  std::atomic<bool> no_additional_data;

  void push(Data const& data)
  {
    std::unique_lock<std::mutex> lock(the_mutex);
    the_queue.push(data);
    lock.unlock();
    the_condition_variable.notify_one();
  }

  bool empty() const
  {
    std::unique_lock<std::mutex> lock(the_mutex);
    return the_queue.empty();
  }

  bool try_pop(Data& popped_value)
  {
    std::unique_lock<std::mutex> lock(the_mutex);
    if(the_queue.empty())
    {
      return false;
    }

    popped_value=the_queue.front();
    the_queue.pop();
    return true;
  }

  bool wait_and_pop(Data& popped_value)
  {
    std::unique_lock<std::mutex> lock(the_mutex);
    while(the_queue.empty())
    {
      if (no_additional_data){
        return false;
      }
      the_condition_variable.wait(lock);
    }
    popped_value=the_queue.front();
    the_queue.pop();
    return true;
  }

  void stop(){
    no_additional_data = true;
    the_condition_variable.notify_all();
  }


};

#endif //EMBED_FAST5_SRC_CONCURRENTQUEUE_HPP_
