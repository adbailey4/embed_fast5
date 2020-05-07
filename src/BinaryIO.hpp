//
// Created by Ryan Lorig-Roach https://github.com/rlorigro/runlength_analysis_cpp/blob/master/inc/BinaryIO.hpp.
//

#ifndef EMBED_FAST5_SRC_BINARYIO_HPP_
#define EMBED_FAST5_SRC_BINARYIO_HPP_

#include <iostream>
#include <ostream>
#include <istream>
#include <string>
#include <vector>
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

using std::ostream;
using std::istream;
using std::string;
using std::cout;
using std::cerr;
using std::vector;
using std::runtime_error;


void write_string_to_binary(ostream& file, string& s){
  ///
  /// Without worrying about size conversions, write any string to a file using ostream.write
  ///

  file.write(reinterpret_cast<const char*>(s.data()), s.size());
}

void read_string_from_binary(istream& file, string& s, uint64_t length){
  ///
  /// Without worrying about size conversions, read any value to a file using ostream.write
  ///

  s.resize(length);
  file.read(const_cast<char *>(reinterpret_cast<const char *>(s.data())), length);
}

template<class T> void write_vector_to_binary(ostream& s, const vector<T>& v){
  ///
  /// Without worrying about size conversions, write any vector to a file using ostream.write
  ///

  s.write(reinterpret_cast<const char*>(v.data()), v.size()*sizeof(T));
}

template<class T> void read_vector_from_binary(istream& s, vector<T>& v, uint64_t length){
  ///
  /// Without worrying about size conversions, read any vector from a file using istream.read
  ///
//  cout << "Reading vector of size: " << sizeof(T)*length << " at position: " << s.tellg() << '\n';

  v.resize(length);
  s.read(reinterpret_cast<char*>(v.data()), sizeof(T)*length);
}

template<class T> void write_value_to_binary(ostream& s, T v){
  ///
  /// Without worrying about size conversions, write any value to a file using ostream.write
  ///

  auto v_temp = v;
  s.write(reinterpret_cast<const char*>(&v_temp), sizeof(T));
}

template<class T> void read_value_from_binary(istream& s, T& v){
  ///
  /// Without worrying about size conversions, read any value from a file using istream.read
  ///
//  cout << "Reading value size of: " << sizeof(T) << " at position: " << s.tellg() << '\n';
  s.read(reinterpret_cast<char*>(&v), sizeof(T));
}

void pread_bytes(int file_descriptor, char* buffer_pointer, size_t bytes_to_read, off_t& byte_index){
  ///
  /// Reimplementation of binary read_bytes(), but with Linux pread, which is threadsafe
  ///

  while (bytes_to_read) {
    const ssize_t byte_count = ::pread(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
    if (byte_count <= 0) {
      throw runtime_error("ERROR " + std::to_string(errno) + " while reading: " + string(::strerror(errno)));
    }
    bytes_to_read -= byte_count;
    buffer_pointer += byte_count;
    byte_index += byte_count;
  }
}

void pread_string_from_binary(int file_descriptor, string& s, uint64_t length, off_t& byte_index){
  ///
  /// Reimplementation of binary read_string_from_binary(), but with Linux pread, which is threadsafe
  ///

  s.resize(length);

  size_t bytes_to_read = length;
  char* buffer_pointer = const_cast<char *>(reinterpret_cast<const char *>(s.data()));

  pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}


template<class T> void pread_value_from_binary(int file_descriptor,  T& v, off_t& byte_index){
  ///
  /// Reimplementation of binary read_value_from_binary(), but with Linux pread, which is threadsafe
  ///

  size_t bytes_to_read = sizeof(T);
  char* buffer_pointer = reinterpret_cast<char*>(&v);

  pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}


template<class T> void pread_vector_from_binary(int file_descriptor, vector<T>& v, uint64_t length, off_t& byte_index){
  ///
  /// Reimplementation of binary read_vector_from_binary(), but with Linux pread, which is threadsafe
  ///

  v.resize(length);

  size_t bytes_to_read = sizeof(T)*length;
  char* buffer_pointer = reinterpret_cast<char*>(v.data());

  pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}


#endif //EMBED_FAST5_SRC_BINARYIO_HPP_
