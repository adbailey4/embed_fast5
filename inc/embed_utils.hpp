//
// Created by Andrew Bailey on 2019-06-18.
//

#ifndef EMBED_FAST5_UTILS_H
#define EMBED_FAST5_UTILS_H

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <string>

using namespace boost::filesystem;
using namespace std;

namespace embed_utils{
  bool are_characters_in_string(string &characters, string &my_string);
  size_t getFilesize(const std::string& filename);
  int64_t convert_to_int(std::string& str_int);
  path make_dir(path output_path);
  bool compareFiles(const std::string& p1, const std::string& p2);
  bool copyDir(path const & source, path const & destination);
  std::vector<std::string> split(std::string &in, char delimiter);
  float convert_to_float(std::string& str_int);
  string sort_string(string &str);
  vector<string> all_lexicographic_recur(string characters, string data, int last, int index);
  vector<string> all_string_permutations(string characters, int length);
  string remove_duplicate_characters(string input_string);
}

#endif //EMBED_FAST5_UTILS_H
