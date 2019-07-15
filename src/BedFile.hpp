//
// Created by Andrew Bailey on 2019-07-12.
//

#ifndef EMBED_FAST5_IMPL_BED_FILE_H_
#define EMBED_FAST5_IMPL_BED_FILE_H_

#include <vector>
#include <string>
#include <boost/unordered_map.hpp>

using namespace std;

struct bed_line{
  uint64_t start;
  uint64_t stop;
  uint64_t coverage;
  uint64_t hits;
  string base;

  uint64_t color1;
  uint64_t color2;
  uint64_t color3;

};


class BedFile {

  boost::unordered_map<std::string, vector<vector<bed_line>, vector<bed_line>>> all_data;

};

#endif //EMBED_FAST5_IMPL_BED_FILE_H_
