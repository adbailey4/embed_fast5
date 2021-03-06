//
// Created by Andrew Bailey on 2019-07-15.
//

#include "AssignmentFile.hpp"
#include "EmbedUtils.hpp"
#include <fstream>
#include <vector>


using namespace std;
using namespace boost::coroutines2;
using namespace embed_utils;


AssignmentFile::~AssignmentFile()
= default;

AssignmentFile::AssignmentFile(const string& input_reads_filename){
  this->file_path = input_reads_filename;
  this->k = 0;
}

/**
Create a push type coroutine for parsing an assignment file

 file format:
GCCTTA	t	83.709275	1.000000

*/
void AssignmentFile::assignment_coroutine(event_kmer_coro::push_type& yield){
  std::ifstream in_file(this->file_path.c_str());
  if(in_file.good()) {
    // read the file
    string line;
    while(getline(in_file, line)) {
      vector<string> fields = split_string(line, '\t');
      yield(eventkmer(fields[0], std::stof(fields[2]), fields[1], std::stof(fields[3])));
    }
  } else  {
    cout << "Error loading file: " << this->file_path << "\n";
  }
  in_file.close();
}

/**
Call push type coroutine to create a pull type coroutine
*/
event_kmer_coro::pull_type AssignmentFile::iterate(){
  event_kmer_coro::pull_type seq {bind(&AssignmentFile::assignment_coroutine, this, std::placeholders::_1)};
  return seq;
}

/**
Get the kmer size for a given file
*/
int64_t AssignmentFile::get_k(){
  std::ifstream in_file(this->file_path.c_str());
  if (in_file.good()) {
    // read the file
    std::string line;
    getline(in_file, line);
    std::vector<std::string> fields = split_string(line, '\t');
    this->k = fields[0].length();
  }
  in_file.close();
  return this->k;
}
