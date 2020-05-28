//
// Created by Andrew Bailey on 2019-06-18.
//

// Embed
#include "EmbedUtils.hpp"

// Boost libraries.
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

// Standard library.
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <functional>
#include <sstream>
#include <cmath>
#include <map>
#include <chrono>
#include <tuple>


using namespace boost::filesystem;
using namespace std;
using namespace std::chrono;


namespace embed_utils {

/**
 * Get the size of a file. https://techoverflow.net/2013/08/21/how-to-get-filesize-using-stat-in-cc/
 * @param filename The name of the file to check size for
 * @return The filesize, or 0 if the file does not exist.
 */
size_t get_file_size(const path &filename) {
  struct stat st{};
  if (stat(filename.c_str(), &st) != 0) {
    return 0;
  }
  return st.st_size;
}

//https://www.systutorials.com/131/convert-string-to-int-and-reverse/
int64_t string_to_int(std::string &str_int) {
  int64_t number;
  std::istringstream iss(str_int);
  iss >> number;
  if (!iss.good()) {
    return number;
  }
}

float string_to_float(const string &str_int) {
  return std::stof(str_int);
}

bool are_characters_in_string(string &characters, string &my_string) {
  for (char &c : characters) {
    if (my_string.find(c) != std::string::npos) {
      return true;
    }
  }
  return false;
}


bool copyDir(
    const boost::filesystem::path &source,
    const boost::filesystem::path &destination
) {
  namespace fs = boost::filesystem;
  try {
    // Check whether the function call is valid
    if (
        !fs::exists(source) ||
            !fs::is_directory(source)
        ) {
      std::cerr << "Source directory " << source.string()
                << " does not exist or is not a directory." << '\n';
      return false;
    }
    if (fs::exists(destination)) {
      std::cerr << "Destination directory " << destination.string()
                << " already exists." << '\n';
      return false;
    }
    // Create the destination directory
    if (!fs::create_directory(destination)) {
      std::cerr << "Unable to create destination directory"
                << destination.string() << '\n';
      return false;
    }
  }
  catch (fs::filesystem_error const &e) {
    std::cerr << e.what() << '\n';
    return false;
  }
  // Iterate through the source directory
  for (
      fs::directory_iterator file(source);
      file != fs::directory_iterator(); ++file
      ) {
    try {
      fs::path current(file->path());
      if (fs::is_directory(current)) {
        // Found directory: Recursion
        if (
            !copyDir(
                current,
                destination / current.filename()
            )
            ) {
          return false;
        }
      } else {
        // Found file: Copy
        fs::copy_file(
            current,
            destination / current.filename()
        );
      }
    }
    catch (fs::filesystem_error const &e) {
      std::cerr << e.what() << '\n';
    }
  }
  return true;
}

//https://stackoverflow.com/questions/6163611/compare-two-files
bool compare_files(const path &p1, const path &p2) {
  throw_assert(is_regular_file(p1), "Path 1 file does not exist:" + p1.string())
  throw_assert(is_regular_file(p2), "Path 2 file does not exist:" + p2.string())

  std::ifstream f1(p1.string(), std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(p2.string(), std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; //file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false; //size mismatch
  }

  //seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

// Split a string into parts based on the delimiter
std::vector<std::string> split_string(std::string &in, char delimiter) {
  std::vector<std::string> out;
  size_t lastPos = 0;
  size_t pos = in.find_first_of(delimiter);

  while (pos != std::string::npos) {
    string split = in.substr(lastPos, pos - lastPos);
    out.push_back(split);
    lastPos = pos + 1;
    pos = in.find_first_of(delimiter, lastPos);
  }
  string split = in.substr(lastPos);
  out.push_back(split);
  return out;
}

//    https://www.geeksforgeeks.org/sort-string-characters/
string sort_string(string &str) {
  sort(str.begin(), str.end());
  return str;
}

/**
Create a directory if it does not exist
@param output_path: path to directory.
@return path to directory that definitely exists
*/
path make_dir(path &output_path) {
  if (!exists(output_path)) {
    create_directory(output_path);
  }
  return output_path;
}

/**
Generator for all permutations of an alphabetically sorted string
source: https://www.geeksforgeeks.org/print-all-permutations-with-repetition-of-characters/

@param characters: alphabetically sorted string.
@param data: characters
@param last: length of string wanted to iterate
@param index: current index to fix
@return vector of strings
*/
vector<string> all_lexicographic_recur(string &characters, string &data, uint64_t last, uint64_t index) {

  //    # One by one fix all characters at the given index and
  //    # recur for the subsequent indexes
  size_t length = characters.length();
  vector<string> kmers;

  for (size_t i = 0; i < length; i++) {

    //    # Fix the ith character at index and if this is not
    //    # the last index then recursively call for higher
    //    # indexes
    data[index] = characters[i];
    //    # If this is the last index then print the string
    //    # stored in data[]
    if (index == last) {
      kmers.push_back(data);

    } else {
      vector<string> tmp_kmers = all_lexicographic_recur(characters, data, last, index + 1);
      for (auto &element : tmp_kmers) {
        kmers.push_back(element);
      }
    }
  }
  return kmers;
}

/**
Creates an alphabetically sorted vector all string permutations of a set of characters in a string

@param characters: characters to use.
@param length: kmer length
@return vector of strings
*/
vector<string> all_string_permutations(string &characters, int &length) {
  assert(length >= 0);
  assert(characters.length() > 0);
  string data;

  for (size_t i = 0; i < length; i++) {
    data += " ";
  }

  characters = remove_duplicate_characters(characters);
  sort_string(characters);

  return all_lexicographic_recur(characters, data, length - 1, 0);

}
/**
Remove duplicate characters. Returns new string

@param input_string: input string to remove characters from
@return new string.
*/
string remove_duplicate_characters(string &input_string) {

  size_t character_map[128] = {0};
  string output_string;
  size_t length = input_string.length();

  for (size_t i = 0; i < length; ++i) {
    char str_char = input_string[i];
    if (character_map[int(str_char)] == 0) {
      output_string += str_char;
    }
    ++character_map[int(str_char)];
  }
  return output_string;
}

/**
Directory coroutine to push paths with a given extension

@param yield: coroutine pushtype
@param directory: path to input directory
@param ext: string for the extension to check files

@return yields a path to a file with extension.

*/
void dir_iterator_coroutine(dir_coro::push_type &yield, path &directory, string &ext) {
  directory_iterator end_itr;
  for (directory_iterator itr(directory); itr != end_itr; ++itr) {
    //        filter for files that are regular, end with ext and are not empty
    if ((is_regular_file(itr->path()) and get_file_size(itr->path().string()) > 0) and
        (ext.empty() or itr->path().extension().string() == ext)) {
      yield(itr->path());
    }
  }
}

/**
 Lists all non empty files in directory with extension. If extension is not passed in then there is no extension
 filtering

@param directory: path to input directory
@param ext: string for the extension to check files

@return yields a path to a file with extension.
*/
dir_coro::pull_type list_files_in_dir(path &directory, string &ext) {
  dir_coro::pull_type file{bind(&dir_iterator_coroutine, std::placeholders::_1, directory, ext)};
  return file;
}

/**
 Create a map of all ambiguous bases so that we can determine what each ambiguous character represents

@return map from base to ambiguous bases.
*/
std::map<string, string> create_ambig_bases() {

  std::map<string, string> ambig_hash;
  ambig_hash.insert(std::pair<string, string>("R", "AG"));
  ambig_hash.insert(std::pair<string, string>("Y", "CT"));
  ambig_hash.insert(std::pair<string, string>("S", "CG"));
  ambig_hash.insert(std::pair<string, string>("W", "AT"));
  ambig_hash.insert(std::pair<string, string>("K", "GT"));
  ambig_hash.insert(std::pair<string, string>("M", "AC"));
  ambig_hash.insert(std::pair<string, string>("B", "CGT"));
  ambig_hash.insert(std::pair<string, string>("D", "AGT"));
  ambig_hash.insert(std::pair<string, string>("H", "ACT"));
  ambig_hash.insert(std::pair<string, string>("V", "ACG"));
  ambig_hash.insert(std::pair<string, string>("X", "ACGT"));
  ambig_hash.insert(std::pair<string, string>("L", "CEO"));
  ambig_hash.insert(std::pair<string, string>("P", "CE"));
  ambig_hash.insert(std::pair<string, string>("Q", "AI"));
  ambig_hash.insert(std::pair<string, string>("f", "AF"));
  ambig_hash.insert(std::pair<string, string>("U", "ACEGOT"));
  ambig_hash.insert(std::pair<string, string>("Z", "JT"));
  ambig_hash.insert(std::pair<string, string>("j", "Tp"));
  ambig_hash.insert(std::pair<string, string>("k", "Gb"));
  ambig_hash.insert(std::pair<string, string>("l", "Gd"));
  ambig_hash.insert(std::pair<string, string>("m", "Ce"));
  ambig_hash.insert(std::pair<string, string>("n", "Th"));
  ambig_hash.insert(std::pair<string, string>("o", "Ai"));

  ambig_hash.insert(std::pair<string, string>("AG", "R"));
  ambig_hash.insert(std::pair<string, string>("CT", "Y"));
  ambig_hash.insert(std::pair<string, string>("CG", "S"));
  ambig_hash.insert(std::pair<string, string>("AT", "W"));
  ambig_hash.insert(std::pair<string, string>("GT", "K"));
  ambig_hash.insert(std::pair<string, string>("AC", "M"));
  ambig_hash.insert(std::pair<string, string>("CGT", "B"));
  ambig_hash.insert(std::pair<string, string>("AGT", "D"));
  ambig_hash.insert(std::pair<string, string>("ACT", "H"));
  ambig_hash.insert(std::pair<string, string>("ACG", "V"));
  ambig_hash.insert(std::pair<string, string>("ACGT", "X"));
  ambig_hash.insert(std::pair<string, string>("CEO", "L"));
  ambig_hash.insert(std::pair<string, string>("CE", "P"));
  ambig_hash.insert(std::pair<string, string>("AI", "Q"));
  ambig_hash.insert(std::pair<string, string>("AF", "f"));
  ambig_hash.insert(std::pair<string, string>("ACEGOT", "U"));
  ambig_hash.insert(std::pair<string, string>("JT", "Z"));
  ambig_hash.insert(std::pair<string, string>("Tp", "j"));
  ambig_hash.insert(std::pair<string, string>("Gb", "k"));
  ambig_hash.insert(std::pair<string, string>("Gd", "l"));
  ambig_hash.insert(std::pair<string, string>("Ce", "m"));
  ambig_hash.insert(std::pair<string, string>("Th", "n"));
  ambig_hash.insert(std::pair<string, string>("Ai", "o"));

  return ambig_hash;
}

/**
 Create a map of all ambiguous bases so that we can determine what each ambiguous character represents

@return map from base to ambiguous bases.
*/
std::map<string, string> create_ambig_bases2(string config_file) {
  std::map<string, string> ambig_hash;
  if (!config_file.empty()){
    char encoding[10];
    char ambig_bases[10];
    char line[100];
    FILE *infile = fopen(config_file.c_str(), "r");
    throw_assert(infile, "Couldn't open " + string(config_file) + " for reading\n")
    int i = 0;
    while(i < 300 && fgets(line, sizeof(line), infile) != nullptr){
      sscanf(line, "%s\t%s", encoding, ambig_bases);
      ambig_hash.insert(std::pair<string, string>(encoding, ambig_bases));
      i++;
    }
  } else {
    ambig_hash = create_ambig_bases();
  }
  return ambig_hash;
}


/**
  Time and execute a function which returns void.

@param a bound function using std::bind
@return string formatted with hours, miutes, seconds to the microsecond
*/
string get_time_string(std::function<void()> bound_function) {
  string output;
  tuple<uint64_t, uint64_t, uint64_t, uint64_t> data = get_time(std::move(bound_function));
  output = "hours: " + to_string(get<0>(data)) + " minutes: " + to_string(get<1>(data)) + " seconds: "
      + to_string(get<2>(data)) + "." + to_string(get<3>(data)) + "\n";
  return output;
}

/**
Time and execute a function which returns void.

@param a bound function using std::bind
@return tuple of uint64_t's [hours, minutes, seconds, microseconds]
*/
tuple<uint64_t, uint64_t, uint64_t, uint64_t> get_time(std::function<void()> bound_function) {
  auto start = high_resolution_clock::now();
  bound_function();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  uint64_t hours = floor(duration.count() / 3600000000);
  uint64_t minutes = floor((duration.count() - (hours * 3600000000)) / 60000000);
  uint64_t seconds = floor((duration.count() - (hours * 3600000000) - (minutes * 60000000)) / 1000000);
  uint64_t microseconds = floor((duration.count() - (hours * 3600000000) - (minutes * 60000000) - (seconds * 1000000)));
  tuple<uint64_t, uint64_t, uint64_t, uint64_t> output(hours, minutes, seconds, microseconds);
  return output;
}

/**
Counts number of lines in a file.

@param file_path: path to file
@return number of lines in file
*/
int64_t lines_in_file(path &file_path) {
  unsigned int number_of_lines = 0;
  FILE *infile = fopen(file_path.c_str(), "r");
  int ch;
  while (EOF != (ch = getc(infile))) {
    if ('\n' == ch) {
      ++number_of_lines;
    }
  }
  fclose(infile);
  return number_of_lines;
}

/**
 * Calculate number of columns in file
 *
 * @param file_path: path to file
 */
uint64_t number_of_columns(const path &file_path, char sep){
  FILE *infile = fopen(file_path.c_str(), "r");
  uint64_t n_col = 0;
  throw_assert((get_file_size(file_path) > 0), "File is empty");
  int ch;
  while (EOF != (ch = getc(infile))) {
    if (sep == ch) {
      ++n_col;
    }
    if ('\n' == ch) {
      break;
    }
  }
  if (n_col != 0){
    ++n_col;
  }
  fclose(infile);
  return n_col;
}

}
