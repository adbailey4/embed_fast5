//
// Created by Andrew Bailey on 2019-06-18.
//

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include "EmbedUtils.hpp"
#include <cmath>

using namespace boost::filesystem;
using namespace std;


namespace embed_utils {

/**
 * Get the size of a file. https://techoverflow.net/2013/08/21/how-to-get-filesize-using-stat-in-cc/
 * @param filename The name of the file to check size for
 * @return The filesize, or 0 if the file does not exist.
 */
    size_t getFilesize(const std::string &filename) {
        struct stat st;
        if (stat(filename.c_str(), &st) != 0) {
            return 0;
        }
        return st.st_size;
    }

//https://www.systutorials.com/131/convert-string-to-int-and-reverse/
    int64_t convert_to_int(std::string &str_int) {
        int64_t number;
        std::istringstream iss(str_int);
        iss >> number;
        if (!iss.good()) {
            return number;
        }
    }

    float convert_to_float(std::string &str_int) {
        std::string::size_type sz;     // alias of size_t
        float number = std::stof(str_int, &sz);
        return number;
    }


    bool are_characters_in_string(string &characters, string &my_string) {
        for (char &c : characters) {
            if (my_string.find(c) != std::string::npos) {
                return true;
            }
        }
        return false;
    }


    path make_dir(path output_path) {
        if (!exists(output_path)) {
            create_directory(output_path);
        }
        return output_path;
    }

    bool copyDir(
            boost::filesystem::path const &source,
            boost::filesystem::path const &destination
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
    bool compareFiles(const std::string &p1, const std::string &p2) {
        std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
        std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

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
    std::vector<std::string> split(std::string &in, char delimiter) {
        std::vector<std::string> out;
        size_t lastPos = 0;
        size_t pos = in.find_first_of(delimiter);

        while (pos != std::string::npos) {
            out.push_back(in.substr(lastPos, pos - lastPos));
            lastPos = pos + 1;
            pos = in.find_first_of(delimiter, lastPos);
        }
        out.push_back(in.substr(lastPos));
        return out;
    }
//    https://www.geeksforgeeks.org/sort-string-characters/
    string sort_string(string &str)
    {
        sort(str.begin(), str.end());
        return str;
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
  vector<string> all_lexicographic_recur(string characters, string data, int last, int index){

  //    # One by one fix all characters at the given index and
  //    # recur for the subsequent indexes
    size_t length = characters.length();
    vector<string> kmers;


    for (size_t i = 0; i < length; i++){

  //    # Fix the ith character at index and if this is not
  //    # the last index then recursively call for higher
  //    # indexes
      data[index] = characters[i];
  //    # If this is the last index then print the string
  //    # stored in data[]
      if (index == last) {
        kmers.push_back(data);

      } else{
        vector<string> tmp_kmers = all_lexicographic_recur(characters, data, last, index + 1);
        for (auto & element : tmp_kmers){
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
  vector<string> all_string_permutations(string characters, int length){
    vector<string> kmers;
    assert(length >= 0);
    assert(characters.length() > 0);
    string data;

    for (size_t i = 0; i < length; i++){
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
  string remove_duplicate_characters(string input_string){

    size_t character_map[128] = {0};
    string output_string;
    size_t length = input_string.length();

    for (size_t i=0; i < length; ++i){
      char str_char = input_string[i];
      if (character_map[int(str_char)] == 0){
        output_string += str_char;
      }
      ++character_map[int(str_char)];
    }
    return output_string;
  }
}

