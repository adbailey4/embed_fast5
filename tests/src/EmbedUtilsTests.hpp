//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_EMBEDUTILSTESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_EMBEDUTILSTESTS_HPP_

// embed lib
#include "EmbedUtils.hpp"
// embed test files
#include "TestFiles.hpp"
// boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace embed_utils;
using namespace test_files;
using namespace boost::filesystem;
using ::testing::ElementsAreArray;


TEST (EmbedUtilsTests, test_lines_in_file){
  Redirect a(true, true);
  path bed_file = TEST_FILES / "bed_files/test.bed";
  EXPECT_EQ(lines_in_file(bed_file), 2);
}

TEST (EmbedUtilsTests, test_filter_emtpy_files){
  Redirect a(true, true);
  path bed_file = TEST_FILES / "bed_files/test.bed";
  path empty_bed = TEST_FILES / "bed_files/empty_test.bed";
  path non_empty_csv = TEST_FILES / "bed_files/other_test.csv";

  vector<string> bed_paths1 = {bed_file.string()};
  vector<string> bed_paths2 = {bed_file.string(), non_empty_csv.string()};
  vector<string> bed_paths3 = {bed_file.string(), empty_bed.string()};
  vector<path> filtered_bed_paths1 = filter_emtpy_files(bed_paths1, ".bed");
  vector<path> filtered_bed_paths2 = filter_emtpy_files(bed_paths2, ".bed");
  vector<path> filtered_bed_paths3 = filter_emtpy_files(bed_paths3, ".bed");

  EXPECT_EQ(filtered_bed_paths1.size(), 1);
  EXPECT_EQ(filtered_bed_paths2.size(), 1);
  EXPECT_EQ(filtered_bed_paths3.size(), 1);
}

TEST (EmbedUtilsTests, test_number_of_columns){
  Redirect a(true, true);
  path bed_file = TEST_FILES / "bed_files/test.bed";
  path empty_bed = TEST_FILES / "bed_files/empty_test.bed";
  path non_empty_csv = TEST_FILES / "bed_files/other_test.csv";
  string test_bed_string = bed_file.string();

  uint64_t n_col1 = number_of_columns(bed_file.string(), '\t');
  uint64_t n_col3 = number_of_columns(non_empty_csv.string(), ',');
  uint64_t n_col4 = number_of_columns(bed_file.string(), ',');
  uint64_t n_col5 = number_of_columns(test_bed_string, '\t');

  ASSERT_THROW(number_of_columns(empty_bed.string(), '\t'), AssertionFailureException);
  EXPECT_EQ(n_col1, 10);
  EXPECT_EQ(n_col3, 2);
  EXPECT_EQ(n_col4, 0);
  EXPECT_EQ(n_col5, 10);
}

TEST (EmbedUtilsTests, test_compare_files){
  Redirect a(true, true);
  path bed_file = TEST_FILES / "bed_files/test.bed";
  path fake_file = TEST_FILES / "bed_files/asdf";
  ASSERT_THROW(compare_files(fake_file, bed_file), AssertionFailureException);
  EXPECT_TRUE(compare_files(bed_file, bed_file));
}

TEST (EmbedUtilsTests, test_Redirect){
  Redirect a(true, true);
  string test1 = "Test";
  cout << test1;
  cerr << test1;
  string output = a.get_cout();
  string error = a.get_cerr();
  EXPECT_EQ(output, test1);
  EXPECT_EQ(error, test1);
}

TEST (EmbedUtilsTests, test_sort_string){
  Redirect a(true, true);
  string something = "something";
  sort_string(something);
  EXPECT_TRUE(something.compare("eghimnost")==0);
}

TEST (EmbedUtilsTests , test_all_lexicographic_recur){
  Redirect a(true, true);
  string alphabet = "AC";
  int length = 2;
  string data = "  ";
  vector<string> kmers = all_lexicographic_recur(alphabet, data, length-1, 0);
  EXPECT_EQ(kmers.size(), 4);
}

TEST (EmbedUtilsTests, test_split){
  Redirect a(true, true);
  string csv = "asdf,asdf,sdf,df";
  string tsv = "asdf\tasdf\tsdf\tdf";
  vector<string> split_answer{"asdf","asdf","sdf","df"};
  vector<string> something2 = split_string(csv, ',');

  ASSERT_THAT(split_answer, ElementsAreArray(something2));

  vector<string> something4 = split_string(tsv, '\t');

  ASSERT_THAT(split_answer, ElementsAreArray(something4));
}

TEST (EmbedUtilsTests, test_all_string_permutations){
  Redirect a(true, true);
  vector<string> correct_kmers = {"CC", "CS", "SC", "SS"};
  string cs = "CS";
  int length = 2;
  vector<string> kmers = all_string_permutations(cs, length);

  for(int i=0; i < 4; ++i){
    EXPECT_EQ(kmers[i], correct_kmers[i]);
  }
  string cc = "CC";

  kmers = all_string_permutations(cc, length);
  correct_kmers = {"CC"};
  for(int i=0; i < 1; ++i){
    EXPECT_EQ(kmers[i], correct_kmers[i]);
  }

}

TEST (EmbedUtilsTests, test_remove_duplicate_characters){
  Redirect a(true, true);
  string test("ATGCCCCCC");
  string output = remove_duplicate_characters(test);
  EXPECT_EQ(output, "ATGC");
}

TEST (EmbedUtilsTests, test_is_character_in_string){
  Redirect a(true, true);
  string query = "atg";
  string target = "qwera";
  string target2 = "qwert";
  string target3 = "qwge";
  string target4 = "qw";
  string empty = "";
  EXPECT_TRUE(are_characters_in_string(query, target));
  EXPECT_TRUE(are_characters_in_string(query, target2));
  EXPECT_TRUE(are_characters_in_string(query, target3));
  EXPECT_FALSE(are_characters_in_string(query, target4));
  EXPECT_FALSE(are_characters_in_string(empty, target4));
  EXPECT_FALSE(are_characters_in_string(empty, empty));
}

TEST (EmbedUtilsTests, test_convert_to_float) {
  Redirect a(true, true);
  string number = "121.212";
  float something = string_to_float(number);
  EXPECT_FLOAT_EQ(121.212, something);
  number = "-121.212";
  something = string_to_float(number);
  EXPECT_FLOAT_EQ(-121.212, something);
}

TEST (EmbedUtilsTests, test_convert_to_int) {
  Redirect a(true, true);
  string number = "121";
  int64_t something = string_to_int(number);
  EXPECT_FLOAT_EQ(121, something);
  number = "-121";
  something = string_to_int(number);
  EXPECT_FLOAT_EQ(-121, something);
}

TEST (EmbedUtilsTests, test_list_files_in_dir) {
  Redirect a(true, true);
  int counter = 0;
  string ext = ".tsv";
  for (auto &i: list_files_in_dir(ASSIGNMENT_DIR, ext)){
    counter += 1;
  }
  EXPECT_EQ(1, counter);
  counter = 0;
  ext = "";
  for (auto &i: list_files_in_dir(ASSIGNMENT_DIR, ext)){
    counter += 1;
  }
  EXPECT_EQ(1, counter);
  counter = 0;
  ext = "fake";
  for (auto &i: list_files_in_dir(ASSIGNMENT_DIR, ext)){
    counter += 1;
  }
  EXPECT_EQ(0, counter);

  counter = 0;
  ext = "";
  dir_coro::pull_type data = list_files_in_dir(ALIGNMENT_DIR, ext);
  while (data){
    path something = data.get();
    counter += 1;
    data();
  }
  EXPECT_EQ(5, counter);
}

TEST (EmbedUtilsTests, test_create_ambig_bases) {
  Redirect a(true, true);
  std::map<string, string> ambig_bases = create_ambig_bases();
  EXPECT_EQ("AF", ambig_bases["f"]);
}

TEST (EmbedUtilsTests, test_create_ambig_bases2) {
  Redirect a(true, true);
  std::map<string, string> ambig_bases = create_ambig_bases2(AMBIG_MOD_FILE.string());
  EXPECT_EQ("Aa", ambig_bases["B"]);
  ambig_bases = create_ambig_bases2("");
  EXPECT_EQ("AF", ambig_bases["f"]);
  ASSERT_THROW(create_ambig_bases2("asdf"), AssertionFailureException);
}


void test_test(int the){
  the += 1;
}

TEST (EmbedUtilsTests, test_get_time) {
  Redirect a(true, true);
  auto funct = bind(test_test, 10);
  string the_time = get_time_string(funct);
  tuple<uint64_t, uint64_t, uint64_t, uint64_t> the_time2 = get_time(funct);
  EXPECT_EQ(0, get<0>(the_time2));
  EXPECT_EQ(0, get<1>(the_time2));
  EXPECT_EQ(0, get<2>(the_time2));
  EXPECT_EQ("hours: 0 minutes: 0", the_time.substr(0, 19));
}



#endif //EMBED_FAST5_TESTS_SRC_EMBEDUTILSTESTS_HPP_
