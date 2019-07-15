//
// Created by Andrew Bailey on 03/15/19.
//

#include "EmbedUtils.hpp"
#include "MaxKmers.hpp"
#include "AlignmentFile.hpp"
#include "scripts/top_kmers.hpp"
#include "scripts/FilterAlignments.hpp"
#include "fast5.hpp"
#include "iostream"
#include "nanopolish_squiggle_read.h"
#include "scripts/embed_fast5.hpp"
#include "nanopolish_read_db.h"
#include "omp.h"
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

using namespace boost::filesystem;
using namespace std;
using namespace embed_utils;


#define ORIGINAL_FAST51 "/tests/test_files/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_mux_scan_MA_821_R9_4_NA12878_11_17_16_95723_ch458_read26_strand.fast5"
#define EMPTY_FAST51 "/tests/test_files/empty_tester.fast5"
#define SIGNAL_FAST51 "/tests/test_files/just_signal.fast5"
#define NO_FAST51 "/tests/test_files/non_existent.fast5"
#define NO_EVENT1 "/tests/test_files/new_individual_files/read_002f9702-c19e-48c2-8e72-9021adbd4a48.fast5"
#define READ_DB1 "/tests/test_files/new_individual_files/new_individual.fastq"
#define R94_FAST51 "tests/test_files/r94_tests/fast5"
#define R94_FASTQ1 "tests/test_files/r94_tests/small_sample.fastq"
#define R94_TEST_DIR1 "tests/test_files/r94_tests"
#define POSITIONS_FILE1 "tests/test_files/positions_tests/CCWGG_ecoli_k12_mg1655.positions"
#define ALIGNMENT_FILE1 "tests/test_files/positions_tests/5cc86bac-79fd-4897-8631-8f1c55954a45.sm.backward.tsv"
#define ALIGNMENT_DIR1 "tests/test_files/positions_tests/"
#define CORRECT_OUTPUT1 "tests/test_files/positions_tests/correct_outputs/"
#define ASSIGNMENT_FILE1 "tests/test_files/assignment_files/d6160b0b-a35e-43b5-947f-adaa1abade28.sm.assignments.tsv"
#define TEST_FILES1 "tests/test_files/"

path HOME = "This is not a path";
path ORIGINAL_FAST5 = ORIGINAL_FAST51;
path EMPTY_FAST5 = EMPTY_FAST51;
path SIGNAL_FAST5 = SIGNAL_FAST51;
path NO_EVENT = NO_EVENT1;
path NO_FAST5 = NO_FAST51;
path READ_DB = READ_DB1;
path R94_FAST5 = R94_FAST51;
path R94_FASTQ = R94_FASTQ1;
path R94_TEST_DIR = R94_TEST_DIR1;
path POSITIONS_FILE = POSITIONS_FILE1;
path ALIGNMENT_FILE = ALIGNMENT_FILE1;
path ALIGNMENT_DIR = ALIGNMENT_DIR1;
path CORRECT_OUTPUT = CORRECT_OUTPUT1;
path ASSIGNMENT_FILE = ASSIGNMENT_FILE1;
path TEST_FILES = TEST_FILES1;


TEST (Fast5AccessTest, isValidFile) {
    EXPECT_TRUE(fast5::File::is_valid_file(ORIGINAL_FAST5.string()));
    EXPECT_TRUE(fast5::File::is_valid_file(SIGNAL_FAST5.string()));
    EXPECT_TRUE(fast5::File::is_valid_file(EMPTY_FAST5.string()));
}

TEST (Fast5AccessTest, hasRequiredFields) {
//    Original Fast5 outputs
    fast5::File f;
    f.open(ORIGINAL_FAST5.string());
    EXPECT_TRUE(f.have_channel_id_params());
    EXPECT_TRUE(f.have_sampling_rate());
    EXPECT_TRUE(f.have_tracking_id_params());
    EXPECT_TRUE(f.have_context_tags_params());
    EXPECT_TRUE(f.have_eventdetection_events());
//    just raw and UniqueGlobalKey
    fast5::File f2;
    f2.open(SIGNAL_FAST5.string());
    EXPECT_TRUE(f2.have_channel_id_params());
    EXPECT_TRUE(f2.have_sampling_rate());
    EXPECT_TRUE(f2.have_tracking_id_params());
    EXPECT_TRUE(f2.have_context_tags_params());
    EXPECT_FALSE(f2.have_eventdetection_events());

//    Empty fast5 file
    fast5::File f3;
    f3.open(EMPTY_FAST5.string());
    EXPECT_FALSE(f3.have_channel_id_params());
    EXPECT_FALSE(f3.have_sampling_rate());
    EXPECT_FALSE(f3.have_tracking_id_params());
    EXPECT_FALSE(f3.have_context_tags_params());
    EXPECT_FALSE(f3.have_eventdetection_events());

}

TEST (Fast5AccessTest, test_copyFile){
    EXPECT_TRUE(is_regular_file(EMPTY_FAST5));
    EXPECT_TRUE(fast5::File::is_valid_file(EMPTY_FAST5.string()));
    if (exists(NO_FAST5)){
        remove(NO_FAST5);
    }

    copy_file(EMPTY_FAST5, NO_FAST5);
    EXPECT_TRUE(is_regular_file(NO_FAST5));
    EXPECT_TRUE(fast5::File::is_valid_file(NO_FAST5.string()));

    remove(NO_FAST5);
    EXPECT_FALSE(exists(NO_FAST5));

}

TEST (Fast5AccessTest, test_addChannelParams) {
//    crete test file
    if (exists(NO_FAST5)){
        remove(NO_FAST5);
    }

    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5.string(), true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5.string());

    auto channel_id_params = original_f.get_channel_id_params();
    emtpy_f.add_channel_id_params(channel_id_params);

    EXPECT_TRUE(emtpy_f.have_channel_id_params());
    EXPECT_TRUE(emtpy_f.have_sampling_rate());
//    remove test file
    remove(NO_FAST5);
}

TEST (Fast5WriteTests, test_addBasecalledGroup) {
//    crete test file
    if (exists(NO_FAST5)){
        remove(NO_FAST5);
    }

    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5.string(), true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5.string());

    bool have_basecall_group = original_f.have_basecall_group();
    EXPECT_TRUE(have_basecall_group);

    auto basecall_groups = original_f.get_basecall_group_list();
    auto basecall_events = original_f.get_basecall_events(0);
    unsigned basecall_group_number = 0;

    emtpy_f.add_basecall_events(basecall_group_number, basecall_groups[0], basecall_events);

    auto fastq = original_f.get_basecall_fastq(0);
    emtpy_f.add_basecall_fastq(basecall_group_number, basecall_groups[0], fastq);

    EXPECT_TRUE(emtpy_f.have_basecall_group());
    EXPECT_EQ(original_f.get_basecall_group_list()[0], basecall_groups[0]);
//    remove test file
    emtpy_f.close();
    original_f.close();
    remove(NO_FAST5);
}

TEST (Fast5WriteTests, test_event_table_to_basecalled_table) {
//    crete test file
    if (exists(NO_FAST5)){
        remove(NO_FAST5);
    }
    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5.string(), true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5.string());

    auto basecall_groups = original_f.get_basecall_group_list();
    unsigned basecall_group_number = 0;

    const string& read_db_path(READ_DB.string());
//    const string& test_read(NO_EVENT.string());
    string read_id("002f9702-c19e-48c2-8e72-9021adbd4a48");


    ReadDB read_db;
    read_db.load(read_db_path);

    SquiggleRead sr(read_id, read_db);
    auto basecall_table = event_table_to_basecalled_table(sr.events[0], 10);

    emtpy_f.add_basecall_events(basecall_group_number, basecall_groups[0], basecall_table);

    EXPECT_TRUE(emtpy_f.have_basecall_group());
    EXPECT_EQ(original_f.get_basecall_group_list()[0], basecall_groups[0]);
//    remove test file
    emtpy_f.close();
    original_f.close();
//    remove(NO_FAST5);
}

TEST (Fast5WriteTests, test_generate_basecall_table) {
//    crete test file
    if (exists(NO_FAST5)){
        remove(NO_FAST5);
    }

    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5.string(), true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5.string());

    auto basecall_groups = original_f.get_basecall_group_list();
    unsigned basecall_group_number = 0;

    string read_db_path((const std::string) READ_DB.string());
    string test_read((const std::string) NO_EVENT.string());
    string read_id("002f9702-c19e-48c2-8e72-9021adbd4a48");


    ReadDB read_db;
    read_db.load(read_db_path);

    SquiggleRead sr(read_id, read_db);
    auto basecall_table = generate_basecall_table(sr);

    emtpy_f.add_basecall_events(basecall_group_number, basecall_groups[0], basecall_table);

    EXPECT_TRUE(emtpy_f.have_basecall_group());
    EXPECT_EQ(original_f.get_basecall_group_list()[0], basecall_groups[0]);

//    remove(NO_FAST5);
}

TEST (Fast5EmbedTests, test_embed_using_readdb) {
//    crete test file
    if (exists(NO_FAST5)){
        remove(NO_FAST5);
    }

    copy_file(NO_EVENT, NO_FAST5);

    const string& read_db_path(READ_DB.string());
    string read_id("002f9702-c19e-48c2-8e72-9021adbd4a48");

    ReadDB read_db;
    read_db.load(read_db_path);

    embed_single_read(read_db, read_id, NO_FAST5.string());

    remove(NO_FAST5);
}

TEST (Fast5EmbedTests, test_r94_embed) {
//    crete test file

    path tempdir = temp_directory_path() / "temp";
    path fast5_dir = tempdir / "fast5";

//    cout << tempdir << "\n";
    path read_db_path = tempdir / R94_FASTQ.filename();
    if (exists(tempdir)){
        remove_all(tempdir);
    }
    copyDir(R94_TEST_DIR, tempdir);
    ReadDB read_db;
    read_db.load(read_db_path.string());
    cd(fast5_dir.string().c_str());

    embed_using_readdb(read_db_path.string(), read_db);

    path some_file = "DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_89607_ch108_read1153_strand.fast5";
    fast5::File original_f;
    original_f.open((fast5_dir / some_file).string());
    EXPECT_EQ(original_f.get_basecall_group_list()[0], "1D_000");

    remove_all(tempdir);

}

TEST (Fast5EmbedTests, test_multiprocess_embed_using_readdb){

    path tempdir = temp_directory_path() / "temp";
    path fast5_dir = tempdir / "fast5";

//    cout << tempdir << "\n";
    path read_db_path = tempdir / R94_FASTQ.filename();
    if (exists(tempdir)){
        remove_all(tempdir);
    }
    copyDir(R94_TEST_DIR, tempdir);
    ReadDB read_db;
    read_db.load(read_db_path.string());
    cd(fast5_dir.string().c_str());

    int num_threads = 5;
    omp_set_num_threads(num_threads);

    #ifndef H5_HAVE_THREADSAFE
        if(num_threads > 1) {
            fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
            fprintf(stderr, "Please recompile libhdf5 or run with -t 1\n");
            exit(1);
        }
    #endif


    multiprocess_embed_using_readdb(read_db_path.string(), read_db);

    path some_file = "DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_89607_ch108_read1153_strand.fast5";
    fast5::File original_f;
    original_f.open((fast5_dir / some_file).string());
    EXPECT_EQ(original_f.get_basecall_group_list()[0], "1D_000");

    remove_all(tempdir);

}

TEST (PositionsFileTests, test_load) {
    PositionsFile pf = PositionsFile(POSITIONS_FILE.string(), 5);
    string contig = "gi_ecoli+";
    EXPECT_TRUE(contains(pf.m_data[contig], 419));
}

TEST (PositionsFileTests, test_is_in) {
    PositionsFile pf = PositionsFile(POSITIONS_FILE.string(), 5);
    string contig = "gi_ecoli+";
    EXPECT_TRUE(pf.is_in(contig, 419));
    EXPECT_TRUE(pf.is_in(contig, 415));
    EXPECT_FALSE(pf.is_in(contig, 414));
}

TEST (AlignmentFileTests, test_get_strand) {
    AlignmentFile af = AlignmentFile(ALIGNMENT_FILE.string());
    af.get_strand();
    EXPECT_EQ("-", af.strand);
}



TEST (AlignmentFileTests, test_filter) {
    path tempdir = temp_directory_path() / "temp";
    path test_output = tempdir / "test_output.assignment.tsv";
    AlignmentFile af = AlignmentFile(ALIGNMENT_FILE.string());
    PositionsFile pf = PositionsFile(POSITIONS_FILE.string(), 5);
    af.filter(&pf, test_output, "");
    std::ifstream in_file(test_output.c_str());
    if (in_file.good()) {
        // read the file
        std::string line;
        while (getline(in_file, line)) {
            std::vector<std::string> fields = embed_utils::split(line, '\t');
            std::size_t found = fields[0].find('C');
            EXPECT_TRUE(found!=std::string::npos);
        }
    }
    in_file.close();
    remove_all(tempdir);
}

TEST (filter_alignments, test_filter_alignment_files) {
  path input_dir = temp_directory_path() / "input" ;
  path output_dir = temp_directory_path() / "output" ;
  copyDir(ALIGNMENT_DIR, input_dir);
  string empty;
//    omp_set_num_threads(2); // Use 2 threads
  filter_alignment_files(input_dir.string(), POSITIONS_FILE.string(), output_dir.string(), empty);

  directory_iterator end_itr;
//    Get all tsvs to process
  for (directory_iterator itr(output_dir); itr != end_itr; ++itr) {

      EXPECT_TRUE(compareFiles(itr->path().string(), (CORRECT_OUTPUT / itr->path().filename()).string()));

  }
  remove_all(input_dir);
  remove_all(output_dir);
}


TEST (AssignmentFileTests, test_iterate_assignment) {
    path input_dir = temp_directory_path() / "input" ;
    path output_dir = temp_directory_path() / "output" ;
    AssignmentFile af = AssignmentFile(ASSIGNMENT_FILE.string());
    float answer = 83.7093;
    for (auto &x: af.iterate()){
      EXPECT_FLOAT_EQ(x.mean, answer);
      break;
    }
}

TEST (AssignmentFileTests, test_get_k) {
  path input_dir = temp_directory_path() / "input" ;
  path output_dir = temp_directory_path() / "output" ;
  AssignmentFile af = AssignmentFile(ASSIGNMENT_FILE.string());
  int64_t answer = af.get_k();
  EXPECT_EQ(6, answer);
}



TEST (MaxKmersTests, test_create_kmers) {
  MaxKmers mk = MaxKmers(10, "ATGC", 5);
  string alphabet = "ATGC";
  vector<string> kmers = mk.create_kmers(alphabet, 5);
  EXPECT_EQ(kmers.size(), 1024);
  kmers = mk.create_kmers(alphabet, 6);
  EXPECT_EQ(kmers.size(), 4096);

}

TEST (MaxKmersTests, test_get_kmer_index) {
  MaxKmers mk = MaxKmers(10, "ATGC", 5);
  string poly_a = "AAAAA";
  int poly_a_index = mk.get_kmer_index(poly_a);
  EXPECT_EQ(poly_a_index, 0);
  string poly_t = "TTTTT";
  int poly_t_index = mk.get_kmer_index(poly_t);
  EXPECT_EQ(poly_t_index, 1023);
}

TEST (MaxKmersTests, test_get_index_kmer) {
  MaxKmers mk = MaxKmers(10, "ATGC", 5);
  string poly_a = mk.get_index_kmer(0);
  EXPECT_EQ(poly_a, "AAAAA");
  string poly_t = mk.get_index_kmer(1023);
  EXPECT_EQ(poly_t, "TTTTT");
}

TEST (MaxKmersTests, test_add_to_heap) {
  MaxKmers mk = MaxKmers(10, "ATGC", 5);
  eventkmer my_kmer = eventkmer("AAAAA", 10, "t", 1.2);
  mk.add_to_heap(my_kmer);
  eventkmer my_kmer2 = eventkmer("AAAAA", 10, "t", 1.4);
#pragma omp parallel for shared(mk, my_kmer2)
  for (int i = 0; i < 10; i++){
    mk.add_to_heap(my_kmer2);
  }
  eventkmer match = mk.kmer_queues[0].top();
  EXPECT_EQ(match.kmer, my_kmer2.kmer);
  EXPECT_EQ(match.strand, my_kmer2.strand);
  EXPECT_EQ(match.mean, my_kmer2.mean);
  EXPECT_EQ(match.prob, my_kmer2.prob);

}

TEST (MaxKmersTests, test_write_to_file) {
  MaxKmers mk = MaxKmers(10, "ATGC", 5);
  for (int i = 10; i > 0; i--){
    eventkmer my_kmer = eventkmer("AAAAA", 10, "t", i);
    mk.add_to_heap(my_kmer);
  }
  if (exists(NO_FAST5)){
    remove(NO_FAST5);
  }
  mk.write_to_file(NO_FAST5);
  EXPECT_TRUE(is_regular_file(NO_FAST5));
  remove(NO_FAST5);

  path log_file = TEST_FILES / "log_file.txt";

  mk.write_to_file(NO_FAST5, log_file);
  EXPECT_TRUE(is_regular_file(NO_FAST5));
  EXPECT_TRUE(is_regular_file(log_file));
  remove(NO_FAST5);
  remove(log_file);


}


TEST (util_functions, test_sort_string){
    string something = "something";
    sort_string(something);
    EXPECT_TRUE(something.compare("eghimnost")==0);
}

TEST (util_functions , test_all_lexicographic_recur){
  string alphabet = "AC";
  int length = 2;
  vector<string> kmers = all_lexicographic_recur(alphabet, "  ", length-1, 0);
  EXPECT_EQ(kmers.size(), 4);

}


TEST (util_functions, test_all_string_permutations){
  vector<string> correct_kmers = {"CC", "CS", "SC", "SS"};
  vector<string> kmers = all_string_permutations("CS", 2);

  for(int i=0; i < 4; ++i){
    EXPECT_EQ(kmers[i], correct_kmers[i]);
  }
  kmers = all_string_permutations("CC", 2);
  correct_kmers = {"CC"};
  for(int i=0; i < 1; ++i){
    EXPECT_EQ(kmers[i], correct_kmers[i]);
  }

}

TEST (util_functions, test_remove_duplicate_characters){
  string output = remove_duplicate_characters("ATGCCCCCC");
  EXPECT_EQ(output, "ATGC");
}

TEST (util_functions, test_is_character_in_string){
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

TEST (util_functions, test_convert_to_float) {
    string number = "121.212";
    string not_number = "asdf";
    float something = convert_to_float(number);
    EXPECT_FLOAT_EQ(121.212, something);
}


int main(int argc, char **argv) {
  H5Eset_auto(0, nullptr, nullptr);
  HOME = argv[1];
  cout << HOME << '\n';
  ORIGINAL_FAST5 = HOME / ORIGINAL_FAST5;
  EMPTY_FAST5 = HOME / EMPTY_FAST5;
  SIGNAL_FAST5 = HOME / SIGNAL_FAST5;
  NO_EVENT = HOME / NO_EVENT;
  NO_FAST5 = HOME / NO_FAST5;
  READ_DB = HOME / READ_DB;
  R94_FAST5 = HOME / R94_FAST5;
  R94_FASTQ = HOME / R94_FASTQ;
  R94_TEST_DIR = HOME / R94_TEST_DIR;
  POSITIONS_FILE = HOME / POSITIONS_FILE;
  ALIGNMENT_FILE = HOME / ALIGNMENT_FILE;
  ALIGNMENT_DIR = HOME / ALIGNMENT_DIR;
  CORRECT_OUTPUT = HOME / CORRECT_OUTPUT;
  ASSIGNMENT_FILE = HOME / ASSIGNMENT_FILE;
  TEST_FILES = HOME / TEST_FILES;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

