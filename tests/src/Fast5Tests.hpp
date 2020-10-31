//
// Created by Andrew Bailey on 5/2/20.
//

#ifndef EMBED_FAST5_TESTS_SRC_FAST5TESTS_HPP_
#define EMBED_FAST5_TESTS_SRC_FAST5TESTS_HPP_

// embed lib
#include "EmbedUtils.hpp"
#include "EmbedFast5.hpp"
// embed test files
#include "TestFiles.hpp"
// boost
#include <boost/filesystem.hpp>
// gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>
// openMP
#include <omp.h>

using namespace test_files;
using namespace embed_utils;
using namespace boost::filesystem;


TEST (Fast5AccessTest, isValidFile) {
  Redirect a(true, true);
  EXPECT_TRUE(fast5::File::is_valid_file(ORIGINAL_FAST5.string()));
  EXPECT_TRUE(fast5::File::is_valid_file(SIGNAL_FAST5.string()));
  EXPECT_TRUE(fast5::File::is_valid_file(EMPTY_FAST5.string()));
}

TEST (Fast5AccessTest, hasRequiredFields) {
  //    Original Fast5 outputs
  Redirect a(true, true);
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
  Redirect a(true, true);
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
  Redirect a(true, true);
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
  Redirect a(true, true);
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
  Redirect a(true, true);
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
  cd(READ_DB_DIR.string().c_str());

  SquiggleRead sr(read_id, read_db);
  auto basecall_table = event_table_to_basecalled_table(sr.events[0], 10);

  emtpy_f.add_basecall_events(basecall_group_number, basecall_groups[0], basecall_table);

  EXPECT_TRUE(emtpy_f.have_basecall_group());
  EXPECT_EQ(original_f.get_basecall_group_list()[0], basecall_groups[0]);
  //    remove test file
  emtpy_f.close();
  original_f.close();
}

TEST (Fast5WriteTests, test_generate_basecall_table) {
  Redirect a(true, true);
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

  cd(READ_DB_DIR.string().c_str());

  ReadDB read_db;
  read_db.load(read_db_path);

  SquiggleRead sr(read_id, read_db);
  auto basecall_table = generate_basecall_table(sr);

  emtpy_f.add_basecall_events(basecall_group_number, basecall_groups[0], basecall_table);

  EXPECT_TRUE(emtpy_f.have_basecall_group());
  EXPECT_EQ(original_f.get_basecall_group_list()[0], basecall_groups[0]);

  remove(NO_FAST5);
}

TEST (Fast5EmbedTests, test_embed_using_readdb) {
  Redirect a(true, true);
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
  Redirect a(true, true);
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
  Redirect a(true, true);
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



#endif //EMBED_FAST5_TESTS_SRC_FAST5TESTS_HPP_
