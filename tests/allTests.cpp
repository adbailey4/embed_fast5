//
// Created by Andrew Bailey on 03/15/19.
//

#include <gtest/gtest.h>
#include "fast5.hpp"
#include "iostream"
#include "nanopolish_squiggle_read.h"
#include "embed_fast5.hpp"
#include "nanopolish_read_db.h"
#include <boost/filesystem.hpp>

#define ORIGINAL_FAST5 "../tests/test_files/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_mux_scan_MA_821_R9_4_NA12878_11_17_16_95723_ch458_read26_strand.fast5"
#define EMPTY_FAST5 "../tests/test_files/empty_tester.fast5"
#define SIGNAL_FAST5 "../tests/test_files/just_signal.fast5"
#define NO_FAST5 "../tests/test_files/non_existent.fast5"
#define NO_EVENT "../tests/test_files/new_individual_files/read_002f9702-c19e-48c2-8e72-9021adbd4a48.fast5"

using namespace boost::filesystem;
using namespace std;

TEST (Fast5AccessTest, isValidFile) {
    EXPECT_TRUE(fast5::File::is_valid_file(ORIGINAL_FAST5));
    EXPECT_TRUE(fast5::File::is_valid_file(SIGNAL_FAST5));
    EXPECT_TRUE(fast5::File::is_valid_file(EMPTY_FAST5));
}

TEST (Fast5AccessTest, hasRequiredFields) {
//    Original Fast5 outputs
    fast5::File f;
    f.open(ORIGINAL_FAST5);
    EXPECT_TRUE(f.have_channel_id_params());
    EXPECT_TRUE(f.have_sampling_rate());
    EXPECT_TRUE(f.have_tracking_id_params());
    EXPECT_TRUE(f.have_context_tags_params());
    EXPECT_TRUE(f.have_eventdetection_events());
//    just raw and UniqueGlobalKey
    fast5::File f2;
    f2.open(SIGNAL_FAST5);
    EXPECT_TRUE(f2.have_channel_id_params());
    EXPECT_TRUE(f2.have_sampling_rate());
    EXPECT_TRUE(f2.have_tracking_id_params());
    EXPECT_TRUE(f2.have_context_tags_params());
    EXPECT_FALSE(f2.have_eventdetection_events());

//    Empty fast5 file
    fast5::File f3;
    f3.open(EMPTY_FAST5);
    EXPECT_FALSE(f3.have_channel_id_params());
    EXPECT_FALSE(f3.have_sampling_rate());
    EXPECT_FALSE(f3.have_tracking_id_params());
    EXPECT_FALSE(f3.have_context_tags_params());
    EXPECT_FALSE(f3.have_eventdetection_events());

}

TEST (Fast5AccessTest, test_copyFile){
    EXPECT_TRUE(is_regular_file(EMPTY_FAST5));
    EXPECT_TRUE(fast5::File::is_valid_file(EMPTY_FAST5));

    copy_file(EMPTY_FAST5, NO_FAST5);
    EXPECT_TRUE(is_regular_file(NO_FAST5));
    EXPECT_TRUE(fast5::File::is_valid_file(NO_FAST5));

    remove(NO_FAST5);
    EXPECT_FALSE(exists(NO_FAST5));

}

TEST (Fast5AccessTest, test_addChannelParams) {
//    crete test file
    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5, true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5);

    auto channel_id_params = original_f.get_channel_id_params();
    emtpy_f.add_channel_id_params(channel_id_params);

    EXPECT_TRUE(emtpy_f.have_channel_id_params());
    EXPECT_TRUE(emtpy_f.have_sampling_rate());
//    remove test file
    remove(NO_FAST5);
}

TEST (Fast5WriteTests, test_addBasecalledGroup) {
//    crete test file
    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5, true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5);

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
    remove(NO_FAST5);
}

TEST (Fast5WriteTests, test_event_table_to_basecalled_table) {
//    crete test file
    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5, true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5);

    auto basecall_groups = original_f.get_basecall_group_list();
    unsigned basecall_group_number = 0;

    string read_db_path("../tests/test_files/new_individual_files/new_individual.fastq");
    string test_read(NO_EVENT);
    string read_id("002f9702-c19e-48c2-8e72-9021adbd4a48");


    ReadDB read_db;
    read_db.load(read_db_path);

    SquiggleRead sr(read_id, read_db);
    auto basecall_table = event_table_to_basecalled_table(sr.events[0], 10);

    emtpy_f.add_basecall_events(basecall_group_number, basecall_groups[0], basecall_table);

    EXPECT_TRUE(emtpy_f.have_basecall_group());
    EXPECT_EQ(original_f.get_basecall_group_list()[0], basecall_groups[0]);
//    remove test file
    remove(NO_FAST5);
}

TEST (Fast5WriteTests, test_generate_basecall_table) {
//    crete test file
    copy_file(EMPTY_FAST5, NO_FAST5);

    fast5::File emtpy_f;
    emtpy_f.open(NO_FAST5, true);

    fast5::File original_f;
    original_f.open(ORIGINAL_FAST5);

    auto basecall_groups = original_f.get_basecall_group_list();
    unsigned basecall_group_number = 0;

    string read_db_path("../tests/test_files/new_individual_files/new_individual.fastq");
    string test_read(NO_EVENT);
    string read_id("002f9702-c19e-48c2-8e72-9021adbd4a48");


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
//    crete test file
    copy_file(NO_EVENT, NO_FAST5);

    const string read_db_path("../tests/test_files/test_readdb/new_individual.fastq");
    string test_read(NO_EVENT);
    string read_id("002f9702-c19e-48c2-8e72-9021adbd4a48");

    ReadDB read_db;
    read_db.load(read_db_path);

    embed_single_read(read_db, read_id, NO_FAST5);

    remove(NO_FAST5);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
