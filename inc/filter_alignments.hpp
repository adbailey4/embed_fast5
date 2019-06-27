//
// Created by Andrew Bailey on 2019-06-07.
//

#ifndef EMBED_FAST5_FILTER_ALIGNMENTS_H
#define EMBED_FAST5_FILTER_ALIGNMENTS_H


#include <map>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/filesystem.hpp>


using namespace std;
using namespace boost;
using namespace boost::icl;
using namespace boost::filesystem;


class PositionsFile
{
public:
    PositionsFile();
    PositionsFile(const string& input_reads_filename, int64_t k);
    ~PositionsFile();

    void load(const string& input_reads_filename, int64_t k);
    bool is_in(string& contig, int64_t position);
    //
    std::map<string, interval_set<int64_t>> m_data;

};


class AlignmentFile
{
public:
    explicit AlignmentFile(const string& input_reads_filename);
    ~AlignmentFile();
    string get_strand();
    int64_t get_k();
    void filter(PositionsFile* pf, boost::filesystem::path& output_file, string bases);
    //
    string strand;
    string file_path;
    int64_t k;
};

int filter_alignments_main(int argc, char** argv);
void filter_alignment_files(string input_reads_dir, const string& positions_file, string output_dir, string bases);

#endif //EMBED_FAST5_FILTER_ALIGNMENTS_H
