//
// Created by Andrew Bailey on 03/17/19.
//

#include <getopt.h>
#include "fast5.hpp"
#include "embed_fast5.hpp"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_emissions.h"
#include "nanopolish_read_db.h"
#include "boost/filesystem.hpp"
#include <omp.h>


using namespace std;

std::vector< fast5::Basecall_Event > event_table_to_basecalled_table(std::vector< SquiggleEvent > et, double start_time){

    std::vector< fast5::Basecall_Event > basecalled_et(et.size());
    for (int i=0; i < et.size() ; i++){
        basecalled_et[i].start = et[i].start_time + start_time;
        basecalled_et[i].length = (double)  et[i].duration;
        basecalled_et[i].mean = (double) et[i].mean;
        basecalled_et[i].stdv = (double) et[i].stdv;
        basecalled_et[i].p_model_state = 0.0;
        basecalled_et[i].move = 0;

    }
    return basecalled_et;

}

std::vector< fast5::Basecall_Event > generate_basecall_table(SquiggleRead& read){

    std::vector<EventAlignment> alignment = read.get_eventalignment_for_1d_basecalls(read.read_sequence, "nucleotide",
                                                                                     read.base_to_event_map,
                                                                                     read.base_model[0]->k, 0, 0);
    const Alphabet* alphabet = read.base_model[0]->pmalphabet;
    //    generate empty basecall table to fill

    fast5::File f_p;
    f_p.open(read.fast5_path);
    auto& sample_read_names = f_p.get_raw_samples_read_name_list();
    if(sample_read_names.empty()) {
        fprintf(stderr, "Error, no raw samples found\n");
        exit(EXIT_FAILURE);
    }
    // we assume the first raw sample read is the one we're after
    std::string sample_read_name = sample_read_names.front();
    double sample_start_time = f_p.get_raw_samples_params(sample_read_name).start_time / f_p.get_sampling_rate();
    // convert event table to basecall table
    auto basecall_table = event_table_to_basecalled_table(read.events[0], sample_start_time);
    int prev_ref_index = alignment[0].ref_position;

    for (auto &i : alignment) {
        basecall_table[i.event_idx].model_state = to_array< char, MAX_K_LEN >(i.ref_kmer);
        uint32_t kmer_rank = alphabet->kmer_rank(i.ref_kmer.c_str(), read.base_model[0]->k);

        //  get prob from model
        float lp_emission = log_probability_match_r9(read, *read.base_model[0], kmer_rank,
                                                     (uint32_t) i.event_idx, (uint8_t) i.strand_idx);

        basecall_table[i.event_idx].p_model_state = exp(lp_emission);
        basecall_table[i.event_idx].move = i.ref_position - prev_ref_index;
        prev_ref_index = i.ref_position;
    }

    int start_event_index = alignment[0].event_idx;
    int end_event_index = alignment[alignment.size()-1].event_idx;

    basecall_table.erase (basecall_table.begin()+end_event_index, basecall_table.end());
    basecall_table.erase (basecall_table.begin(), basecall_table.begin() + start_event_index);

    return basecall_table;
}


void embed_single_read(const ReadDB& read_db, std::string read_id, std::string fast5_path){
    try
    {
        std::string fastq_sequence;
        //        make sure the file exists
        if (boost::filesystem::is_regular_file(fast5_path)){

            std::string read_sequence = read_db.get_read_sequence(read_id);
            //            make sure the read sequence is long enough to process
            if (read_sequence.length() > 10) {
                SquiggleRead sr(read_id, read_db);
                cout << sr.events[0].empty() << "\n";
                if (!sr.events[0].empty()) {
                    auto data = generate_basecall_table(sr);
                    fast5::File fast5_file;
                    fast5_file.open(fast5_path, true);
                    std::string path_1 = fast5::File::basecall_events_path("1D_000", 0);
                    if (!fast5_file.exists(path_1)) {
                        cout << "running: " << fast5_path << "\n";
                        fast5_file.add_basecall_events(0, "1D_000", data);
                    } else {
                        cout << "passing: " << fast5_path << "\n";

                    }
                    fast5_file.close();
                } else {
                    cout << "FAILED SR" << fast5_path << "\n";
                }
            } else{
                cout << "Too Short: " << fast5_path << "\n";

            }
        }
    }
    catch (int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
    }

}

//void embed_single_read(const ReadDB& read_db, std::string read_id, std::string fast5_path){
//    std::string read_sequence = read_db.get_read_sequence(read_id);
//    if (read_sequence.length() > 10){
//        SquiggleRead sr(read_id, read_db);
//        if (!sr.events[0].empty()){
//            auto data = generate_basecall_table(sr);
//
//            fast5::File fast5_file;
//            fast5_file.open(fast5_path, true);
//            auto basecall_groups = fast5_file.get_basecall_group_list();
//            std::string path_1 = fast5::File::basecall_events_path(basecall_groups[0], 0);
//            if (fast5_file.exists(path_1)) {
//                cout << path_1 << '\n';
//            }
//            else {
//                fast5_file.add_basecall_events(0, basecall_groups[0], data);
//            }
//
//        } else{
//            cout <<  '\n';
//
//        }
//    }
//}

void multiprocess_embed_using_readdb(const std::string& input_reads_filename, const ReadDB& read_db){
//        int N = 10;
//        float a[N], b[N], c[N];
//        int i;
//
//        /* Initialize arrays a and b */
//        for (i = 0; i < N; i++) {
//            a[i] = i * 2.0;
//            b[i] = i * 3.0;
//        }
//
//        /* Compute values of array c = a+b in parallel. */
//#pragma omp parallel shared(a, b, c) private(i)
//        {
//#pragma omp for
//            for (i = 0; i < N; i++) {
//                c[i] = a[i] + b[i];
//                printf ("%f\n", c[1]);
//            }
//        }


    omp_set_num_threads(2); // Use 4 threads for all consecutive parallel regions


        // generate input filenames
        std::string m_indexed_reads_filename = input_reads_filename + ".index";
        std::string in_filename = m_indexed_reads_filename + ".readdb";
        //
        std::ifstream in_file(in_filename.c_str());
        if(in_file.good()) {
            // read the database
            std::string line;
            std::vector<std::string> lines;

            // Read the file
            while(getline(in_file, line)) {
                lines.push_back(line);
            }
            in_file.close();
#pragma omp parallel for
            for(auto it = lines.begin(); it < lines.end(); it++) {

                std::vector<std::string> fields = split(*it, '\t');
                cout << "Hello, world." << *it << "\n";

                static std::string name = "";
                static std::string path = "";
                if (fields.size() == 2) {
                    name = fields[0];
                    path = fields[1];

                    embed_single_read(read_db, name, path);
                }
            }
    }

}

// basically a replica of ReadDB::load but I want access to the private data
void embed_using_readdb(const std::string& input_reads_filename, const ReadDB& read_db)
{

    // generate input filenames
    std::string m_indexed_reads_filename = input_reads_filename + ".index";
    std::string in_filename = m_indexed_reads_filename + ".readdb";
    //
    std::ifstream in_file(in_filename.c_str());
    if(in_file.good()) {
        // read the database
        std::string line;
        while(getline(in_file, line)) {
            std::vector<std::string> fields = split(line, '\t');

            static std::string name = "";
            static std::string path = "";
            if(fields.size() == 2) {
                name = fields[0];
                path = fields[1];
                embed_single_read(read_db, name, path);
            }
        }
    }
}



//
// Getopt
//
#define SUBPROGRAM "embed"
#define EMBED_VERSION "0.0.1"
#define THIS_NAME "run"
#define PACKAGE_BUGREPORT2 "None"

static const char *EMBED_FAST5_VERSION_MESSAGE =
        SUBPROGRAM " Version " EMBED_VERSION "\n";

static const char *EMBED_FAST5_USAGE_MESSAGE =
        "Usage: " THIS_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
        "Classify nucleotides as methylated or not.\n"
        "\n"
        "  -v, --verbose                        display verbose output\n"
        "      --version                        display version\n"
        "      --help                           display this help and exit\n"
        "  -r, --reads=FILE                     the ONT reads are in fasta FILE\n"
        "  -i, --read_id=ID                     read id of file to process\n"
        "  -f, --fast5=FILE                     path to Fast5 file\n"
        "\nReport bugs to " PACKAGE_BUGREPORT2 "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string read_id;
    static std::string fast5;
    static int num_threads = 1;
}

static const char* shortopts = "r:i:f:vn";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
        { "verbose",          no_argument,       nullptr, 'v' },
        { "reads",            required_argument, nullptr, 'r' },
        { "read_id",          optional_argument, nullptr, 'i' },
        { "fast5",            optional_argument, nullptr, 'f' },
        { "help",             no_argument,       nullptr, OPT_HELP },
        { "version",          no_argument,       nullptr, OPT_VERSION },
        { NULL, 0, NULL, 0 }
};

void parse_embed_main_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'i': arg >> opt::read_id; break;
            case 'f': arg >> opt::fast5; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << EMBED_FAST5_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << EMBED_FAST5_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::reads_file.empty()) {
        std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << EMBED_FAST5_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int embed_fast5_main(int argc, char** argv)
{
    parse_embed_main_options(argc, argv);
    ReadDB read_db;
    read_db.load(opt::reads_file);

//    cout << read_db.get_read_sequence("00087eac-80d1-4328-b755-24f4b192c7d0") << '\n';
//    cout << read_db.get_read_sequence("0037480f-176e-44b7-bd11-c23124837596") << '\n';
//    cout << read_db.get_read_sequence("0052794f-313f-4b93-9063-7fdb0f79d23c") << '\n';
//    cout << read_db.get_read_sequence("002f9702-c19e-48c2-8e72-9021adbd4a48") << '\n';
//    cout << read_db.get_read_sequence("002625e8-9a2e-411c-99a7-ee3fa6b9eef1") << '\n';
//    cout << read_db.get_read_sequence("002625e8-9a2e-411c-99a7-ee3fa6b9eef1") << '\n';

    embed_using_readdb(opt::reads_file, read_db);

//    SquiggleRead sr(opt::read_id, read_db);
//
//    auto data = generate_basecall_table(sr);
//
//    fast5::File fast5_file;
//    fast5_file.open(opt::fast5, true);
//    auto basecall_groups = fast5_file.get_basecall_group_list();
//    fast5_file.add_basecall_events(0, basecall_groups[0], data);
//
    return EXIT_SUCCESS;
}

// taskManager run -c '/home/ubuntu/embed_fast5/build/main_cpp embed -r ../rel3-nanopore-wgs-3574887596-FAB43577_2.fastq' && taskManager run -c 'python /home/ubuntu/modification_detection_pipeline/multi_fast5_pipeline/run_deep_mod_r9.py --config /home/ubuntu/modification_detection_pipeline/multi_fast5_pipeline/r9_pipeline.config.json' -r r94_ucsc_deepmod_1000_reads_04_11_19.txt
// taskManager run -c 'python /home/ubuntu/modification_detection_pipeline/multi_fast5_pipeline/run_deep_mod_r9.py --config /home/ubuntu/modification_detection_pipeline/multi_fast5_pipeline/r9_pipeline.config.json' --to andbaile@ucsc.edu -r