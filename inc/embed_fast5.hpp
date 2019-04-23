//
// Created by Andrew Bailey on 03/18/19.
//

#include "fast5.hpp"
#include "nanopolish_squiggle_read.h"

#ifdef _WIN32
#include <direct.h>
// MSDN recommends against using getcwd & chdir names
#define cwd _getcwd
#define cd _chdir
#else
#include "unistd.h"
#define cwd getcwd
#define cd chdir
#endif


std::vector< fast5::Basecall_Event > event_table_to_basecalled_table(std::vector< SquiggleEvent > et,
                                                                     double start_time);

std::vector< fast5::Basecall_Event > generate_basecall_table(SquiggleRead& read);

int embed_fast5_main(int argc, char** argv);


void embed_using_readdb(const std::string& input_reads_filename, const ReadDB& read_db);


void embed_single_read(const ReadDB& read_db, std::string read_id, std::string fast5_path);

void multiprocess_embed_using_readdb(const std::string& input_reads_filename, const ReadDB& read_db);

template<class T, size_t N, class V>
std::array<T, N> to_array(const V& v)
{
    assert(v.size() < N);
    std::array<T, N> d;
    using std::begin; using std::end;
    std::copy( begin(v), end(v), begin(d) ); // this is the recommended way
    d[v.size()] = '\0';
    return d;
}
