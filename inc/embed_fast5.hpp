//
// Created by Andrew Bailey on 03/18/19.
//

#include "fast5.hpp"
#include "nanopolish_squiggle_read.h"


std::vector< fast5::Basecall_Event > event_table_to_basecalled_table(std::vector< SquiggleEvent > et);

std::vector< fast5::Basecall_Event > generate_basecall_table(SquiggleRead& read);

int embed_fast5_main(int argc, char** argv);

template<class T, size_t N, class V>
std::array<T, N> to_array(const V& v)
{
    assert(v.size() <= N);
    std::array<T, N> d;
    using std::begin; using std::end;
    std::copy( begin(v), end(v), begin(d) ); // this is the recommended way
    return d;
}
