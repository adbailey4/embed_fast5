//
// Created by Andrew Bailey on 2019-06-06.
//

#ifndef EMBED_PYBIND_API_H
#define EMBED_PYBIND_API_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "LoadVariantPaths.hpp"
#include "TopKmers.hpp"

using namespace pybind11;


PYBIND11_MODULE(bindings, module) {
  module.doc() = R"pbdoc(
        Embed Wrappers:
        - LoadVariantPaths
        - generate_master_kmer_table
    )pbdoc";

  class_<LoadVariantPaths>(module, "LoadVariantPaths")
      .def(init<const std::string &, const std::string &, bool, uint64_t >(),
          "Load variants based on ",
          pybind11::arg("positions_file_path"),
          pybind11::arg("full_sa_dir"),
          pybind11::arg("rna") = false,
          pybind11::arg("num_locks") = 100)
      .def("write_per_read_calls",
          &LoadVariantPaths::write_per_read_calls,
           "Write read_id, path_id, path, and nucleotide path for each read",
          pybind11::arg("output_path"))
      .def("write_per_path_counts",
          &LoadVariantPaths::write_per_path_counts,
           "For each path_id/ path/ nucleotide path write the number of reads falling into that path",
           pybind11::arg("output_path"))
           ;

  module.def("generate_master_kmer_table", &generate_master_kmer_table_wrapper, R"pbdoc(
  - Generate the master assignment table by parsing assignment files and outputting the top n kmers to a files
 @param assignment_dir: path to assignment files directory
 @param output_dir: path to output directory where new builtAssignment.tsv will be written
 @param heap_size: number of max kmers to keep for each kmer
 @param alphabet: alphabet used to generate kmers
 @param n_threads: set number of threads to use: default 2
 @param verbose: print file names: default false

    )pbdoc",
    pybind11::arg("event_table_files"),
    pybind11::arg("output_dir"),
    pybind11::arg("heap_size"),
    pybind11::arg("alphabet"),
    pybind11::arg("min_prob") = 0.0,
    pybind11::arg("n_threads") = 2,
    pybind11::arg("verbose") = false);

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  module.attr("__version__") = "dev";
#endif
}

#endif
