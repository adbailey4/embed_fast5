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

int add1(int i, int j) {
    return i + j;
}

using namespace pybind11;


PYBIND11_MODULE(bindings, module) {

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

  module.doc() = R"pbdoc(
        Embed Wrappers:
        - add
        - subtract
        - LoadVariantPaths
        - generate_master_assignment_table
    )pbdoc";

  module.def("add", &add1, R"pbdoc(
        Add two numbers
    )pbdoc");

  module.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
    )pbdoc");

  module.def("generate_master_assignment_table", &generate_master_assignment_table, R"pbdoc(
  - Generate the master assignment table by parsing assignment files and outputting the top n kmers to a files
 @param assignment_dir: path to assignment files directory
 @param output_dir: path to output directory where new builtAssignment.tsv will be written
 @param heap_size: number of max kmers to keep for each kmer
 @param alphabet: alphabet used to generate kmers

    )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  module.attr("__version__") = "dev";
#endif
}

#endif
