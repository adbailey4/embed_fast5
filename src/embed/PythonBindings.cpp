//
// Created by Andrew Bailey on 2019-06-06.
//

#ifndef EMBED_PYBIND_API_H
#define EMBED_PYBIND_API_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "LoadVariantPaths.hpp"

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
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: cmake_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

  module.def("add", &add1, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

  module.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  module.attr("__version__") = "dev";
#endif
}

#endif
