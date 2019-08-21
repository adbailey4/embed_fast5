//
// Created by Andrew Bailey on 2019-06-06.
//

#ifndef EMBED_PYBIND_API_H
#define EMBED_PYBIND_API_H

#include <pybind11/pybind11.h>
//#include "python_scripts/EmbedFast5.hpp"

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;


PYBIND11_MODULE(bindings, m) {
  m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: cmake_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

  m.def("add", &add, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

  m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}

#endif
