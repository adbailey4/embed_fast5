//
// Created by Andrew Bailey on 2019-06-06.
//

#include <pybind11/pybind11.h>
#include "scripts/embed_fast5.hpp"


namespace py = pybind11;


int add(int i, int j) {
    return i + j;
}


PYBIND11_MODULE(some_embed2, m) {
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

//    m.def("multiprocess_embed_using_readdb", &multiprocess_embed_using_readdb, "A function which embeds reads");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

//PYBIND11_MODULE(embed, m) {
//    m.doc() = "Embed reads"; // optional module docstring
//
//    m.def("multiprocess_embed_using_readdb", &multiprocess_embed_using_readdb, "A function which embeds reads");
//}
