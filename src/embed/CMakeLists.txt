## pybind
if(BUILD_SHARED_LIBS)
    message(STATUS "Building Python Bindings")
    pybind11_add_module(bindings PythonBindings.cpp)
    target_link_libraries(bindings PRIVATE embedlib)
else(NOT BUILD_SHARED_LIBS)
    message(STATUS "Skipping Python Bindings: Cannot build using static libs")
    pybind11_add_module(bindings PythonBindings.cpp)
    target_link_libraries(bindings PRIVATE embedlib)
endif(BUILD_SHARED_LIBS)

