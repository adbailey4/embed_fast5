## pybind
pybind11_add_module(bindings PythonBindings.cpp )
if(BUILD_DYNAMIC_LIBRARY)
    target_link_libraries(bindings PRIVATE embedDynamicLibrary)
else()
    target_link_libraries(bindings PRIVATE embedStaticLibrary)
endif()
