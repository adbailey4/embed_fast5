cmake_minimum_required(VERSION 3.7)

project(embedDynamicLibrary)

add_definitions(-std=c++11)
add_definitions(-Wall -Wconversion -Wno-unused-result)
add_definitions(-DBOOST_NO_CXX11_SCOPED_ENUMS)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
    add_definitions(-fsanitize=address)
    set(CMAKE_VERBOSE_MAKEFILE ON)
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# Native build.
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

## Definitions needed to eliminate dependency on the boost system library.
#add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
#add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Define our library.
add_library(embedDynamicLibrary SHARED $<TARGET_OBJECTS:embed_objlib>)
target_link_libraries(embedDynamicLibrary embed_objlib)
# set name
set_target_properties(embedDynamicLibrary PROPERTIES OUTPUT_NAME "embed")
set_target_properties(embedDynamicLibrary PROPERTIES PREFIX "")
set_target_properties(embedDynamicLibrary PROPERTIES DEFINE_SYMBOL "")
# Install the shared library into the bin directory.
install(TARGETS embedDynamicLibrary DESTINATION embed-install/bin)
