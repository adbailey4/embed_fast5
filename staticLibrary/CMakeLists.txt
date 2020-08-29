cmake_minimum_required(VERSION 3.7)

project(embedStaticLibrary)

# C++ dialect.
set (CMAKE_CXX_STANDARD 11)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    #    add_definitions(-g)
    add_definitions(-O0)
#    add_definitions(-fsanitize=address)
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
    add_definitions(-DBOOST_NO_CXX11_SCOPED_ENUMS)
endif(BUILD_DEBUG)

# Native build.
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

# Build id.
add_definitions(-DBUILD_ID=${BUILD_ID})

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Define our library.
add_library(embedStaticLibrary STATIC $<TARGET_OBJECTS:embed_objlib>)
target_link_libraries(embedStaticLibrary embed_objlib)

# set name
set_target_properties(embedStaticLibrary PROPERTIES OUTPUT_NAME "embed")
set_target_properties(embedStaticLibrary PROPERTIES PREFIX "")
set_target_properties(embedStaticLibrary PROPERTIES DEFINE_SYMBOL "")
# Install the shared library into the bin directory.
install(TARGETS embedStaticLibrary DESTINATION embed-install/bin)
