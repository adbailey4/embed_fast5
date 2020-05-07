############################################################################################################
cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_FLAGS "-O3 -g -Wall -std=c++11 -fopenmp -fsigned-char -DBOOST_NO_CXX11_SCOPED_ENUMS")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_FLAGS_DEBUG  "-O0 -g -Wall --pedantic -fopenmp -funroll-loops -DDEBUG")
project(embed_fast5)
############################################################################################################
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -Wall -std=c++11 -fopenmp -fsigned-char -DBOOST_NO_CXX11_SCOPED_ENUMS")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-3")
message("CMakeInstall Prefix: ${CMAKE_INSTALL_PREFIX}")
############################################################################################################
# GTEST INSTALL
############################################################################################################
if(APPLE)
    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
endif(APPLE)
# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()
# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)
# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()
message(gtest_SOURCE_DIR: ${gtest_SOURCE_DIR})
############################################################################################################
# BOOST LIBRARY
############################################################################################################
#set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.69.0 COMPONENTS system date_time filesystem context iostreams coroutine REQUIRED)
message(Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS})
message(Boost_LIBRARIES: ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})
############################################################################################################
# ZLIB LIBRARY
############################################################################################################
find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()
message(Zlib: ${ZLIB_LIBRARIES})
############################################################################################################
# DL LIBRARY
############################################################################################################
find_library(LIBDL_LIBRARY NAMES dl)
message(LibDL: ${LIBDL_LIBRARY})
############################################################################################################
# M LIBRARY
############################################################################################################
find_library(LIBM_LIBRARY NAMES m)
message(LibM: ${LIBM_LIBRARY})
############################################################################################################
# pthread LIBRARY
############################################################################################################
find_library(LIBPTHREAD_LIBRARY NAMES pthread)
message(LibPthread: ${LIBPTHREAD_LIBRARY})
############################################################################################################
# BZ2 LIBRARY
############################################################################################################
find_library(LIBBZ2_LIBRARY NAMES bz2)
message(LibBZ2: ${LIBBZ2_LIBRARY})
############################################################################################################
# CURL LIBRARY
############################################################################################################
find_library(LIBCURL_LIBRARY NAMES curl)
message(LibCURL: ${LIBCURL_LIBRARY})
############################################################################################################
# LZMA LIBRARY
############################################################################################################
find_library(LIBLZMA_LIBRARY NAMES lzma)
message(LibLZMA: ${LIBLZMA_LIBRARY})
############################################################################################################
# HDF LIBRARY
############################################################################################################
find_library(HDF5_LIBRARIES libhdf5.a
        HINTS ${NANOPOLISH_HOME}/lib
        ${CMAKE_INSTALL_PREFIX}/lib
        /usr/lib
        ${LIBDIR})

find_path(HDF5_INCLUDE_DIRS hdf5.h
        HINTS ${HDF5_LIBRARIES}
        "${NANOPOLISH_HOME}/include"
        ${CMAKE_INSTALL_PREFIX}/include
        ${INCLUDEDIR}
        /usr/include)
message(HDF5_INCLUDE_DIRS: "${HDF5_INCLUDE_DIRS}")
message(HDF5_LIBRARIES: "${HDF5_LIBRARIES}")
include_directories(${HDF5_INCLUDE_DIRS})
############################################################################################################
# EIGEN LIBRARY
############################################################################################################
find_path(EIGEN3_INCLUDE_DIR Eigen HINTS
        ${CMAKE_INSTALL_PREFIX}/include
        ${NANOPOLISH_HOME}
        ${INCLUDEDIR}
        PATH_SUFFIXES eigen3 eigen)

include_directories(${EIGEN3_INCLUDE_DIR})
message(EIGEN3_INCLUDE_DIR: "${EIGEN3_INCLUDE_DIR}")
############################################################################################################
# HTSLIB LIBRARY
############################################################################################################
find_library(HTSLIB_LIBRARY NAMES libhts.a
        HINTS "${NANOPOLISH_HOME}/htslib"
        ${CMAKE_INSTALL_PREFIX}/lib
        /usr/lib
        ${LIBDIR})

find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
        HINTS ${HTSLIB_LIBRARY}
        ${CMAKE_INSTALL_PREFIX}/include
        "${NANOPOLISH_HOME}/htslib"
        /usr/include
        ${INCLUDEDIR})

message(HTSLIB_LIBRARIES: "${HTSLIB_LIBRARY}")
message(HTSLIB_INCLUDE_DIRS: "${HTSLIB_INCLUDE_DIR}")
include_directories(${HTSLIB_INCLUDE_DIRS})
############################################################################################################
## OPENSSL_LIB LIBRARY
#############################################################################################################
#find_library(OPENSSL_LIB libcrypto.a
#        HINTS "${NANOPOLISH_HOME}/lib"
#        ${CMAKE_INSTALL_PREFIX}/lib
#        /usr/local/opt/openssl/lib
#        /usr/lib
#        ${LIBDIR})
#message(OPENSSL_LIB: ${OPENSSL_LIB})
#############################################################################################################
### Curl LIBRARY
##############################################################################################################
#
#find_library(LIBCURL_LIBRARY NAMES curl)
#message(LIBCURL_LIB: ${LIBCURL_LIBRARY})

############################################################################################################
# NANOPOLISH LIBRARY
############################################################################################################
set(NANOPOLISH_HOME ${CMAKE_SOURCE_DIR}/submodules/nanopolish)
# nanopolish includes
include_directories(${NANOPOLISH_HOME}/src)
include_directories(${NANOPOLISH_HOME}/src/alignment)
include_directories(${NANOPOLISH_HOME}/src/common)
include_directories(${NANOPOLISH_HOME}/src/hmm)
include_directories(${NANOPOLISH_HOME}/src/main)
include_directories(${NANOPOLISH_HOME}/src/pore_model)
include_directories(${NANOPOLISH_HOME}/src/test)
include_directories(${NANOPOLISH_HOME}/src/thirdparty)
include_directories(${NANOPOLISH_HOME}/src/thirdparty/scrappie)

set(nanopolish_source
        ${NANOPOLISH_HOME}/src/nanopolish_call_variants.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_call_methylation.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_extract.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_getmodel.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_haplotype.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_index.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_methyltrain.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_phase_reads.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_raw_loader.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_read_db.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_variant_db.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_train_poremodel_from_basecalls.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_squiggle_read.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_scorereads.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_polya_estimator.cpp
        ${NANOPOLISH_HOME}/src/nanopolish_vcf2fasta.cpp
        ${NANOPOLISH_HOME}/src/hmm/nanopolish_duration_model.cpp
        ${NANOPOLISH_HOME}/src/hmm/nanopolish_profile_hmm.cpp
        ${NANOPOLISH_HOME}/src/hmm/nanopolish_profile_hmm_r7.cpp
        ${NANOPOLISH_HOME}/src/hmm/nanopolish_profile_hmm_r9.cpp
        ${NANOPOLISH_HOME}/src/hmm/nanopolish_transition_parameters.cpp
        ${NANOPOLISH_HOME}/src/common/fs_support.cpp
        ${NANOPOLISH_HOME}/src/common/logsum.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_alphabet.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_bam_processor.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_bam_utils.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_common.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_fast5_io.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_iupac.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_klcs.cpp
        ${NANOPOLISH_HOME}/src/common/nanopolish_variant.cpp
        ${NANOPOLISH_HOME}/src/alignment/nanopolish_alignment_db.cpp
        ${NANOPOLISH_HOME}/src/alignment/nanopolish_anchor.cpp
        ${NANOPOLISH_HOME}/src/alignment/nanopolish_eventalign.cpp
        ${NANOPOLISH_HOME}/src/pore_model/nanopolish_poremodel.cpp
        ${NANOPOLISH_HOME}/src/pore_model/nanopolish_pore_model_set.cpp
        ${NANOPOLISH_HOME}/src/pore_model/nanopolish_model_names.cpp
        ${NANOPOLISH_HOME}/src/training_core.cpp
        ${NANOPOLISH_HOME}/src/thirdparty/scrappie/scrappie_common.c
        ${NANOPOLISH_HOME}/src/thirdparty/scrappie/event_detection.c
        ${NANOPOLISH_HOME}/src/thirdparty/stdaln.c)

set(NANOPOLISH ${NANOPOLISH_HOME}/src/main/nanopolish.cpp)
set(NP_LIBRARY ${HDF5_LIBRARIES} ${HTSLIB_LIBRARY} ${LIBCURL_LIBRARY} ${ZLIB_LIBRARIES} ${LIBLZMA_LIBRARY} ${LIBBZ2_LIBRARY} ${ZLIB_LIBRARIES} ${LIBM_LIBRARY} ${LIBDL_LIBRARY} ${LIBPTHREAD_LIBRARY})
# nanopolish main
#add_executable(nanopolish ${NANOPOLISH} ${nanopolish_source})
#target_link_libraries(nanopolish ${NP_LIBRARY})
############################################################################################################
# my library
set(EMBED_LIBRARY ${Boost_LIBRARIES} ${NP_LIBRARY})
# my includes
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/fast5/include)
include_directories(${CMAKE_SOURCE_DIR}/fast5/include/fast5)
# my sources
set(TEST_CPP tests/src/allTests.cpp
        tests/src/BinaryFileTests.hpp
        tests/src/ConcurrentQueueTests.hpp
        tests/src/Fast5Tests.hpp
        tests/src/EmbedUtilsTests.hpp
        tests/src/TopKmersTests.hpp
        tests/src/TestFiles.hpp
        tests/src/FileTests.hpp
        tests/src/MaxKmersTests.hpp
        tests/src/VariantPathTests.hpp
        tests/src/MarginalizeVariantsTests.hpp
        tests/src/FilterAlignmentsTests.hpp
        tests/src/PerPositionKmersTests.hpp)

set(MAIN_CPP src/main.cpp)
set(SOURCE
        src/EmbedFast5.cpp src/EmbedFast5.hpp
        src/TopKmers.cpp src/TopKmers.hpp
        src/EmbedUtils.cpp src/EmbedUtils.hpp
        src/SignalAlignToBed.cpp src/SignalAlignToBed.hpp
        src/PositionsFile.cpp src/PositionsFile.hpp
        src/AlignmentFile.cpp src/AlignmentFile.hpp
        src/FilterAlignments.cpp src/FilterAlignments.hpp
        src/AssignmentFile.cpp src/AssignmentFile.hpp
        src/MaxKmers.cpp src/MaxKmers.hpp
        src/MarginalizeVariants.cpp src/MarginalizeVariants.hpp
        src/VariantCall.cpp src/VariantCall.hpp
        src/SplitByRefPosition.cpp src/SplitByRefPosition.hpp
        src/PerPositonKmers.cpp src/PerPositonKmers.hpp
        src/ReferenceHandler.cpp src/ReferenceHandler.hpp
        src/FolderHandler.cpp src/FolderHandler.hpp
        src/VariantPath.cpp src/VariantPath.hpp
        src/LoadVariantPaths.cpp src/LoadVariantPaths.hpp
        src/ConcurrentQueue.hpp
        src/BinaryIO.hpp
        src/BinaryEventWriter.cpp src/BinaryEventWriter.hpp
        src/BinaryEventReader.cpp src/BinaryEventReader.hpp)
############################################################################################################
# main
add_executable(embed_main ${MAIN_CPP} ${nanopolish_source} ${SOURCE})
target_link_libraries(embed_main ${EMBED_LIBRARY})
# pybind
add_subdirectory(submodules/pybind11)
pybind11_add_module(bindings src/embed/PythonBindings.cpp ${nanopolish_source} ${SOURCE})
target_link_libraries(bindings PRIVATE ${EMBED_LIBRARY} ${HDF5_LIBRARIES})

############################################################################################################

# testing
add_executable(test_embed ${TEST_CPP} ${SOURCE} ${nanopolish_source})
target_link_libraries(test_embed ${HDF5_LIBRARIES} gtest_main gmock_main ${EMBED_LIBRARY})
enable_testing()
# cpp tests
add_test(cpp_tests test_embed ${CMAKE_SOURCE_DIR})
# python tests
message(PYTHON_EXECUTABLE: "${PYTHON_EXECUTABLE}")
add_test(NAME test_split_multi_read COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_split_multi_read.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME test_embedding_helpers COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_embedding_helpers.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# make check command for testing
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
