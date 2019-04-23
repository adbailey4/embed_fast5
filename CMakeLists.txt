cmake_minimum_required(VERSION 3.7)
project(embed_fast5)

set(CMAKE_C_STANDARD 99)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -fopenmp -fsigned-char -DBOOST_NO_CXX11_SCOPED_ENUMS")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-3")


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


# Boost libraries
#add_subdirectory(submodules/boost-cmake)
set(Boost_USE_MULTITHREADED      OFF)

set(NANOPOLISH_HOME ${CMAKE_SOURCE_DIR}/submodules/nanopolish)

find_package(Boost COMPONENTS system date_time filesystem iostreams REQUIRED)
#find_package(hdf5 1.8 REQUIRED HINTS
#        ${CMAKE_INSTALL_PREFIX}/include
#        )

find_library(HDF5_LIBRARIES libhdf5.a
        HINTS "${NANOPOLISH_HOME}/lib"
        ${CMAKE_INSTALL_PREFIX}/lib
        /usr/lib)

find_path(HDF5_INCLUDE_DIRS hdf5.h
        HINTS ${HDF5_LIBRARIES}
        "${NANOPOLISH_HOME}/include"
        ${CMAKE_INSTALL_PREFIX}/include
        /usr/include)


message("CMakeInstall Prefix: ${CMAKE_INSTALL_PREFIX}")
message(Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS})
message(HDF5_INCLUDE_DIRS: "${HDF5_INCLUDE_DIRS}")
message(HDF5_LIBRARIES: "${HDF5_LIBRARIES}")

find_path(EIGEN3_INCLUDE_DIR Eigen HINTS
        ${CMAKE_INSTALL_PREFIX}/include
        ${NANOPOLISH_HOME}
        PATH_SUFFIXES eigen3 eigen)

message(EIGEN3_INCLUDE_DIR: "${EIGEN3_INCLUDE_DIR}")

find_library(HTSLIB_LIBRARY NAMES libhts.a
        HINTS "${NANOPOLISH_HOME}/htslib"
        ${CMAKE_INSTALL_PREFIX}/include)

find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
        HINTS ${HTSLIB_LIBRARY}
        ${CMAKE_INSTALL_PREFIX}/include
        "${NANOPOLISH_HOME}/htslib"
        /usr/include)


set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY} )
set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR} )

message(HTSLIB_LIBRARIES: "${HTSLIB_LIBRARIES}")
message(HTSLIB_INCLUDE_DIRS: "${HTSLIB_INCLUDE_DIRS}")

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${HTSLIB_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${Eigen3_INCLUDE_DIRS})

link_directories(/usr/local/lib)
link_directories(/usr/lib)

# my includes
include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/fast5/include)
include_directories(${CMAKE_SOURCE_DIR}/fast5/include/fast5)

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
include_directories(/usr/local/include)
include_directories(/usr/include)



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

set(SOURCE_CPP
        impl/embed_fast5.cpp)

set(CPP_HEADERS
        inc/embed_fast5.hpp)

set(TEST_CPP
        tests/allTests.cpp)

set(MAIN_CPP main.cpp)

set(NANOPOLISH
        ${NANOPOLISH_HOME}/src/main/nanopolish.cpp)

set(NP_LIBRARY -lz -lm -ldl -lpthread -lbz2 -lcurl -llzma ${HDF5_LIBRARIES} ${HTSLIB_LIBRARY} )
set(EMBED_LIBRARY ${NP_LIBRARY} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARIES} )


# nanopolish main
add_executable(nanopolish ${NANOPOLISH} ${nanopolish_source})
target_link_libraries(nanopolish ${NP_LIBRARY})

# main and tests
add_executable(main_cpp ${MAIN_CPP} ${nanopolish_source} ${SOURCE_CPP} ${CPP_HEADERS})
target_link_libraries(main_cpp ${EMBED_LIBRARY})

add_executable(test_cpp ${TEST_CPP} ${SOURCE_CPP} ${CPP_HEADERS} ${nanopolish_source})
target_link_libraries(test_cpp ${HDF5_LIBRARIES} gtest_main ${EMBED_LIBRARY})

# testing

enable_testing()
add_test(cpp_tests test_cpp ${CMAKE_SOURCE_DIR})
add_test(NAME python_tests COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/embed/test/test_split_multi_read.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
