############################################################################################################
cmake_minimum_required(VERSION 3.7)
project(embed_fast5)
############################################################################################################
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#message("CMAKE_MODULE_PATH Prefix: ${CMAKE_MODULE_PATH}")
set(CMAKE_INSTALL_PREFIX .)
message(STATUS "CMakeInstall Prefix: ${CMAKE_INSTALL_PREFIX}")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

############################################################################################################
# Find dependencies
############################################################################################################
include(AddDependencies)
include(AddGoogleTest)
add_subdirectory(submodules/pybind11)
find_package(nanopolish 0.11.1 COMPONENTS nanopolish_objlib nanopolish_static_lib nanopolishlib REQUIRED)
if (NOT nanopolish_FOUND)
    add_subdirectory(submodules/nanopolish)
else()
    message("nanopolish FOUND")
endif()

############################################################################################################
set(MAIN_CPP ${PROJECT_SOURCE_DIR}/src/main.cpp)
set(SOURCE
        ${PROJECT_SOURCE_DIR}/src/EmbedFast5.cpp ${PROJECT_SOURCE_DIR}/src/EmbedFast5.hpp
        ${PROJECT_SOURCE_DIR}/src/TopKmers.cpp ${PROJECT_SOURCE_DIR}/src/TopKmers.hpp
        ${PROJECT_SOURCE_DIR}/src/EmbedUtils.cpp ${PROJECT_SOURCE_DIR}/src/EmbedUtils.hpp
        ${PROJECT_SOURCE_DIR}/src/SignalAlignToBed.cpp ${PROJECT_SOURCE_DIR}/src/SignalAlignToBed.hpp
        ${PROJECT_SOURCE_DIR}/src/PositionsFile.cpp ${PROJECT_SOURCE_DIR}/src/PositionsFile.hpp
        ${PROJECT_SOURCE_DIR}/src/AlignmentFile.cpp ${PROJECT_SOURCE_DIR}/src/AlignmentFile.hpp
        ${PROJECT_SOURCE_DIR}/src/FilterAlignments.cpp ${PROJECT_SOURCE_DIR}/src/FilterAlignments.hpp
        ${PROJECT_SOURCE_DIR}/src/AssignmentFile.cpp ${PROJECT_SOURCE_DIR}/src/AssignmentFile.hpp
        ${PROJECT_SOURCE_DIR}/src/MaxKmers.cpp ${PROJECT_SOURCE_DIR}/src/MaxKmers.hpp
        ${PROJECT_SOURCE_DIR}/src/MarginalizeVariants.cpp ${PROJECT_SOURCE_DIR}/src/MarginalizeVariants.hpp
        ${PROJECT_SOURCE_DIR}/src/VariantCall.cpp ${PROJECT_SOURCE_DIR}/src/VariantCall.hpp
        ${PROJECT_SOURCE_DIR}/src/SplitByRefPosition.cpp ${PROJECT_SOURCE_DIR}/src/SplitByRefPosition.hpp
        ${PROJECT_SOURCE_DIR}/src/ReferenceHandler.cpp ${PROJECT_SOURCE_DIR}/src/ReferenceHandler.hpp
        ${PROJECT_SOURCE_DIR}/src/FolderHandler.cpp ${PROJECT_SOURCE_DIR}/src/FolderHandler.hpp
        ${PROJECT_SOURCE_DIR}/src/VariantPath.cpp ${PROJECT_SOURCE_DIR}/src/VariantPath.hpp
        ${PROJECT_SOURCE_DIR}/src/LoadVariantPaths.cpp ${PROJECT_SOURCE_DIR}/src/LoadVariantPaths.hpp
        ${PROJECT_SOURCE_DIR}/src/PerPositionKmers.hpp
        ${PROJECT_SOURCE_DIR}/src/ConcurrentQueue.hpp
        ${PROJECT_SOURCE_DIR}/src/BinaryIO.hpp
        ${PROJECT_SOURCE_DIR}/src/BinaryEventWriter.hpp
        ${PROJECT_SOURCE_DIR}/src/BinaryEventReader.hpp
        ${PROJECT_SOURCE_DIR}/src/BaseKmer.hpp
        ${PROJECT_SOURCE_DIR}/src/EventDataHandler.hpp
        ${PROJECT_SOURCE_DIR}/src/PositionsKmerDistributions.cpp
        ${PROJECT_SOURCE_DIR}/src/PositionsKmerDistributions.hpp
        ${PROJECT_SOURCE_DIR}/src/AmbigModel.hpp)

############################################################################################################
option(BUILD_STATIC_LIBRARY "Build the static library." ON)
option(BUILD_DYNAMIC_LIBRARY "Build the dynamic library." ON)
option(BUILD_STATIC_EXECUTABLE "Build the shared executable." ON)
option(BUILD_DYNAMIC_EXECUTABLE "Build the dynamic executable." OFF)

# Write out what we will actually build.
message(STATUS "BUILD_STATIC_LIBRARY is " ${BUILD_STATIC_LIBRARY})
message(STATUS "BUILD_STATIC_EXECUTABLE is " ${BUILD_STATIC_EXECUTABLE})
message(STATUS "BUILD_DYNAMIC_LIBRARY is " ${BUILD_DYNAMIC_LIBRARY})
message(STATUS "BUILD_DYNAMIC_EXECUTABLE is " ${BUILD_DYNAMIC_EXECUTABLE})

if(BUILD_STATIC_EXECUTABLE AND BUILD_DYNAMIC_EXECUTABLE)
    message(FATAL_ERROR "BUILD_STATIC_EXECUTABLE and BUILD_DYNAMIC_EXECUTABLE cannot both be set to ON. Exiting" )
endif()

if(NOT BUILD_STATIC_EXECUTABLE AND NOT BUILD_DYNAMIC_EXECUTABLE)
    message(FATAL_ERROR "BUILD_STATIC_EXECUTABLE and BUILD_DYNAMIC_EXECUTABLE cannot both be set to OFF. Exiting" )
endif()
# this is the "object library" target: compiles the sources only once
add_library(embed_objlib OBJECT ${SOURCE})
# my includes
target_link_libraries(embed_objlib PUBLIC
        Boost::system
        Boost::date_time
        Boost::filesystem
        Boost::context
        Boost::iostreams
        Boost::coroutine)

if (nanopolish_FOUND)
    target_link_libraries(embed_objlib PUBLIC nanopolish::nanopolish_static_lib)
else()
    target_link_libraries(embed_objlib PUBLIC nanopolishlib)
endif()

            target_include_directories(embed_objlib PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/fast5/include
        ${CMAKE_SOURCE_DIR}/fast5/include/fast5)
# shared libraries need PIC
set_property(TARGET embed_objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

# Add the subdirectories we need.
if(BUILD_STATIC_LIBRARY)
    add_subdirectory(staticLibrary)
endif(BUILD_STATIC_LIBRARY)
if(BUILD_DYNAMIC_LIBRARY)
    add_subdirectory(dynamicLibrary)
endif(BUILD_DYNAMIC_LIBRARY)
############################################################################################################
add_subdirectory(src/embed)
enable_testing()
add_subdirectory(tests)

add_executable(embedExecutable ${MAIN_CPP})
set_target_properties(embedExecutable
        PROPERTIES
        OUTPUT_NAME "embed_main"
        SUFFIX "")

if(BUILD_DYNAMIC_EXECUTABLE)
    target_link_libraries(embedExecutable PRIVATE embedDynamicLibrary)
else()
    target_link_libraries(embedExecutable PRIVATE embedStaticLibrary)
endif()

install(TARGETS embedExecutable DESTINATION embed-install/bin)

