############################################################################################################
cmake_minimum_required(VERSION 3.8)
project(embed_fast5 VERSION 0.0.5)
############################################################################################################

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
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
include(cmake/AddDependencies.cmake)
############################################################################################################
# Define builds
############################################################################################################
option(BUILD_STATIC_LIBRARY "Build the static library." ON)
option(BUILD_DYNAMIC_LIBRARY "Build the dynamic library." ON)
option(BUILD_STATIC_EXECUTABLE "Build the shared executable." ON)
option(BUILD_DYNAMIC_EXECUTABLE "Build the dynamic executable." ON)

# Write out what we will actually build.
message(STATUS "BUILD_STATIC_LIBRARY is " ${BUILD_STATIC_LIBRARY})
message(STATUS "BUILD_STATIC_EXECUTABLE is " ${BUILD_STATIC_EXECUTABLE})
message(STATUS "BUILD_DYNAMIC_LIBRARY is " ${BUILD_DYNAMIC_LIBRARY})
message(STATUS "BUILD_DYNAMIC_EXECUTABLE is " ${BUILD_DYNAMIC_EXECUTABLE})

############################################################################################################
# Create object library
############################################################################################################
add_library(embed_objlib OBJECT
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

target_link_libraries(embed_objlib PUBLIC ${embed_LINK_LIBRARIES})
target_include_directories(embed_objlib
        PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/embed>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fast5/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fast5/include/fast5>
        )
# shared libraries need PIC
set_target_properties(embed_objlib PROPERTIES
        POSITION_INDEPENDENT_CODE 1)

target_compile_features(embed_objlib PRIVATE cxx_std_11)
target_compile_options(embed_objlib PRIVATE
        -Wall
        #        -Werror
        #        -Wno-error=maybe-uninitialized
        -Wconversion
        -Wno-unused-result
        -Wpedantic
        $<$<CONFIG:RELEASE>:-g0>
        $<$<CONFIG:RELEASE>:-O3>
        $<$<CONFIG:DEBUG>:-ggdb3>
        $<$<CONFIG:DEBUG>:-O0>
        $<$<CONFIG:DEBUG>:-ggdb3>
        )
target_compile_definitions(embed_objlib PRIVATE
        $<$<CONFIG:RELEASE>:NDEBUG>
        $<$<CONFIG:RELEASE>:BOOST_DISABLE_ASSERTS>
        )

############################################################################################################
# Create Static library
############################################################################################################
if(BUILD_STATIC_LIBRARY)
    add_library(embedStaticLibrary STATIC $<TARGET_OBJECTS:embed_objlib>)
    target_link_libraries(embedStaticLibrary embed_objlib)
    set_target_properties(embedStaticLibrary PROPERTIES
            OUTPUT_NAME "embedStaticLibrary"
            PREFIX ""
            DEFINE_SYMBOL ""
            VERSION ${PROJECT_VERSION})
endif(BUILD_STATIC_LIBRARY)

############################################################################################################
# Create dynamic library
############################################################################################################
if(BUILD_DYNAMIC_LIBRARY)
    add_library(embedDynamicLibrary SHARED $<TARGET_OBJECTS:embed_objlib>)
    target_link_libraries(embedDynamicLibrary embed_objlib)
    set_target_properties(embedDynamicLibrary PROPERTIES
            OUTPUT_NAME "embedDynamicLibrary"
            PREFIX ""
            DEFINE_SYMBOL ""
            VERSION ${PROJECT_VERSION})
endif(BUILD_DYNAMIC_LIBRARY)

############################################################################################################
# Create static executable
############################################################################################################
if (BUILD_STATIC_EXECUTABLE)
    add_executable(embed_main ${PROJECT_SOURCE_DIR}/src/main.cpp)
    target_link_libraries(embed_main PUBLIC embedStaticLibrary)
    target_include_directories(embed_main
            PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/fast5/include
            ${CMAKE_CURRENT_SOURCE_DIR}/fast5/include/fast5)
endif(BUILD_STATIC_EXECUTABLE)
############################################################################################################
# Create dynamic executable
############################################################################################################
if (BUILD_DYNAMIC_EXECUTABLE)
    add_executable(embed_main_dynamic ${PROJECT_SOURCE_DIR}/src/main.cpp)
    target_link_libraries(embed_main_dynamic PUBLIC embedDynamicLibrary)
endif(BUILD_DYNAMIC_EXECUTABLE)
############################################################################################################
# Python bindings
############################################################################################################
add_subdirectory(src/embed)
############################################################################################################
# Install
############################################################################################################
install(TARGETS embedStaticLibrary embedDynamicLibrary embed_main embed_main_dynamic embed_objlib
        EXPORT embedTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
############################################################################################################
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AmbigModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AlignmentFile.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AmbigModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AssignmentFile.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BaseKmer.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BinaryEventReader.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BinaryEventWriter.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BinaryIO.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ConcurrentQueue.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/EmbedFast5.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/EmbedUtils.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/EventDataHandler.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/FilterAlignments.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/FolderHandler.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/LoadVariantPaths.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/MarginalizeVariants.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/MaxKmers.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/PerPositionKmers.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/PositionsFile.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/PositionsKmerDistributions.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ReferenceHandler.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/SignalAlignToBed.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/SplitByRefPosition.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/TopKmers.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/VariantCall.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/VariantPath.hpp
        DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/embed
        )
############################################################################################################
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        cmake/embedConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

install(EXPORT embedTargets
        FILE embedTargets.cmake
        NAMESPACE embed::
        DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake/embed
        )

configure_file(cmake/embedConfig.cmake
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/embedConfig.cmake"
        COPYONLY
        )

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/embedConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/embedConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/embed
        )
############################################################################################################
# export to package registry
############################################################################################################
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE embed)

############################################################################################################
# uninstall
############################################################################################################
if(NOT TARGET uninstall)
    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake
            IMMEDIATE @ONLY)
    add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

############################################################################################################
# testing
############################################################################################################
enable_testing()
add_subdirectory(tests)

