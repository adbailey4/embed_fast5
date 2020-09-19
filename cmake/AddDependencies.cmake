############################################################################################################
message(CHECK_START "Finding Embed Dependencies")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

if(BUILD_SHARED_LIBS)
    message(STATUS "BUILDING SHARED LIBS")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(Boost_USE_STATIC_LIBS        OFF)
    set(Boost_USE_MULTITHREADED      OFF)
    set(Boost_USE_STATIC_RUNTIME    OFF)

else(NOT BUILD_SHARED_LIBS)
    message(STATUS "BUILDING STATIC LIBS")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(Boost_USE_STATIC_LIBS        ON)
    set(Boost_USE_MULTITHREADED      OFF)
    set(Boost_USE_STATIC_RUNTIME    ON)

endif(BUILD_SHARED_LIBS)

############################################################################################################
# GOOGLE TEST LIBRARY
############################################################################################################
set(gtest_hide_internal_symbols ON)
include(${PROJECT_SOURCE_DIR}/cmake/AddGoogleTest.cmake)
############################################################################################################
# BOOST LIBRARY
############################################################################################################
#find_package(Boost 1.69.0 COMPONENTS system date_time filesystem context iostreams coroutine thread atomic REQUIRED)
find_package(Boost 1.65.1 COMPONENTS system date_time filesystem context iostreams coroutine thread atomic REQUIRED)

############################################################################################################
# NANOPOLISH LIBRARY
############################################################################################################
#find_package(nanopolish 0.11.1 COMPONENTS nanopolish_objlib nanopolish_static_lib nanopolishlib)
if (NOT nanopolish_FOUND)
    message(STATUS "adding nanopolish subdirectory:")
    add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/nanopolish)
else()
    message(STATUS "nanopolish FOUND: ${nanopolish_CONFIG}")
endif()

############################################################################################################
# Pybind LIBRARY
############################################################################################################
find_package(PythonLibs REQUIRED)
add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/pybind11)
############################################################################################################
set(embed_LINK_LIBRARIES
        Boost::system
        Boost::date_time
        Boost::filesystem
        Boost::context
        Boost::iostreams
        Boost::coroutine
        Boost::thread
        Boost::atomic)

if (nanopolish_FOUND)
    set(nanopolish_LIB nanopolish::nanopolishlib)
else()
    set(nanopolish_LIB nanopolishlib)
endif()


############################################################################################################
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "all components found")
