############################################################################################################
message(CHECK_START "Finding Embed Dependencies")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
############################################################################################################
# GOOGLE TEST LIBRARY
############################################################################################################
include(${PROJECT_SOURCE_DIR}/cmake/AddGoogleTest.cmake)
############################################################################################################
# BOOST LIBRARY
############################################################################################################
#set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.69.0 COMPONENTS system date_time filesystem context iostreams coroutine REQUIRED)
############################################################################################################
# NANOPOLISH LIBRARY
############################################################################################################
find_package(nanopolish 0.11.1 COMPONENTS nanopolish_objlib nanopolish_static_lib nanopolishlib REQUIRED)
if (NOT nanopolish_FOUND)
    add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/nanopolish)
else()
    message(STATUS "nanopolish FOUND: ${nanopolish_CONFIG}")
endif()
############################################################################################################
# Pybind LIBRARY
############################################################################################################
add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/pybind11)
############################################################################################################
set(embed_LINK_LIBRARIES
        Boost::system
        Boost::date_time
        Boost::filesystem
        Boost::context
        Boost::iostreams
        Boost::coroutine)
if (nanopolish_FOUND)
    list(APPEND embed_LINK_LIBRARIES nanopolish::nanopolish_static_lib)
else()
    list(APPEND embed_LINK_LIBRARIES nanopolish)
endif()

############################################################################################################
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "all components found")
