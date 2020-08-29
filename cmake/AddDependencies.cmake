message(CHECK_START "Finding Boost Package")
############################################################################################################
# BOOST LIBRARY
############################################################################################################
#set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.69.0 COMPONENTS system date_time filesystem context iostreams coroutine REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "boost package not found")
else()
    message(CHECK_PASS "found: ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")
endif()
