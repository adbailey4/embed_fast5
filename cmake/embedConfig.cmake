############################################################################################################
message(CHECK_START "Finding Embed Dependencies")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
############################################################################################################
get_filename_component(embed_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
# TODO CONF_INCLUDE_DIRS is not working
set(embed_INCLUDE_DIRS "@CONF_INCLUDE_DIRS@")
message("${embed_INCLUDE_DIRS}")
include(CMakeFindDependencyMacro)
# Same syntax as find_package
#set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)
find_dependency(Boost 1.69.0 COMPONENTS system date_time filesystem context iostreams coroutine thread atomic REQUIRED)
find_dependency(nanopolish 0.11.1 REQUIRED nanopolish_static_lib nanopolishlib nanopolish_objlib)
############################################################################################################
include("${embed_CMAKE_DIR}/embedTargets.cmake")
############################################################################################################
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "all components found")


