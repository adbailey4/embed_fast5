############################################################################################################
# testing
############################################################################################################
# my sources
set(TEST_CPP
        ${PROJECT_SOURCE_DIR}/tests/src/allTests.cpp
        ${PROJECT_SOURCE_DIR}/tests/src/BinaryIOTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/ConcurrentQueueTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/Fast5Tests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/EmbedUtilsTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/TopKmersTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/TestFiles.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/FileTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/MaxKmersTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/VariantPathTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/MarginalizeVariantsTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/FilterAlignmentsTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/PerPositionKmersTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/BaseKmerTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/BinaryEventTests.hpp
        ${PROJECT_SOURCE_DIR}/tests/src/AmbigModelTests.hpp)

add_executable(test_embed ${TEST_CPP})
if(BUILD_DYNAMIC_LIBRARY)
    target_link_libraries(test_embed PRIVATE embedDynamicLibrary)
else()
    target_link_libraries(test_embed PRIVATE embedStaticLibrary)
endif()

############################################################################################################
# cpp tests
add_gtest(test_embed ${CMAKE_SOURCE_DIR})
# python tests
message(STATUS PYTHON_EXECUTABLE: "${PYTHON_EXECUTABLE}")
add_test(NAME test_split_multi_read COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_split_multi_read.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME test_embedding_helpers COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_embedding_helpers.py
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
#install(TARGETS test_embed DESTINATION embed-install/bin)
