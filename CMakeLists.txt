cmake_minimum_required(VERSION 3.23)
project(merrill2)

set(CMAKE_CXX_STANDARD 20)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(MERRILL2_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(MERRILL2_SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(MERRILL2_THIRD_PARTY_DIR "${CMAKE_SOURCE_DIR}/third-party")
set(MERRILL2_TEST_SRC_DIR "${CMAKE_SOURCE_DIR}/test-src")

include_directories(${MERRILL2_INCLUDE_DIR})

if (${MULTIPRECISION})

    message(STATUS "Building a multiprecision version of merrill2.")

    # Try to load the multiprecision library.
    find_package(MPFR REQUIRED)

    # Print some messages.
    message(STATUS "Found MPFR, details follow ...")
    message(STATUS "MPFR_INCLUDES: ${MPFR_INCLUDES}")
    message(STATUS "MPFR_LIBRARIES: ${MPFR_LIBRARIES}")
    message(STATUS "MPFR_VERSION: ${MPFR_VERSION}")

    # Add mpreal.
    set(MPREAL_INCLUDE_DIR "${MERRILL2_THIRD_PARTY_DIR}/mpreal-v4.2.0/include")
    message(STATUS "MPREAL_INCLUDE_DIR: ${MPREAL_INCLUDE_DIR}")

    include_directories(${MPFR_INCLUDES})
    include_directories(${MPREAL_INCLUDE_DIR})

    add_compile_definitions(WITH_MULTIPRECISION=1)

endif()

if (${BUILD_TESTS})

    message(STATUS "Building tests.")

    set(MERRILL2_TEST_SRC_DIR "${CMAKE_SOURCE_DIR}/test-src")
    set(MERRILL2_TEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR/test-include}")
    set(CATCH_INCLUDE_DIR "${MERRILL2_THIRD_PARTY_DIR}/catch2-v2.13.9/include")

    include_directories(${CATCH_INCLUDE_DIR})

    add_subdirectory(${MERRILL2_TEST_SRC_DIR})

endif()

if (${DEBUG_MESSAGES})

    message(STATUS "Building with debug messages enabled.")
    add_compile_definitions(DEBUG_MESSAGES=1)

endif()

add_executable(merrill2 main.cpp)
target_link_libraries(merrill2
        ${MPFR_LIBRARIES})
