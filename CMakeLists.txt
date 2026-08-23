cmake_minimum_required(VERSION 3.5)
project(SSHASH-LITE)

set(CMAKE_CXX_STANDARD 20)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (SSHASH_QUERY_VERBOSE_OUTPUT)
  MESSAGE(STATUS "SSHash query verbose output")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D SSHASH_QUERY_VERBOSE_OUTPUT")
endif()

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces")

  # Flags for PTHash:
  if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")        # for multi-threading

  if (SSHASH_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

# create a static lib
add_library(sshashlite_static STATIC
  include/gz/zip_stream.cpp
  include/dictionary.cpp
  include/info.cpp
  include/builder/build.cpp
)

# build with the lib
add_executable(sshash-lite src/sshash-lite.cpp)
target_link_libraries(sshash-lite
  sshashlite_static
  z
)

include(CTest)
if (BUILD_TESTING)
  add_executable(bigreedy-ms-regression tests/bigreedy_ms_regression.cpp)
  target_link_libraries(bigreedy-ms-regression
    sshashlite_static
    z
  )
  add_test(NAME bigreedy-ms-regression COMMAND bigreedy-ms-regression)
endif()
