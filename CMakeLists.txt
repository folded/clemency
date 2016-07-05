cmake_minimum_required(VERSION 2.6)

add_subdirectory(gloop)
include_directories(${gloop_SOURCE_DIR}/include)

project(clemency)

set(clemency_VERSION_MAJOR 1)
set(clemency_VERSION_MINOR 0)
set(clemency_VERSION_PATCH 0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${clemency_SOURCE_DIR}/cmake/Modules/")

set(CLEMENCY_VERSION ${clemency_VERSION_MAJOR}.${clemency_VERSION_MINOR}.${clemency_VERSION_PATCH})

option(BUILD_COVERAGE                    "Compile with gcov"                                 OFF)

enable_language(CXX)

if (BUILD_COVERAGE)
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_BUILD_TYPE DEBUG)
  IF (CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
  ENDIF()
endif(BUILD_COVERAGE)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.40 REQUIRED COMPONENTS program_options)

if (Boost_FOUND)
  if (!Boost_PROGRAM_OPTIONS_FOUND)
    message(ERROR "Boost program_options not found")
  endif ()

  include_directories(${Boost_INCLUDE_DIRS})
else (Boost_FOUND)
  message(ERROR "Boost not found")
endif (Boost_FOUND)

find_package(OpenCL REQUIRED)
find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

configure_file (
  "${clemency_SOURCE_DIR}/include/clemency/cmake-config.h.in"
  "${clemency_BINARY_DIR}/include/clemency/config.h"
)
include_directories(${clemency_BINARY_DIR}/include)

add_definitions(-DCMAKE_BUILD)

add_subdirectory(include) 
add_subdirectory(lib)
add_subdirectory(src) 
add_subdirectory(tests)
