cmake_minimum_required(VERSION 2.8)
project(minCut)

set(CMAKE_CXX_STANDARD 20)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_C_COMPILER /opt/scorep/bin/scorep-mpicc)
    set(CMAKE_CXX_COMPILER /opt/scorep/bin/scorep-mpicxx)
else()
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_CXX_COMPILER mpiCC)
endif()

set(GCC_COVERAGE_COMPILE_FLAGS " -O3 -funroll-loops -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(minCut main.cpp Graph.cpp Graph.hpp State.hpp State.cpp Solution.cpp Solution.hpp)
