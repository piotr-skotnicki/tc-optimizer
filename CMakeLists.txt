# TC Optimizing Compiler

project(tc C CXX)

cmake_minimum_required(VERSION 2.8)

include_directories(
    isl/include
    pet/include
    barvinok)

add_executable(tc
    src/main.cpp
    src/utility.cpp
    src/scop.cpp
    src/options.cpp
    src/debug.cpp
    src/tiling.cpp
    src/codegen.cpp
    src/ast.cpp
    src/slicing.cpp
    src/regular_tiling.cpp
    src/stencil_tiling.cpp
    src/merge_tiling.cpp
    src/correction_tiling.cpp
    src/correction_inv_tiling.cpp
    src/scheduling.cpp
    src/lex_scheduling.cpp
    src/sfs_scheduling.cpp
    src/free_scheduling.cpp
    src/dynamic_free_scheduling.cpp
    src/serial_codegen.cpp
    src/omp_cpu_codegen.cpp
    src/transitive_closure.cpp
    src/tarjan_transitive_closure.cpp
    src/floyd_warshall_transitive_closure.cpp
    src/iterative_transitive_closure.cpp
    src/tile_statistics.cpp
    src/for_decorator.cpp
    src/timer.cpp)

include(ExternalProject)

ExternalProject_Add(isl
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/isl
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/isl
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/isl/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/isl
    BUILD_COMMAND ${MAKE})

ExternalProject_Add(pet
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pet
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/pet
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/pet/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/pet
    BUILD_COMMAND ${MAKE})

ExternalProject_Add(barvinok
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/barvinok
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/barvinok
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/barvinok/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/barvinok
    BUILD_COMMAND ${MAKE})

add_dependencies(tc
    isl
    pet
    barvinok)

find_library(glpk NAMES glpk libglpk)
find_library(ntl NAMES ntl libntl)
find_library(gmp NAMES gmp libgmp)

if(NOT glpk)
    message(FATAL_ERROR "glpk library not found")
endif()

if(NOT ntl)
    message(FATAL_ERROR "ntl library not found")
endif()

if(NOT gmp)
    message(FATAL_ERROR "gmp library not found")
endif()

target_link_libraries(tc
    ${CMAKE_CURRENT_BINARY_DIR}/barvinok/lib/libbarvinok.a
    ${CMAKE_CURRENT_BINARY_DIR}/isl/lib/libisl.a
    ${CMAKE_CURRENT_BINARY_DIR}/pet/lib/libpet.so
    ${CMAKE_CURRENT_BINARY_DIR}/barvinok/lib/libpolylibgmp.a
    ${glpk}
    ${ntl}
    ${gmp})
