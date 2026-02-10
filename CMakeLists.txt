cmake_minimum_required(VERSION 3.20)
project(GeomKernel VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library
add_library(GeomKernel 
    src/MCTables.cpp
)
target_include_directories(GeomKernel PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Test / Example Executable
add_executable(kernel_test src/main.cpp)
target_link_libraries(kernel_test PRIVATE GeomKernel)

add_executable(mvp_test src/mvp_test.cpp)
target_link_libraries(mvp_test PRIVATE GeomKernel)

# Phase 2 Test
add_executable(phase2_test src/phase2_test.cpp)
target_link_libraries(phase2_test PRIVATE GeomKernel)

# Phase 3 Benchmark
add_executable(phase3_benchmark src/phase3_benchmark.cpp)
target_link_libraries(phase3_benchmark PRIVATE GeomKernel)
