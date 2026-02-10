cmake_minimum_required(VERSION 3.20)
project(GeomKernel VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library
add_library(GeomKernel INTERFACE)
target_include_directories(GeomKernel INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Test / Example Executable
add_executable(kernel_test src/main.cpp)
target_link_libraries(kernel_test PRIVATE GeomKernel)
