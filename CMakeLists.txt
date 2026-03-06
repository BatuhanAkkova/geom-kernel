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

# Consolidated Tests
add_executable(test_primitives src/test_primitives.cpp)
target_link_libraries(test_primitives PRIVATE GeomKernel)

add_executable(test_booleans src/test_booleans.cpp)
target_link_libraries(test_booleans PRIVATE GeomKernel)

add_executable(test_smooth_booleans src/test_smooth_booleans.cpp)
target_link_libraries(test_smooth_booleans PRIVATE GeomKernel)

add_executable(test_tpms src/test_tpms.cpp)
target_link_libraries(test_tpms PRIVATE GeomKernel)

add_executable(test_core_utils src/test_core_utils.cpp)
target_link_libraries(test_core_utils PRIVATE GeomKernel)

add_executable(test_performance src/test_performance.cpp)
target_link_libraries(test_performance PRIVATE GeomKernel)

# Gallery Showcase
add_executable(gallery_showcase src/gallery_showcase.cpp)
target_link_libraries(gallery_showcase PRIVATE GeomKernel)
