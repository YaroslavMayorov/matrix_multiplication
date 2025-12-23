cmake_minimum_required(VERSION 3.10)
project(MatrixMultiplication)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(include)

add_executable(matrix_algo src/main.cpp src/Algorithms.cpp src/Matrix.cpp)

find_library(ACCELERATE_FRAMEWORK Accelerate)
if(ACCELERATE_FRAMEWORK)
    target_link_libraries(matrix_algo ${ACCELERATE_FRAMEWORK})
    target_compile_definitions(matrix_algo PRIVATE HAS_ACCELERATE)
endif()
