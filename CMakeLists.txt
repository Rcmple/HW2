cmake_minimum_required(VERSION 3.28)
project(cpp_project_2)

set(CMAKE_CXX_STANDARD 20)

add_executable(cpp_project_2 my_fluid.cpp)

target_compile_definitions(cpp_project_2 PRIVATE "TYPES=FLOAT,DOUBLE,FAST_FIXED(32,5)" "SIZES=S(5,4),S(7,4)")

target_include_directories(cpp_project_2 PRIVATE ${CMAKE_SOURCE_DIR}/mnt/d/Rcmple/cpp_project_2/FixedClass.h)