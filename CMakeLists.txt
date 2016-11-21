cmake_minimum_required(VERSION 2.6)
project(tests)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(tests ${PROJECT_SOURCE_DIR}/src/tests/main.cpp)
