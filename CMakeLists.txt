cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(run_sophia ${SOURCES})
