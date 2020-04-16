cmake_minimum_required(VERSION 3.5)
project(sohpianext)

set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(include)

add_library(sophianext SHARED
	src/sophia_data.cpp
	src/sophia_interface.cpp
)

add_executable(run_sophia src/execute_sophia.cpp)
target_link_libraries(run_sophia sophianext)
