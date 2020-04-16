cmake_minimum_required(VERSION 3.5)
project(sohpianext)

set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(include)

add_library(sophianext SHARED
	src/sophia_data.cpp
	src/sophia_interface.cpp
)

# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
option(ENABLE_TESTING "Build tests and enable test target" OFF)
if(ENABLE_TESTING)
	enable_testing()

        include_directories(lib/gtest/googletest/include)
        add_subdirectory(lib/gtest/googletest)
        if(APPLE)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
        endif(APPLE)
        
        add_executable(testFirst test/testFirst.cpp)
        target_link_libraries(testFirst sophianext gtest gtest_main)
        add_test(testFirst testFirst)
endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Executables
# ----------------------------------------------------------------------------

add_executable(run_sophia src/execute_sophia.cpp)
target_link_libraries(run_sophia sophianext)
