cmake_minimum_required(VERSION 3.5)
project(sohpianext)

set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(include)

add_library(sophianext SHARED
	src/sophia_data.cpp
	src/sophia_interface.cpp
)

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via Pybind11" ON)
find_package(PythonInterp)
find_package(PythonLibs)

if(ENABLE_PYTHON AND PYTHONLIBS_FOUND)
        include_directories(${PYTHON_INCLUDE_DIRS})
        add_subdirectory(lib/pybind11)

        pybind11_add_module(pysophia
                python/interface.cpp
        )

        set_target_properties(pysophia PROPERTIES
                                        PREFIX ""
                                        OUTPUT_NAME "pysophia"
                                        LINKER_LANGUAGE CXX)
        target_link_libraries(pysophia PRIVATE sophianext)

	execute_process(
	  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
		    from distutils import sysconfig as sc
		    print(sc.get_python_lib(prefix='', plat_specific=True))"
	  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
	  OUTPUT_STRIP_TRAILING_WHITESPACE)
        install(TARGETS pysophia LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}")
endif(ENABLE_PYTHON AND PYTHONLIBS_FOUND)

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

	# python tests
        if(ENABLE_PYTHON AND PYTHONLIBS_FOUND)
		CONFIGURE_FILE(test/testPythonInterface.py.in testPythonInterface.py)
		add_test(testPythonInterface ${PYTHON_EXECUTABLE} testPythonInterface.py)
	endif(ENABLE_PYTHON AND PYTHONLIBS_FOUND)
endif(ENABLE_TESTING)

# ----------------------------------------------------------------------------
# Executables
# ----------------------------------------------------------------------------

add_executable(run_sophia src/execute_sophia.cpp)
target_link_libraries(run_sophia sophianext)
