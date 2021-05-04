cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME MeshUtilityCPP)
project (${PROJECT_NAME} CXX)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

if (MSVC)
    find_package(pybind11 REQUIRED PATHS "C:\\vcpkg\\installed\\x64-windows")
endif (MSVC)
if (UNIX)
    # set two variables below
    set(PYTHON_LIBRARY "~/anaconda3/lib/libpython3.8.so" CACHE FILEPATH "")
    set(PYTHON_EXECUTABLE "~/anaconda3/bin/python3.8" CACHE FILEPATH "")
    #find_package(pybind11 REQUIRED PATHS "~/anaconda3/share/cmake/pybind11/")
    find_package(pybind11 REQUIRED PATHS "~/anaconda3/lib/python3.8/site-packages/pybind11/share/cmake")
endif (UNIX)

add_subdirectory("./pyisocurve")
add_subdirectory("./pygeodesic")
