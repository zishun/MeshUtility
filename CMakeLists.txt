cmake_minimum_required(VERSION 3.4...3.18)
project (MeshUtilityCPP CXX)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory("./ext/pybind11")

add_subdirectory("./src/OpenMesh")
add_subdirectory("./src/pyisocurve")
add_subdirectory("./src/pygeodesic")
add_subdirectory("./src/pyremesh")
add_subdirectory("./src/pyshapeop")
