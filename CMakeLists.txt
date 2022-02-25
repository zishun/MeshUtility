cmake_minimum_required(VERSION 3.4...3.18)
project (MeshUtilityCPP CXX)

add_subdirectory("./ext/pybind11")

add_subdirectory("./src/OpenMesh")
add_subdirectory("./src/pyisocurve")
add_subdirectory("./src/pygeodesic")
add_subdirectory("./src/pyremesh")
add_subdirectory("./src/pyshapeop")
