cmake_minimum_required(VERSION 3.4...3.18)

set(PROJECT_NAME pyremesh)
project (${PROJECT_NAME} CXX)

set(SRC_FILES
    pyremesh.cpp
    IncrementalRemesher.h
    IncrementalRemesher.cpp
    PolylineNetwork.h
    PolylineNetwork.cpp
    PolylineBranch.h
    PolylineBranch.cpp
)

pybind11_add_module(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/eigen/")
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/OpenMesh-9.0.0/")
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/libigl/include/")
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCore)

if (MSVC)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
endif (MSVC)

if (UNIX)
    #find_package(Threads REQUIRED)
endif (UNIX)
