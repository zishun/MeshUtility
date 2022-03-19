cmake_minimum_required(VERSION 3.4...3.18)

set(PROJECT_NAME utils)
project (${PROJECT_NAME} CXX)

set(SRC_FILES
    utils.cpp
)

pybind11_add_module(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/eigen/")
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/libigl/include/")
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCore)

if (MSVC)
endif (MSVC)

if (UNIX)
endif (UNIX)
