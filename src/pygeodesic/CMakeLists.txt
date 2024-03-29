cmake_minimum_required(VERSION 3.4...3.18)

set(PROJECT_NAME pygeodesic)
project (${PROJECT_NAME} CXX)

set(SRC_FILES
    main.cpp
    mmp_wrapper.h
    mmp_wrapper.cpp
    fm.h
    fm.cpp
)

pybind11_add_module(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/geodesic/src/")
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/eigen/")
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/OpenMesh-9.0.0/")
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCore)

if (MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
endif (MSVC)

if (UNIX)
    #if (CMAKE_BUILD_TYPE MATCHES "Release")
    #endif (CMAKE_BUILD_TYPE MATCHES "Release")
endif (UNIX)
