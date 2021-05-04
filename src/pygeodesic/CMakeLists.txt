cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME pygeodesic)
project (${PROJECT_NAME} CXX)

set(SRC_FILES
    main.cpp
    mmp_wrapper.h
    mmp_wrapper.cpp
)

pybind11_add_module(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/geodesic/src/")

if (MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
    target_include_directories(${PROJECT_NAME} PRIVATE
       C:/vcpkg/installed/x64-windows/include
       C:/vcpkg/installed/x64-windows/include/eigen3/
    )
    if (CMAKE_BUILD_TYPE MATCHES "Release")
       target_link_libraries(${PROJECT_NAME} PUBLIC
           C:/vcpkg/installed/x64-windows/lib/OpenMeshCore.lib
       )
    endif (CMAKE_BUILD_TYPE MATCHES "Release")
endif (MSVC)

if (UNIX)
    target_include_directories(${PROJECT_NAME} PRIVATE
       /usr/include/eigen3/
    )
    if (CMAKE_BUILD_TYPE MATCHES "Release")
       target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCore)
    endif (CMAKE_BUILD_TYPE MATCHES "Release")
endif (UNIX)
