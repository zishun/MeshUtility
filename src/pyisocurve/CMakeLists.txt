cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME pyisocurve)
project (${PROJECT_NAME} CXX)

set(SRC_FILES
    main.cpp
    isocurve.h
    isocurve.cpp
)

pybind11_add_module(${PROJECT_NAME} ${SRC_FILES})

if (MSVC)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
    target_include_directories(${PROJECT_NAME} PRIVATE
        C:/vcpkg/installed/x64-windows/include/
        C:/vcpkg/installed/x64-windows/include/eigen3/
    )
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
        target_link_libraries(${PROJECT_NAME} PUBLIC
            C:/vcpkg/installed/x64-windows/lib/OpenMeshCore.lib
            )
    endif (CMAKE_BUILD_TYPE MATCHES "Release")
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        target_link_libraries(${PROJECT_NAME} PUBLIC
            C:/vcpkg/installed/x64-windows/debug/lib/OpenMeshCored.lib
            )
    endif (CMAKE_BUILD_TYPE MATCHES "Debug")
endif (MSVC)
if (UNIX)
    #find_package(Threads REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE
        /usr/include/eigen3/
        )
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCore)
endif (UNIX)
