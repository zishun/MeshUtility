cmake_minimum_required(VERSION 3.4...3.18)

set(PROJECT_NAME OpenMeshCore)
project (${PROJECT_NAME} CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SRC_FILES
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/Config.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/EigenVectorT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/LoopSchemeMaskT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/MathDefs.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/NormalConeT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/NormalConeT_impl.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/Plane3d.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/QuadricT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/Vector11T.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/VectorT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Geometry/VectorT_inc.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/ArrayItems.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/ArrayKernel.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/ArrayKernelT_impl.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/AttribKernelT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/Attributes.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/BaseKernel.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/BaseMesh.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/Casts.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/CirculatorsT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/DefaultPolyMesh.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/DefaultTriMesh.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/FinalMeshItemsT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/Handles.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/IteratorsT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/PolyConnectivity.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/PolyConnectivity_inline_impl.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/PolyMeshT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/PolyMeshT_impl.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/SmartHandles.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/SmartRange.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/Status.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/Tags.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/Traits.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/TriConnectivity.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/TriMeshT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/TriMeshT_impl.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/compiler.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/config.h
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/config.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/mostream.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/omstream.hh
    # ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/OpenMeshDLLMacros.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/AutoPropertyHandleT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/BaseProperty.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/color_cast.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/Endian.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/GenProg.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/HandleToPropHandle.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/Noncopyable.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/Property.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/PropertyContainer.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/PropertyManager.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/RandomNumberGenerator.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/SingletonT.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/SingletonT_impl.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/typename.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/vector_cast.hh
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/vector_traits.hh

    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/ArrayKernel.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/BaseKernel.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/PolyConnectivity.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Mesh/TriConnectivity.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/System/omstream.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/BaseProperty.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/Endian.cc
    ../../ext/OpenMesh-9.0.0/OpenMesh/Core/Utils/RandomNumberGenerator.cc
)

add_library(${PROJECT_NAME} STATIC
            ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE "../../ext/OpenMesh-9.0.0/")

if (MSVC)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
    add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
endif (MSVC)
if (UNIX)
endif (UNIX)
