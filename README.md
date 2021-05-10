# Mesh Utility

A collection of python utilities for mesh processing.


## Functions
- ```colormap_vertex_color```: assign vertex color to visualize a scalar field defined on mesh.
- ```isocurve```: extract isocurves on a scalar field defined on a manifold triangular mesh.
- ```mesh_cut```: cut a mesh along a vertex chain.
- ```mesh_split```: split a mesh by inserting new vertices defined on mesh edges.
- ```sphere_cvt```: iteratively approximate centroidal Voronoi tessellation (CVT) on the unit sphere (kind of uniform sampling).
- ```write_obj_lines```: write polyline as a [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements), which can be open with MeshLab.
- ```pyisocurve``` (C++ based): almost the same as ```isocurve``` above.
- ```pygeodesic``` (C++ based): exact geodesic for triangular meshes.

## Build
1. Clone the repo: ```git clone --recurse-submodules https://github.com/zishun/MeshUtility.git```
2. Build C++ modules. if you do not need them (see the list above), set ```USE_CPP = False``` in [```__init__.py```](https://github.com/zishun/MeshUtility/blob/main/__init__.py#L1) and skip this step.
    1. Install required C++ libraries
        - pybind11
        - eigen3
        - openmesh
    2. Build using ```src\CMakeLists.txt```.
        - On Windows, you may copy ```OpenMeshCore.dll``` to ```src\lib\```.
3. Install required Python modules. See ```requirements.txt```.
