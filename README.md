# Mesh Utility

A collection of python utilities for mesh processing.


## Functions
- ```colormap_vertex_color```: assign vertex color to visualize a scalar field defined on mesh.
- ```isocurve```, ```pyisocurve```: extract isocurves on a scalar field defined on a manifold triangular mesh.
- ```mesh_cut```: cut a mesh along a vertex chain.
- ```mesh_split```: split a mesh by inserting new vertices defined on mesh edges.
- ```sphere_cvt``` iteratively approximate centroidal Voronoi tessellation (CVT) on the unit sphere.
- ```write_obj_lines```: write polyline as a [Wavefront_.obj_file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements), which can be open with MeshLab.

## Build
1. Clone the repo: ```git clone --recurse-submodules https://github.com/zishun/MeshUtility.git```
2. Install required C++ libraries
    - pybind11
    - eigen3
    - openmesh
3. Build C++ modules using ```src\CMakeLists.txt```.
    - On Windows, you may copy ```OpenMeshCore.dll``` to ```src\lib\```.
4. Install required Python modules. see ```requirements.txt```

