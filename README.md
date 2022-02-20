# Mesh Utility

A collection of python utilities for mesh processing.

A simple tutorial is available [here](https://zishun.github.io/projects/MeshUtility/).


## Functions
- ```colormap_vertex_color```: assign vertex color to visualize a scalar field defined on mesh.
- ```isocurve```: extract isocurves on a scalar field defined on a manifold triangular mesh.
- ```mesh_cut```: cut a mesh along a vertex chain.
- ```mesh_split```: split a mesh by inserting new vertices defined on mesh edges.
- ```sphere_cvt```: iteratively approximate centroidal Voronoi tessellation (CVT) on the unit sphere (kind of uniform sampling).
- ```split_connected_components```: split connected components.
- ```write_obj_lines```: write polyline as a [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements), which can be open with MeshLab.
- ```pyisocurve``` (C++ based): almost the same as ```isocurve``` above.
- ```pygeodesic``` (C++ based): exact geodesic for triangular meshes.


## Dependencies

Python:
* numpy: pip
* openmesh: pip
* networkx: pip
* matplotlib: pip
* scipy: pip
* meshplot: only for jupyter notebook 
    * ```conda install -c conda-forge meshplot```

C++:
* eigen: submodule.
* OpenMesh: included.
* libigl: submodule.
* geodesic: submodule. Danil Kirsanov's implementation of MMP algorithm.


## Build
1. Clone the repo: ```git clone --recurse-submodules https://github.com/zishun/MeshUtility.git```
2. ```python setup.py install```
