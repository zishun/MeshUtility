[![Build](https://github.com/zishun/MeshUtility/actions/workflows/wheels.yml/badge.svg?branch=main)](https://github.com/zishun/MeshUtility/actions/workflows/wheels.yml)
# MeshUtility

> A collection of python utilities for mesh processing.

Install with
```
pip install meshutility==0.0.2
```

A simple tutorial is available [here](https://zishun.github.io/projects/MeshUtility/) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/zishun/meshutility/blob/main/tests/meshutility-tutorial.ipynb).


## Functions & Modules
- ```colormap_vertex_color```: assign vertex color to visualize a scalar field defined on mesh.
- ```cut_along_curve```: cut a mesh along a vertex chain.
- ```ff_graph```: face-face graph of mesh.
- ```get_scalar_field_on_resampled_points```: given a scalar field defined on a triangular mesh, get field values on resampled points.
- ```pygeodesic``` module (C++): geodesic for triangular meshes.
    - exact geodesic by Danil Kirsanov.
    - fast marching, may use different speed on each face.
- ```pyisocurve``` module (C++): extract isocurves on a scalar field defined on a manifold triangular mesh.
- ```pyremesh``` module (C++): incremental isotropic remeshing.
- ```pyshapeop``` module (C++): a partial binding of ShapeOp.
- ```read_obj_lines```: read polyline from a [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements).
- ```remove_unreferenced_vertices```: remove unreferenced vertices.
- ```sphere_cvt```: iteratively approximate centroidal Voronoi tessellation (CVT) on the unit sphere (kind of uniform sampling).
- ```split_connected_components```: split connected components.
- ```split_mesh, split_mesh_complete```: split a mesh by inserting new vertices defined on mesh edges.
- ```write_obj_lines```: write polyline as a Wavefront .obj file that can be open with MeshLab.
- ```vv_graph```: vertex-vertex graph of mesh.


## Build from Source

see [```build.md```](https://github.com/zishun/MeshUtility/blob/main/build.md)
