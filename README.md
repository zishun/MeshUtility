[![Build](https://github.com/zishun/MeshUtility/actions/workflows/wheels.yml/badge.svg?branch=main)](https://github.com/zishun/MeshUtility/actions/workflows/wheels.yml)
# MeshUtility

> A collection of python utilities for mesh processing.

Install with
```
pip install meshutility==0.0.1
```

A simple tutorial is available [here](https://zishun.github.io/projects/MeshUtility/).


## Functions & Modules
- ```colormap_vertex_color```: assign vertex color to visualize a scalar field defined on mesh.
- ```IsoCurve``` module: extract isocurves on a scalar field defined on a manifold triangular mesh.
- ```mesh_cut```: cut a mesh along a vertex chain.
- ```mesh_split```: split a mesh by inserting new vertices defined on mesh edges.
- ```read_obj_lines```: read polyline from a [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements).
- ```sphere_cvt```: iteratively approximate centroidal Voronoi tessellation (CVT) on the unit sphere (kind of uniform sampling).
- ```split_connected_components```: split connected components.
- ```write_obj_lines```: write polyline as a Wavefront .obj file that can be open with MeshLab.
- ```pygeodesic``` module (C++): exact geodesic for triangular meshes.
- ```pyisocurve``` module (C++): almost the same as ```IsoCurve``` above.
- ```pyremesh``` module (C++): incremental isotropic remeshing.
- ```pyshapeop``` module (C++): a partial binding of ShapeOp.


## Build from Source

see [```build.md```](https://github.com/zishun/MeshUtility/blob/main/build.md)
