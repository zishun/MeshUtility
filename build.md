# Setup MeshUtility

## TL;DR
```
pip install numpy scipy networkx matplotlib openmesh
git clone --recurse-submodules https://github.com/zishun/MeshUtility.git
python setup.py install
```

## Dependencies

<details>
<summary>
<b>
Python:
</b>
</summary>

* numpy
* scipy
* networkx
* matplotlib
* openmesh

</details>


Install with 
```shell
pip install numpy scipy networkx matplotlib openmesh
```

<details>
<summary>
<b>
C++:
</b>
</summary>
* Eigen
* OpenMesh
* libigl
* geodesic: Danil Kirsanov's implementation of MMP algorithm.
* ShapeOp
</details>

All have been included here directly or as submodules.

## Build
1. Clone the repo: ```git clone --recurse-submodules https://github.com/zishun/MeshUtility.git```
2. ```python setup.py install```
