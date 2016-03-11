# Projective Dynamics experiments.


## Selecting a Backend

- Eigen3: Simply run `make`

- CUDA: Run `BACKEND=cuda make`

- OpenCL: Run `BACKEND=opencl make`

The CUDA and OpenCL backends require [ViennaCL](http://viennacl.sourceforge.net/) to be installed.


## Blender Export

The file `blender/export.py` provides an add-on for exporting *triangular* meshes from Blender
to .json file. This file can be loaded in projective dynamics by running
`./pd 0 0 10 filename.json`.

The spring constraints are formed by edges of the triangle mesh and attachment constraints come
from the vertex group selected for pinning the cloth in Blender cloth simulation properties.
