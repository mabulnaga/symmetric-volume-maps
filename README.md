# Symmetric Volume Maps

An algorithm for finding dense correspondences between volumes represented as tetrahedral meshes. This method can map between volumes that are close or far from isometries, and is flexible to handle landmark-based initializations, or initializations with a boundary surface map or volumetric map. There is also an interactive tool for manually selecting landmarks. This code is based on the paper "Symmetric Volume Maps" by S. Mazdak Abulnaga, Oded Stein, Polina Golland, and Justin Solomon, [eprint arXiV:2202.02568](https://arxiv.org/abs/2202.02568).

![alt text](https://github.com/mabulnaga/symmetric-volume-maps/blob/main/symmetric-volume-map-teaser.png)

## Requirements
- MATLAB
- MATLAB Package:
    - [geometry processing toolbox](https://github.com/alecjacobson/gptoolbox)
- Fast 3x3 SVD Computation:
    - Linux and Windows: [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html)
    - Mac: [ARFF toolbox](https://github.com/dpa1mer/arff)
- (Optional) GPU with CUDA compute capability

## Organization
- The main mapping algorithm is `symmetric_volume_map.m`. We have include three demos demonstrating how to run the code with different initializations.
- The `SVD` folder contains C++ code for fast 3x3 SVD computation taken from the [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) library.
- The `GPU_projection` folder contains CUDA code for tetrahedron projection used for faster convergence. The code is adopted from the code base of Li et al., ["Interactive All-Hex Meshing via Cuboid Decomposition" (2021)](https://github.com/lingxiaoli94/interactive-hex-meshing/blob/main/README.md).
- The `helpers`, `utils` and `energy` folders contain MATLAB code used in the algorithm.
- The `data` folder contains examples of volumetric meshes as `.VTK` files, manually marked landmarks, `landmarks.mat`, and initialized maps using ["Reversible Harmonic Maps between Discrete Surfaces"](https://arxiv.org/abs/1801.02453) by Ezuz et al. (2019). Tetrahedral meshes were generating using ["Fast Tetrahedral Meshing in the Wild"](https://arxiv.org/abs/1908.03581) by Hu et al. (2020); [(code)](https://github.com/wildmeshing/fTetWild).

## Installation
- SVD mex code:
    - (Linux and Windows): modify the path to your installation of Eigen in `SVD/call_mex_Eigen.m` and run ```cd SVD; call_mex_Eigen.m``` 
    - (Mac): Build `batchop_cpu` in ARFF: ```cd src/batchop; mexbuild /path/to/tbb/include```
    - If you encounter issues with Eigen, you may use `batchop_cpu` on a Windows or Linux machine. Simply modify line 15 of `helpers/compute_signed_SVD_batch.m`
- gptoolbox: install with Cmake, following these [instructions](https://github.com/alecjacobson/gptoolbox/blob/master/mex/README.md).
- GPU projection (optional): modify the path to point to your CUDA installation in `GPU_projection/mex_build.m` and run ```cd GPU_projection; mex_build.m```

## Usage
Our method can compute maps when initialized by either a coarse set of landmarks, or a surface or volume map. We require pairs of tetrahedral meshes in .VTK format. 

There are three demos provided to get started:
- demo_initialized_volume_map.m: Mapping when provided with an initial volumetric map;
- demo_landmark.m: Mapping with only a landmark-based initialization; and
- demo_interactive_landmark.m: Using our interactive tool to manually select landmarks.

We are currently working on a demo to convert a surface map initialization to a volume map.

## Development
Please contact Mazdak Abulnaga, abulnaga@mit.edu.

## Citing and Paper
If you use this method or some parts of the code, please consider citing our paper: [eprint arXiV:2202.02568](https://arxiv.org/abs/2202.02568)
```
@article{abulnaga2023symmetric,
    title = {Symmetric Volume Maps: Order-invariant Volumetric Mesh Correspondence with Free Boundary},
    author={Abulnaga, S Mazdak and Stein, Oded and Golland, Polina and Solomon, Justin},
    journal={ACM Trans. Graph.},
    volume = {42},
    number = {3},
    articleno = {25},
    year={2023}}
```
