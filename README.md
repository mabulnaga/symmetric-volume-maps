# Symmetric-Volume-Maps

An algorithm for finding correspondences between volumes represented as tetrahedral meshes. 

![alt text](https://github.com/mabulnaga/placenta-flattening/blob/master/flattening_flowchart.png)

### Requirements
- MATLAB
- MATLAB Packages:
    - [geometry processing toolbox](https://github.com/alecjacobson/gptoolbox)
- [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) for fast 3x3 SVD computation
- (Optional) GPU with CUDA compute capability

### Installation
To install the fast 3x3 SVD mex code, modify "SVD/call_mex_Eigen.m" and point to your installation of Eigen. Note that we've encountered issues with this library in Mac. If that applies to you, instead install the MATLAB package [ARFF](https://github.com/dpa1mer/arff), and compile "batchop_cpu". There is a check for Mac OS in the appropriate file, "helpers/compute_signed_SVD_batch.m". You may modify this file if wanting to use batchop instead of Eigen.

If using a GPU, mexbuild the mexcuda file in "GPU_projection/mexbuild.m". You will need to modify this file to point to your CUDA library.

### Usage
There are three demos provided.
- demo_initialized_volume_map.m: A demo when already provided with an initial volumetric map;
- demo_landmark.m: Mapping with only a landmark-based initialization; and
- demo_interactive_landmark.m: Using our interactive tool to manually select landmarks.

We are currently working on a demo to convert a surface map initialization to a volume map.

### Development
Please contact Mazdak Abulnaga, abulnaga@mit.edu.

### Citing and Paper
If you use this method or some parts of the code, please consider citing our paper: [eprint arXiV:2202.02568](https://arxiv.org/abs/2202.02568)
```
@article{abulnaga2022symmetric,
    title={Symmetric Volume Maps},
    author={Abulnaga, S Mazdak and Stein, Oded and Golland, Polina and Solomon, Justin},
    journal={arXiv preprint arXiv:2202.02568},
    year={2022}}
```
