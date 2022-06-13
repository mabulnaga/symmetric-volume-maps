#include "mex.h"
#include "matrix.h"
#include "gpu/mxGPUArray.h"
#include "tet_proj.cuh"

static bool mexInitialized;
static mxArray * gpuCanary;
static cudaStream_t stream;
static void uninit();
static char const * const errInputId = "qtimes:InvalidInput";
static char const * const errCudaId = "qtimes:CudaError";
static char const * const errCudaMsg = "qtimes encountered a CUDA error.";

//initialize stream
static bool init() {
    if (!mexInitialized || !mxGPUIsValidGPUData(gpuCanary)) {
        // Initialize the MATLAB GPU API if not already initialized.
        if (mxInitGPU() != MX_GPU_SUCCESS) {
            return false;
        }

        const size_t one = 1;
        cudaStreamCreate(&stream); // Note that a pointer must be passed to `cudaCreateStream`.

        mxGPUArray * canary = mxGPUCreateGPUArray(1, &one, mxDOUBLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
        gpuCanary = mxGPUCreateMxArrayOnGPU(canary);
        mxGPUDestroyGPUArray(canary);
        mexMakeArrayPersistent(gpuCanary);
        mexInitialized = true;
        mexAtExit(uninit);
    }
    return true;
}

//destroy stream
static void uninit() {
    if (mexInitialized) {
        mxDestroyArray(gpuCanary);
        cudaStreamDestroy(stream); // Note that a value, not a pointer, is passed to `cudaDestroyStream`.
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{   //inputs:
    // points: dim x n_p
    // v0: dim x n_t
    // e_mat:  dim x 3xn_t (vertex i - vertex 0)
    // e_dot_mat: 3x3xn_t
    // e_dot_inv_mat: 3x3xn_t
    // is_degenerate: n_tx1
    // Initialize solver
    if (!init()) {
        mexErrMsgIdAndTxt(errCudaId, errCudaMsg);
    }

    if (nrhs < 2) {
        mexErrMsgIdAndTxt(errInputId, "Must provide an operation and input matrix.");
    }


	// Throw an error if the 1st input is not a GPU array.
    if (!mxIsGPUArray(prhs[0]) || !mxGPUIsValidGPUData(prhs[0])) {
        mexErrMsgIdAndTxt(errInputId, "Inputs to tet proj must be of type double gpuArray.");
    }

    // Unwrap input to an mxGPUArray (must be real double). Repeat for e_mat, E_dot_mat, e_dot_inv_mat
    const mxGPUArray *points = mxGPUCreateFromMxArray(prhs[0]);
    const mxGPUArray *v0 = mxGPUCreateFromMxArray(prhs[1]);
    const mxGPUArray *e_mat = mxGPUCreateFromMxArray(prhs[2]);
    const mxGPUArray *e_dot_mat = mxGPUCreateFromMxArray(prhs[3]);
    const mxGPUArray *e_dot_inv_mat = mxGPUCreateFromMxArray(prhs[4]);
    const mxGPUArray *is_degenerate =mxGPUCreateFromMxArray(prhs[5]);
    // do this for all of the inputs
    if (mxGPUGetClassID(v0) != mxDOUBLE_CLASS || mxGPUGetComplexity(v0) != mxREAL) {
        mexErrMsgIdAndTxt(errInputId, "Inputs to tet proj must be of type double gpuArray.");
    }
    //repeat for each input
    size_t nDimv0 = mxGPUGetNumberOfDimensions(v0);
    if (nDimv0 != 2) {
        mexErrMsgIdAndTxt(errInputId, "batchop operates on 2D or 3D arrays only.");
    }
    const size_t * dim_points = mxGPUGetDimensions(points);
    const size_t * dim_v0 = mxGPUGetDimensions(v0);
    const size_t num_tets = dim_v0[1];
    const size_t num_points = dim_points[1];
    const size_t dim = dim_v0[0];
    const size_t * dim_e_mat = mxGPUGetDimensions(e_mat);

    if(dim_e_mat[2] != num_tets){
        mexErrMsgIdAndTxt(errInputId, "3D matrix has incorrect dimensions.");
    }
    
    size_t num_threads = 128;
    dim3 num_blocks(num_points);
    size_t shared_size = num_threads * (sizeof(double) + sizeof(size_t));

    // create outputs. Repeat for all. Use int for the outputindices. mxINT32_CLASS
    size_t dimResultDists[] = {num_points, 1};
    size_t dimResultWeights[] = {num_points,3};
    mxGPUArray * result_dists = mxGPUCreateGPUArray(2,dimResultDists, mxDOUBLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray * result_idxs = mxGPUCreateGPUArray(2,dimResultDists, mxINT32_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    mxGPUArray * result_weights = mxGPUCreateGPUArray(2,dimResultWeights, mxDOUBLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);

    // call function
    // TODO: FIGURE OUT A BETTER WAY TO INITIALIZE THE DIMENSION (NOT HARD CODE 6). add stream back in.
    tet_kernels::GeneralizedTetrahedronProjectionKernel<6, double>
                    <<<num_blocks, num_threads, shared_size, stream>>>(
                            static_cast<const double *>(mxGPUGetDataReadOnly(points)),
                            num_tets,
                            static_cast<const double *>(mxGPUGetDataReadOnly(v0)),
                            static_cast<const double *>(mxGPUGetDataReadOnly(e_mat)),
                            static_cast<const double *>(mxGPUGetDataReadOnly(e_dot_mat)),
                            static_cast<const double *>(mxGPUGetDataReadOnly(e_dot_inv_mat)),
                            static_cast<const double*>(mxGPUGetDataReadOnly(is_degenerate)),
                            static_cast<double *>(mxGPUGetData(result_dists)), 
                            static_cast<int *>(mxGPUGetData(result_idxs)),
                            static_cast<double *>(mxGPUGetData(result_weights)));
                        

    plhs[0] = mxGPUCreateMxArrayOnGPU(result_dists);
    plhs[1] = mxGPUCreateMxArrayOnGPU(result_idxs);
    plhs[2] = mxGPUCreateMxArrayOnGPU(result_weights);
    mxGPUDestroyGPUArray(result_dists);
    mxGPUDestroyGPUArray(result_idxs);
    mxGPUDestroyGPUArray(result_weights);

//     mxFree((void *)dim_points);
//     mxFree((void *)nDimv0);
//     mxFree((void *)dim_e_mat);
    mxGPUDestroyGPUArray(points);
    mxGPUDestroyGPUArray(v0);
    mxGPUDestroyGPUArray(e_mat);
    mxGPUDestroyGPUArray(e_dot_mat);
    mxGPUDestroyGPUArray(e_dot_inv_mat);
    mxGPUDestroyGPUArray(is_degenerate);

    
}