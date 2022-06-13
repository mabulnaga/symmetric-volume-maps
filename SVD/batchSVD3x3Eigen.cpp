#include "mex.h"
#include <iostream>
#include <omp.h>
#include <cstring>
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inMatrix; // n x 3 x 3
    double *allU, *allS, *allV; // output matrices

    if(nrhs!=1) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","One input required.");
    if(nlhs!=3) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Three outputs required.");
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    
    if (mxGetNumberOfDimensions(prhs[0]) != 3 || mxGetDimensions(prhs[0])[0] != 3 ||
        mxGetDimensions(prhs[0])[1] != 3)
        mexErrMsgTxt("Input should be 3 x 3 x n");
    
    size_t n = mxGetDimensions(prhs[0])[2];
    
    inMatrix = mxGetPr(prhs[0]);
    
    /* create the output matrix */
    size_t outSize[3] = {3,3,n};
    plhs[0] = mxCreateNumericArray(3, outSize, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, outSize, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(3, outSize, mxDOUBLE_CLASS, mxREAL);
    
    allU = mxGetPr(plhs[0]);
    allS = mxGetPr(plhs[1]);
    allV = mxGetPr(plhs[2]);
    
    memset(allS, 0, 9*n*sizeof(double)); // did I do this right?
    
    #pragma omp parallel for 
    for (int i = 0; i < n; i++) {
        //mexPrintf("%g\n", omp_get_thread_num() );
        double *M = inMatrix + 9*i;
        
        Matrix3d mtx { {M[0], M[3], M[6]}, {M[1], M[4], M[7]}, {M[2], M[5], M[8]} }; // C++11 is magic
        JacobiSVD<Matrix3d,HouseholderQRPreconditioner> svd(mtx, ComputeFullU | ComputeFullV);
        
        Matrix3d UU = svd.matrixU();
        Vector3d SS = svd.singularValues();
        Matrix3d VV = svd.matrixV();
        
        double *U = allU + 9*i, *S = allS + 9*i, *V = allV + 9*i;
        U[0] = UU(0,0); U[3] = UU(0,1); U[6] = UU(0,2);
        U[1] = UU(1,0); U[4] = UU(1,1); U[7] = UU(1,2);
        U[2] = UU(2,0); U[5] = UU(2,1); U[8] = UU(2,2);
        S[0] = SS(0); // diagonal elements...
        S[4] = SS(1);
        S[8] = SS(2);
        V[0] = VV(0,0); V[3] = VV(0,1); V[6] = VV(0,2);
        V[1] = VV(1,0); V[4] = VV(1,1); V[7] = VV(1,2);
        V[2] = VV(2,0); V[5] = VV(2,1); V[8] = VV(2,2);
    }
}
