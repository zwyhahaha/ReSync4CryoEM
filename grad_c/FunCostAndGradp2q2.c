#include <mex.h>
#include <stdio.h>
#include <math.h>
/*
 *
 * Computational function that update .
 **********************************************************
 * function [GradRi,ObjFun] = GradObjFun(Rmat, C)
K      = size(Rmat,3);
GradRi = zeros(3,3,K);
ObjFun = 0;
% tic;
for i =1:K
    for j = 1:K
        tmp = (Rmat(:,:,i) *C(:,i,j)-Rmat(:,:,j) *C(:,j,i));
        GradRi(:,:,i) = GradRi(:,:,i) + tmp *C(:,i,j)';
        ObjFun = ObjFun + (norm(tmp,2))^2;
    end
end
 **********************************************************
 * This is a MEX-file for MATLAB. **/

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    mwSize i,j,k,K;
    double *C, *S, *R, *B;
    double t;
    
    double *z; 
    double a[3]; 
    const mwSize *dims0,*dims1;
    mwSize     ndim0,ndim1;
    mwSize pos; 
    
    /* Check for proper number of arguments. */
    if(nrhs!=2) {
        mexErrMsgTxt("Two inputs required.");
    } else if(nlhs>3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    ndim0 = mxGetNumberOfDimensions(prhs[0]);
    dims0 = mxGetDimensions(prhs[0]);
    K     = dims0[2];
    
    ndim1 = mxGetNumberOfDimensions(prhs[1]);
    dims1 = mxGetDimensions(prhs[1]);
    
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateNumericArray(ndim0, dims0, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateNumericArray(ndim1, dims1, mxDOUBLE_CLASS, mxREAL);
  
    /* Assign pointers to each input and output. */
    R       = mxGetPr(prhs[0]);
    C       = mxGetPr(prhs[1]);
    S       = mxGetPr(plhs[0]);
    z       = mxGetPr(plhs[1]);
    B       = mxGetPr(plhs[2]);
    
    // initial z and S
    z[0] = 0; 
    pos  = 0; 
    mwSize i6,j6,i2K,j2K,i2,j2;
    for (i=0; i<K; i++){
        i6 = i*6; 
        for (k=0; k<3; k++){
            S[i6+k]   = 0; 
            S[i6+k+3] = 0;
        }
    }
                    
    for (i=0; i<K; i++){
        i6  = i*6;   i2 = i*2; 
        i2K = i*2*K; 
        for (j=0; j<K; j++){
            j6  = j*6;    j2 = j*2; 
            j2K = j*2*K; 
            for (k=0; k<3; k++){ 
                a[k] = R[i6+k]*C[j2K+i2]+R[i6+k+3]*C[j2K+i2+1];
                a[k] = a[k] - R[j6+k]*C[i2K+j2] - R[j6+k+3]*C[i2K+j2+1];
                z[0]    +=  a[k]*a[k]; 
                S[i6+k] = S[i6+k] + a[k] * C[j2K+i2]; 
                S[i6+k+3] = S[i6+k+3] + a[k] * C[j2K+i2+1];
                B[j2K+i2+k] = a[k];
            }
        }
    }
}
