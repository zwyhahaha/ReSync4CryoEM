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
 * This is a MEX-file for MATLAB.
 *
 */



void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    mwSize ii,jj,i,j,k,n,K,K2;
    double *C, *S, *R, *B, *Re, *filter1, *filter2;
    double t,fa;
    
    double *z; 
    double a[3],rs[6]; 
    const mwSize *dims0,*dims1,*dims2,*dims3;
    mwSize     ndim0,ndim1,ndim2,ndim3;
    mwSize pos;
    
    /* Check for proper number of arguments. */
    if(nrhs!=4) {
        mexErrMsgTxt("Four inputs required.");
    } else if(nlhs>4) {
        mexErrMsgTxt("Too many output arguments.");
    }
    ndim0 = mxGetNumberOfDimensions(prhs[0]);
    dims0 = mxGetDimensions(prhs[0]);
    
    ndim1 = mxGetNumberOfDimensions(prhs[1]);
    dims1 = mxGetDimensions(prhs[1]);
    
    ndim2 = mxGetNumberOfDimensions(prhs[2]);
    dims2 = mxGetDimensions(prhs[2]);
    
    ndim3 = mxGetNumberOfDimensions(prhs[3]);
    dims3 = mxGetDimensions(prhs[3]);
    
    n     = dims0[2]; //size of C[2]
    K     = dims2[1]; //size of filter1
    K2    = dims3[1]; //size of filter2
    
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateNumericArray(ndim0, dims0, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateNumericArray(ndim1, dims1, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(ndim0, dims0, mxDOUBLE_CLASS, mxREAL);
  
    /* Assign pointers to each input and output. */
    R       = mxGetPr(prhs[0]);
    C       = mxGetPr(prhs[1]);
    filter1 = mxGetPr(prhs[2]);
    filter2 = mxGetPr(prhs[3]);
    S       = mxGetPr(plhs[0]); //subgrad
    z       = mxGetPr(plhs[1]); //function value
    B       = mxGetPr(plhs[2]);
    Re      = mxGetPr(plhs[3]); //Remannian subgrad
    
    //initial z, S, Re
    z[0] = 0; 
    pos  = 0; 
    mwSize i9,j9,i2K,j2K,i2,j2,k3;
    for (i=0; i<n; i++){
        i9 = i*9; 
        for (k=0; k<3; k++){
            S[i9+k]   = 0; 
            S[i9+k+3] = 0;
            S[i9+k+6] = 0;
            Re[i9+k]  = 0; 
            Re[i9+k+3]= 0;
            Re[i9+k+6]= 0;
        }
    }
       
    for (ii=0; ii<K; ii++){
        i   = filter1[ii]-1;
        i9  = i*9;   i2 = i*2; 
        i2K = i*2*K; 
        for (jj=0; jj<K2; jj++){
            j   = filter2[jj]-1;
            j9  = j*9;    j2 = j*2; 
            j2K = j*2*K; 
            fa = 0; 
            for (k=0; k<3; k++){ 
                a[k] = R[i9+k]*C[j2K+i2]+R[i9+k+3]*C[j2K+i2+1];
                a[k] = a[k] - R[j9+k]*C[i2K+j2] - R[j9+k+3]*C[i2K+j2+1];
                fa  +=  a[k]*a[k]; 
                
                B[j2K+i2+k] = a[k];
            }
            fa    = sqrt(fa); 
            z[0] += fa;
            fa    = sqrt(fa * fa + 1e-16); 
            for (k=0; k<3; k++){                
                S[i9+k] = S[i9+k] + a[k] * C[j2K+i2]/fa; 
                S[i9+k+3] = S[i9+k+3] + a[k] * C[j2K+i2+1]/fa; 
            }
        }
    }
    
    for (ii=0; ii<K; ii++){
        i  = filter1[ii]-1;
        i9 = i*9;
        for (k=0; k<3; k++){
            k3      = 3*k;
            rs[k]   = R[i9+k3]*S[i9]+R[i9+k3+1]*S[i9+1]+R[i9+k3+2]*S[i9+2];
            rs[k+3] = R[i9+k3]*S[i9+3]+R[i9+k3+1]*S[i9+1+3]+R[i9+k3+2]*S[i9+2+3];
        }
        for (k=0; k<3; k++){
            Re[i9+k]   = Re[i9+k]+R[i9+k+3]*rs[1]+R[i9+k+6]*rs[2]-R[i9+k+3]*rs[3];
            Re[i9+k+3] = Re[i9+k+3]+R[i9+k]*rs[3]-R[i9+k]*rs[1]+R[i9+k+6]*rs[5];
            Re[i9+k+6] = Re[i9+k+6]-R[i9+k]*rs[2]-R[i9+k+3]*rs[5];
        }
    }
}
