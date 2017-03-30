/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _GpuConv1D_cu
#define _GpuConv1D_cu

#include <stdio.h>

#include "GaussFunction.h"
#include "ScalarRadialKernel.h"

#include "GpuConv1D.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Conv1D :
////////////////////////////////////////////////////////////////////////////////////////////////////

// Thread kernel: computation of \f$ \gamma_i = \sum_j K(x_i,y_j)\beta_j for index i given by thread id.
template < typename TYPE, int DIMPOINT, int DIMVECT, class KER  >
__global__ void GpuConv1DOnDevice(KER Ker,
                                      TYPE *x, TYPE *y, TYPE *beta, TYPE *gamma,
                                      int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ char SharedData_char[];
    TYPE* const SharedData = reinterpret_cast<TYPE*>(SharedData_char);

    TYPE xi[DIMPOINT], gammai[DIMVECT];
    if(i<nx)  // we will compute gammai only if i is in the range
    {
        // load xi from device global memory
        for(int k=0; k<DIMPOINT; k++)
            xi[k] = x[i*DIMPOINT+k];
        for(int k=0; k<DIMVECT; k++)
            gammai[k] = 0.0f;
    }

    for(int jstart = 0, tile = 0; jstart < ny; jstart += blockDim.x, tile++)
    {
        int j = tile * blockDim.x + threadIdx.x;
        if(j<ny) // we load yj and betaj from device global memory only if j<ny
        {
            int inc = DIMPOINT + DIMVECT;
            for(int k=0; k<DIMPOINT; k++)
                SharedData[threadIdx.x*inc+k] = y[j*DIMPOINT+k];
            for(int k=0; k<DIMVECT; k++)
                SharedData[threadIdx.x*inc+DIMPOINT+k] = beta[j*DIMVECT+k];
        }
        __syncthreads();
        
        if(i<nx) // we compute gammai only if needed
        {
            TYPE *yj, *betaj;
            yj = SharedData;
            betaj = SharedData + DIMPOINT;
            int inc = DIMPOINT + DIMVECT;
            for(int jrel = 0; jrel < blockDim.x && jrel<ny-jstart; jrel++, yj+=inc, betaj+=inc)
		Ker.Eval(gammai,xi,yj,betaj);
    	}
	__syncthreads();
    }

    // Save the result in global memory.
    if(i<nx)
        for(int k=0; k<DIMVECT; k++)
            gamma[i*DIMVECT+k] = gammai[k];
}



template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuEvalConv1D(KER Ker, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny)
{
    // Data on the device.
    TYPE* x_d;
    TYPE* y_d;
    TYPE* beta_d;
    TYPE* gamma_d;

    // Allocate arrays on device.
    cudaMalloc((void**)&x_d, sizeof(TYPE)*(nx*DIMPOINT));
    cudaMalloc((void**)&y_d, sizeof(TYPE)*(ny*DIMPOINT));
    cudaMalloc((void**)&beta_d, sizeof(TYPE)*(ny*DIMVECT));
    cudaMalloc((void**)&gamma_d, sizeof(TYPE)*(nx*DIMVECT));

    // Send data from host to device.
    cudaMemcpy(x_d, x_h, sizeof(TYPE)*(nx*DIMPOINT), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y_h, sizeof(TYPE)*(ny*DIMPOINT), cudaMemcpyHostToDevice);
    cudaMemcpy(beta_d, beta_h, sizeof(TYPE)*(ny*DIMVECT), cudaMemcpyHostToDevice);

    // Compute on device.
    dim3 blockSize;
    blockSize.x = 192; // number of threads in each block
    dim3 gridSize;
    gridSize.x =  nx / blockSize.x + (nx%blockSize.x==0 ? 0 : 1);
	
	GpuConv1DOnDevice<TYPE,DIMPOINT,DIMVECT,KER>
		<<<gridSize,blockSize,blockSize.x*(DIMVECT+DIMPOINT)*sizeof(TYPE)>>>
			(Ker, x_d, y_d, beta_d, gamma_d, nx, ny);

    // block until the device has completed
    cudaThreadSynchronize();

    // Send data from device to host.
    cudaMemcpy(gamma_h, gamma_d, sizeof(TYPE)*(nx*DIMVECT),cudaMemcpyDeviceToHost);

    // Free memory.
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(beta_d);
    cudaFree(gamma_d);

    return 0;
}


template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuEvalConv1D(TYPE sigma, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny)
{
	
	return GpuEvalConv1D < TYPE, DIMPOINT, DIMVECT, ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> > >
		(ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> >(GaussFunction<TYPE>(sigma)),
			x_h, y_h, beta_h, gamma_h, nx, ny);
}





////////////////////////////////////////////////////////////////////////////////////////////////////
// Grad1 Conv1D :
////////////////////////////////////////////////////////////////////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
__global__ void GpuGrad1Conv1DOnDevice(KER Ker,
        TYPE *alpha, TYPE *x, TYPE *y, TYPE *beta, TYPE *gamma,
        int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ char SharedData_char[];
    TYPE* const SharedData = reinterpret_cast<TYPE*>(SharedData_char);

    TYPE xi[DIMPOINT], alphai[DIMVECT], gammai[DIMPOINT];
    if(i<nx)  // we will compute gammai only if i is in the range
    {
        // load xi and alphai from device global memory
        for(int k=0; k<DIMPOINT; k++)
            xi[k] = x[i*DIMPOINT+k];
        for(int k=0; k<DIMVECT; k++)
            alphai[k] = alpha[i*DIMVECT+k];
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] = 0.0f;
    }

    for(int jstart = 0, tile = 0; jstart < ny; jstart += blockDim.x, tile++)
    {
        int j = tile * blockDim.x + threadIdx.x;
        if(j<ny) // we load yj and betaj from device global memory only if j<ny
        {
            int inc = DIMPOINT + DIMVECT;
            for(int k=0; k<DIMPOINT; k++)
                SharedData[threadIdx.x*inc+k] = y[j*DIMPOINT+k];
            for(int k=0; k<DIMVECT; k++)
                SharedData[threadIdx.x*inc+DIMPOINT+k] = beta[j*DIMVECT+k];
        }
        __syncthreads();
        if(i<nx) // we compute gammai only if i is in the range
        {
            TYPE *yj, *betaj;
            yj = SharedData;
            betaj = SharedData + DIMPOINT;
            int inc = DIMPOINT + DIMVECT;
            for(int jrel = 0; jrel < blockDim.x && jrel<ny-jstart; jrel++, yj+=inc, betaj+=inc)
	            Ker.Grad1(gammai,alphai,xi,yj,betaj);
        }
        __syncthreads();
    }

    // Save the result in global memory.
    if(i<nx)
        for(int k=0; k<DIMPOINT; k++)
            gamma[i*DIMPOINT+k] = gammai[k];
}



template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGrad1Conv1D(KER Ker, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny)
{

    // Data on the device.
    TYPE* x_d;
    TYPE* y_d;
    TYPE* alpha_d;
    TYPE* gamma_d;
    TYPE* beta_d;

    // Allocate arrays on device.
    cudaMalloc((void**)&x_d, sizeof(TYPE)*(nx*DIMPOINT));
    cudaMalloc((void**)&y_d, sizeof(TYPE)*(ny*DIMPOINT));
    cudaMalloc((void**)&alpha_d, sizeof(TYPE)*(nx*DIMVECT));
    cudaMalloc((void**)&beta_d, sizeof(TYPE)*(ny*DIMVECT));
    cudaMalloc((void**)&gamma_d, sizeof(TYPE)*(nx*DIMPOINT));

    // Send data from host to device.
    cudaMemcpy(x_d, x_h, sizeof(TYPE)*(nx*DIMPOINT), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y_h, sizeof(TYPE)*(ny*DIMPOINT), cudaMemcpyHostToDevice);
    cudaMemcpy(alpha_d, alpha_h, sizeof(TYPE)*(nx*DIMVECT), cudaMemcpyHostToDevice);
    cudaMemcpy(beta_d, beta_h, sizeof(TYPE)*(ny*DIMVECT), cudaMemcpyHostToDevice);

    // compute on device.
    dim3 blockSize;
    blockSize.x = 192; // number of threads in each block
    dim3 gridSize;
    gridSize.x =  nx / blockSize.x + (nx%blockSize.x==0 ? 0 : 1);

    GpuGrad1Conv1DOnDevice<TYPE,DIMPOINT,DIMVECT,KER>
		<<<gridSize,blockSize,blockSize.x*(DIMPOINT+DIMVECT)*sizeof(TYPE)>>>
			(Ker, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);

    // block until the device has completed
    cudaThreadSynchronize();

    // Send data from device to host.
    cudaMemcpy(gamma_h, gamma_d, sizeof(TYPE)*(nx*DIMPOINT),cudaMemcpyDeviceToHost);

    // Free memory.
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(alpha_d);
    cudaFree(gamma_d);
    cudaFree(beta_d);

    return 0;
}



template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGrad1Conv1D(TYPE sigma, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny)
{
	return GpuGrad1Conv1D < TYPE, DIMPOINT, DIMVECT, ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> > >
		(ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> >(GaussFunction<TYPE>(sigma)),
			alpha_h, x_h, y_h, beta_h, gamma_h, nx, ny);
}






////////////////////////////////////////////////////////////////////////////////////////////////////
// GradDiff Conv1D :
////////////////////////////////////////////////////////////////////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
__global__ void GpuGradDiffConv1DOnDevice(KER Ker,
        TYPE *x, TYPE *beta, TYPE *eta, TYPE *gamma,
        int nx)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ char SharedData_char[];
    TYPE* const SharedData = reinterpret_cast<TYPE*>(SharedData_char);

    TYPE xi[DIMPOINT], betai[DIMVECT], etai[DIMPOINT], gammai[DIMPOINT];
    if(i<nx)  // we will compute gammai only if i is in the range
    {
        // load xi, etai, betai from device global memory
        for(int k=0; k<DIMPOINT; k++)
            xi[k] = x[i*DIMPOINT+k];
        for(int k=0; k<DIMVECT; k++)
            betai[k] = beta[i*DIMVECT+k];
        for(int k=0; k<DIMPOINT; k++)
            etai[k] = eta[i*DIMPOINT+k];
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] = 0.0f;
    }

    for(int jstart = 0, tile = 0; jstart < nx; jstart += blockDim.x, tile++)
    {
        int j = tile * blockDim.x + threadIdx.x;
        if(j<nx) // we load xj, etaj and betaj from device global memory only if j<nx
        {
            int inc = 2 * DIMPOINT + DIMVECT;
            for(int k=0; k<DIMPOINT; k++)
                SharedData[threadIdx.x*inc+k] = x[j*DIMPOINT+k];
            for(int k=0; k<DIMVECT; k++)
                SharedData[threadIdx.x*inc+DIMPOINT+k] = beta[j*DIMVECT+k];
            for(int k=0; k<DIMPOINT; k++)
                SharedData[threadIdx.x*inc+DIMPOINT+DIMVECT+k] = eta[j*DIMPOINT+k];
        }
        __syncthreads();
        if(i<nx) // we compute gammai only if i is in the range
        {
            TYPE *xj, *betaj, *etaj;
            xj = SharedData;
            betaj = SharedData + DIMPOINT;
            etaj = SharedData + DIMPOINT + DIMVECT;
            int inc = 2 * DIMPOINT + DIMVECT;
            for(int jrel = 0; jrel < blockDim.x && jrel<nx-jstart; jrel++, xj+=inc, betaj+=inc, etaj+=inc)
                Ker.GradDiff(gammai, xi, xj, betai, betaj, etai, etaj);
        }
        __syncthreads();
    }

    // Save the result in global memory.
    if(i<nx)
        for(int k=0; k<DIMPOINT; k++)
            gamma[i*DIMPOINT+k] = gammai[k];
}



template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGradDiffConv1D(KER Ker,
        TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h,
         int nx)
{

    // Data on the device.
    TYPE* x_d;
    TYPE* beta_d;
    TYPE* gamma_d;
    TYPE* eta_d;

    // Allocate arrays on device.
    cudaMalloc((void**)&x_d, sizeof(TYPE)*(nx*DIMPOINT));
    cudaMalloc((void**)&beta_d, sizeof(TYPE)*(nx*DIMVECT));
    cudaMalloc((void**)&eta_d, sizeof(TYPE)*(nx*DIMPOINT));
    cudaMalloc((void**)&gamma_d, sizeof(TYPE)*(nx*DIMPOINT));

    // Send data from host to device.
    cudaMemcpy(x_d, x_h, sizeof(TYPE)*(nx*DIMPOINT), cudaMemcpyHostToDevice);
    cudaMemcpy(beta_d, beta_h, sizeof(TYPE)*(nx*DIMVECT), cudaMemcpyHostToDevice);
    cudaMemcpy(eta_d, eta_h, sizeof(TYPE)*(nx*DIMPOINT), cudaMemcpyHostToDevice);

    // compute on device.
    dim3 blockSize;
    blockSize.x = 192; // number of threads in each block
    dim3 gridSize;
    gridSize.x =  nx / blockSize.x + (nx%blockSize.x==0 ? 0 : 1);

    GpuGradDiffConv1DOnDevice<TYPE,DIMPOINT,DIMVECT,KER>
        <<<gridSize,blockSize,blockSize.x*(2*DIMPOINT+DIMVECT)*sizeof(TYPE)>>>
            (Ker, x_d, beta_d, eta_d, gamma_d, nx);

    // block until the device has completed
    cudaThreadSynchronize();

    // Send data from device to host.
    cudaMemcpy(gamma_h, gamma_d, sizeof(TYPE)*(nx*DIMPOINT),cudaMemcpyDeviceToHost);

    // Free memory.
    cudaFree(x_d);
    cudaFree(eta_d);
    cudaFree(beta_d);
    cudaFree(gamma_d);

    return 0;
}



template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradDiffConv1D(TYPE sigma, TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h, int nx)
{
	return GpuGradDiffConv1D < TYPE, DIMPOINT, DIMVECT, ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> > >
		(ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> >(GaussFunction<TYPE>(sigma)),
			x_h, beta_h, eta_h, gamma_h, nx);
}





////////////////////////////////////////////////////////////////////////////////////////////////////
// Diff Conv1D :
////////////////////////////////////////////////////////////////////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
__global__ void GpuDiffConv1DOnDevice(KER Ker,
        TYPE *x, TYPE *beta, TYPE *eta, TYPE *gamma,
        int nx)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ char SharedData_char[];
    TYPE* const SharedData = reinterpret_cast<TYPE*>(SharedData_char);

    TYPE xi[DIMPOINT], etai[DIMPOINT], gammai[DIMVECT];
    if(i<nx)  // we will compute gammai only if i is in the range
    {
        // load xi, etai from device global memory
        for(int k=0; k<DIMPOINT; k++)
            xi[k] = x[i*DIMPOINT+k];
        for(int k=0; k<DIMPOINT; k++)
            etai[k] = eta[i*DIMPOINT+k];
        for(int k=0; k<DIMVECT; k++)
            gammai[k] = 0.0f;
    }

    for(int jstart = 0, tile = 0; jstart < nx; jstart += blockDim.x, tile++)
    {
        int j = tile * blockDim.x + threadIdx.x;
        if(j<nx) // we load xj, betaj and etaj from device global memory only if j<nx
        {
            int inc = 2 * DIMPOINT + DIMVECT;
            for(int k=0; k<DIMPOINT; k++)
                SharedData[threadIdx.x*inc+k] = x[j*DIMPOINT+k];
            for(int k=0; k<DIMVECT; k++)
                SharedData[threadIdx.x*inc+DIMPOINT+k] = beta[j*DIMVECT+k];
            for(int k=0; k<DIMPOINT; k++)
                SharedData[threadIdx.x*inc+DIMPOINT+DIMVECT+k] = eta[j*DIMPOINT+k];
        }
        __syncthreads();
        if(i<nx) // we compute gammai only if i is in the range
        {
            TYPE *xj, *betaj, *etaj;
            xj = SharedData;
            betaj = SharedData + DIMPOINT;
            etaj = SharedData + DIMPOINT + DIMVECT;
            int inc = 2 * DIMPOINT + DIMVECT;
            for(int jrel = 0; jrel < blockDim.x && jrel<nx-jstart; jrel++, xj+=inc, betaj+=inc, etaj+=inc)
                Ker.Diff(gammai, xi, xj, betaj, etai, etaj);
        }
        __syncthreads();
    }

    // Save the result in global memory.
    if(i<nx)
        for(int k=0; k<DIMVECT; k++)
            gamma[i*DIMVECT+k] = gammai[k];
}



template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuDiffConv1D(KER Ker,
        TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h,
        int nx)
{

    // Data on the device.
    TYPE* x_d;
    TYPE* beta_d;
    TYPE* gamma_d;
    TYPE* eta_d;

    // Allocate arrays on device.
    cudaMalloc((void**)&x_d, sizeof(TYPE)*(nx*DIMPOINT));
    cudaMalloc((void**)&beta_d, sizeof(TYPE)*(nx*DIMVECT));
    cudaMalloc((void**)&eta_d, sizeof(TYPE)*(nx*DIMPOINT));
    cudaMalloc((void**)&gamma_d, sizeof(TYPE)*(nx*DIMVECT));

    // Send data from host to device.
    cudaMemcpy(x_d, x_h, sizeof(TYPE)*(nx*DIMPOINT), cudaMemcpyHostToDevice);
    cudaMemcpy(beta_d, beta_h, sizeof(TYPE)*(nx*DIMVECT), cudaMemcpyHostToDevice);
    cudaMemcpy(eta_d, eta_h, sizeof(TYPE)*(nx*DIMPOINT), cudaMemcpyHostToDevice);

    // compute on device.
    dim3 blockSize;
    blockSize.x = 192; // number of threads in each block
    dim3 gridSize;
    gridSize.x =  nx / blockSize.x + (nx%blockSize.x==0 ? 0 : 1);

    GpuDiffConv1DOnDevice<TYPE,DIMPOINT,DIMVECT,KER>
        <<<gridSize,blockSize,blockSize.x*(2*DIMPOINT+DIMVECT)*sizeof(TYPE)>>>
            (Ker, x_d, beta_d, eta_d, gamma_d, nx);

    // block until the device has completed
    cudaThreadSynchronize();

    // Send data from device to host.
    cudaMemcpy(gamma_h, gamma_d, sizeof(TYPE)*(nx*DIMVECT),cudaMemcpyDeviceToHost);

    // Free memory.
    cudaFree(x_d);
    cudaFree(eta_d);
    cudaFree(beta_d);
    cudaFree(gamma_d);

    return 0;
}



template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuDiffConv1D(TYPE sigma, TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h, int nx)
{
	return GpuDiffConv1D < TYPE, DIMPOINT, DIMVECT, ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> > >
		(ScalarRadialKernel<TYPE,DIMPOINT,DIMVECT,GaussFunction<TYPE> >(GaussFunction<TYPE>(sigma)), 
			x_h, beta_h, eta_h, gamma_h, nx);
}





// http://www.parashift.com/c++-faq-lite/separate-template-fn-defn-from-decl.html
#define DECLARE_Conv1DS(TYPE,DIMPOINT,DIMVECT) \
	template int GaussGpuEvalConv1D<TYPE,DIMPOINT,DIMVECT>(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int, int); \
	template int GaussGpuGrad1Conv1D<TYPE,DIMPOINT,DIMVECT>(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int); \
	template int GaussGpuGradDiffConv1D<TYPE,DIMPOINT,DIMVECT>(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int); \
	template int GaussGpuDiffConv1D<TYPE,DIMPOINT,DIMVECT>(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);
#define DECLARE_Conv1DS_ALLDIMS_FOR(TYPE) \
	DECLARE_Conv1DS(TYPE,1,1) \
	DECLARE_Conv1DS(TYPE,2,1) \
	DECLARE_Conv1DS(TYPE,2,2) \
	DECLARE_Conv1DS(TYPE,2,4) \
	DECLARE_Conv1DS(TYPE,3,1) \
	DECLARE_Conv1DS(TYPE,3,3) \
	DECLARE_Conv1DS(TYPE,3,6)
DECLARE_Conv1DS_ALLDIMS_FOR(float)
DECLARE_Conv1DS_ALLDIMS_FOR(double)




#endif /* _GpuConv1D_cu */
