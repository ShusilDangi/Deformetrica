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

#ifndef _CUDAExactKernel_txx
#define _CUDAExactKernel_txx

#include "CUDAExactKernel.h"

#include "../../lib/cuda_convolutions/GpuConv1D.h"

#include "writeMatrixDLM.txx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
CUDAExactKernel<TScalar, PointDim>
::CUDAExactKernel(const CUDAExactKernel& o)
 {
	Superclass::m_Sources = o.m_Sources;
	Superclass::m_Weights = o.m_Weights;
	this->SetKernelWidth(o.GetKernelWidth());

	if (o.IsModified())
		this->SetModified();
	else
		this->UnsetModified();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
typename CUDAExactKernel<TScalar, PointDim>::VNLMatrixType
CUDAExactKernel<TScalar, PointDim>
::Convolve(const VNLMatrixType& X) {

	if (this->GetSources().rows() != this->GetWeights().rows())
		throw std::runtime_error("Sources and weights count mismatch");

	int DimVect = this->GetWeights().columns();
	VNLMatrixType gamma(X.rows(), DimVect, 0.0);
	TScalar *xp, *yp, *betap, *gammap;
	xp     = (TScalar*)X.data_block();
	yp     = this->GetSources().data_block();
	betap  = this->GetWeights().data_block();
	gammap = gamma.data_block();

	if (DimVect == PointDim)
		GaussGpuEvalConv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				xp, yp, betap, gammap, X.rows(), this->GetSources().rows());
	else if (DimVect == 2*PointDim)
		GaussGpuEvalConv1D<TScalar, PointDim, (2*PointDim)>(this->GetKernelWidth(),
				xp, yp, betap, gammap, X.rows(), this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::Convolve(X) - Invalid number of rows of beta !");

	return gamma;
}



template <class TScalar, unsigned int PointDim>
typename CUDAExactKernel<TScalar, PointDim>::VNLMatrixType
CUDAExactKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X, const VNLMatrixType& alpha)
 {
	int DimVect = this->GetWeights().columns();
	VNLMatrixType gamma(X.rows(), DimVect, 0.0);
	TScalar *xp, *yp, *alphap, *betap, *gammap;
	xp     = (TScalar*)X.data_block();
	yp     = this->GetSources().data_block();
	alphap = (TScalar*)alpha.data_block();
	betap  = this->GetWeights().data_block();
	gammap = gamma.data_block();

	if (DimVect == PointDim)
		GaussGpuGrad1Conv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				alphap, xp, yp, betap, gammap, X.rows(), this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::ConvolveGradient(X, alpha) - Problem with the number of rows of beta !");

	return gamma;
 }



template <class TScalar, unsigned int PointDim>
typename CUDAExactKernel<TScalar, PointDim>::VNLMatrixType
CUDAExactKernel<TScalar, PointDim>
::ConvolveSpecialHessian(const VNLMatrixType& xi)
{
	if (this->GetSources().rows() != this->GetWeights().rows())
		throw std::runtime_error("Sources and weights count mismatch");
	if (this->GetSources().rows() != xi.rows())
		throw std::runtime_error("Y and Xi count mismatch");

	int DimVect = this->GetWeights().columns();
	VNLMatrixType gamma(xi.rows(), DimVect, 0.0);
	TScalar *xi_p, *y_p, *beta_p, *gamma_p;
	xi_p    = (TScalar*)xi.data_block();
	y_p     = this->GetSources().data_block();
	beta_p  = this->GetWeights().data_block();
	gamma_p = gamma.data_block();

	if (DimVect == PointDim)
		GaussGpuGradDiffConv1D<TScalar, PointDim, PointDim>(this->GetKernelWidth(),
				y_p, beta_p, xi_p, gamma_p, this->GetSources().rows());
	else
		throw std::runtime_error("In CUDAExactKernel::ConvolveSpecialHessian(xi) - Invalid number of rows of beta !");

	return gamma * (- 0.5);
}



#endif /* _CUDAExactKernel_txx */
