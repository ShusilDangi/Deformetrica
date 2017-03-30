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

#ifndef _ExactKernel_txx

#include "ExactKernel.h"

#include <exception>
#include <stdexcept>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
ExactKernel<TScalar, PointDim>
::ExactKernel(const ExactKernel& o)
 {
	Superclass::m_Sources = o.m_Sources;
	Superclass::m_Weights = o.m_Weights;
	this->SetKernelWidth(o.GetKernelWidth());

	if (o.IsModified())
		this->SetModified();
	else
		this->UnsetModified();
 }



template <class TScalar, unsigned int PointDim>
ExactKernel<TScalar, PointDim>*
ExactKernel<TScalar, PointDim>
::Clone() const
 {
	return new ExactKernel(*this);
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VNLMatrixType
ExactKernel<TScalar, PointDim>
::Convolve(const VNLMatrixType& X)
 {
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	VNLMatrixType V(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			TScalar Kij = this->EvaluateKernel(X, Y, i, j);

			for(unsigned int k=0; k< weightDim; k++)
				V(i,k) += W(j,k) * Kij;
		}
	}

	return V;
 }



template <class TScalar, unsigned int PointDim>
std::vector<typename ExactKernel<TScalar, PointDim>::VNLMatrixType>
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X)
 {
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	std::vector<VNLMatrixType> gradK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLMatrixType Gi(weightDim, PointDim, 0.0);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VNLVectorType g = this->EvaluateKernelGradient(X, Y, i, j);
			for (unsigned int k = 0; k < weightDim; k++)
			{
				TScalar Wjk = W(j,k);
				for (unsigned int l = 0; l < PointDim; l++)
				{
					Gi(k,l) +=  g[l]*Wjk;
				}
			}
		}

		gradK.push_back(Gi);
	}

	return gradK;
}



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VNLMatrixType
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X, const VNLMatrixType& alpha)
 {
	std::vector<VNLMatrixType> convolveGradient = this->ConvolveGradient(X);

	VNLMatrixType result(X.rows(), PointDim, 0);
	for (unsigned int j = 0; j < X.rows(); j++)
		result.set_row(j, convolveGradient[j].transpose() * alpha.get_row(j) );

	return result;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VNLMatrixType
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X, unsigned int dim)
 {
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");
	if (dim >= PointDim || dim < 0)
		throw std::runtime_error("dimension index out of bounds");

	unsigned int weightDim = W.columns();

	VNLMatrixType gradK(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLVectorType xi = X.get_row(i);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VNLVectorType yj = Y.get_row(j);
			VNLVectorType wj = W.get_row(j);

			VNLVectorType g = this->EvaluateKernelGradient(xi, yj);
			gradK.set_row(i, gradK.get_row(i) + g[dim]*wj);
		}
	}

	return gradK;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VNLVectorType
ExactKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X, unsigned int k, unsigned int dp)
 {
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	if (k >= weightDim)
		throw std::runtime_error("Invalid weight index");

	if (dp >= PointDim)
		throw std::runtime_error("Invalid derivative direction");

	unsigned int numPoints = Y.rows();

	VNLVectorType gradK(numPoints, 0);

	for (unsigned int i = 0; i < numPoints; i++)
	{
		VNLVectorType xi = X.get_row(i);

		gradK[i] = 0;

		for (unsigned int j = 0; j < numPoints; j++)
		{
			VNLVectorType yj = Y.get_row(j);
			VNLVectorType wj = W.get_row(j);

			VNLVectorType g = this->EvaluateKernelGradient(xi, yj);

			gradK[i] += wj[k] * g[dp];
		}
	}

	return gradK;
 }



template <class TScalar, unsigned int PointDim>
std::vector< std::vector<typename ExactKernel<TScalar, PointDim>::VNLMatrixType> >
ExactKernel<TScalar, PointDim>
::ConvolveHessian(const VNLMatrixType& X)
 {
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	std::vector< std::vector<VNLMatrixType> > hessK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		std::vector<VNLMatrixType> Hi;
		for (unsigned int k = 0; k < weightDim; k++)
			Hi.push_back(VNLMatrixType(PointDim, PointDim, 0.0));

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VNLVectorType wj = W.get_row(j);

			VNLMatrixType H = this->EvaluateKernelHessian(X, Y, i, j);
			for (unsigned int k = 0; k < weightDim; k++)
				Hi[k] = Hi[k] + H*wj[k];
		}

		hessK.push_back(Hi);
	}

	return hessK;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VNLMatrixType
ExactKernel<TScalar, PointDim>
::ConvolveHessian(const VNLMatrixType& X, unsigned int row, unsigned int col)
 {
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");
	if (row >= PointDim || col >= PointDim)
		throw std::runtime_error("Dimension index out of bounds");

	unsigned int weightDim = W.columns();

	VNLMatrixType hessK(X.rows(), weightDim, 0.0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLVectorType xi = X.get_row(i);

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VNLVectorType yj = Y.get_row(j);
			VNLVectorType wj = W.get_row(j);

			VNLMatrixType H = this->EvaluateKernelHessian(xi, yj);
			hessK.set_row(i, hessK.get_row(i) + wj*H(row,col) );
		}
	}

	return hessK;
 }



template <class TScalar, unsigned int PointDim>
typename ExactKernel<TScalar, PointDim>::VNLVectorType
ExactKernel<TScalar, PointDim>
::ConvolveHessian(const VNLMatrixType& X, unsigned int k,
		unsigned int dp, unsigned int dq)
		{
	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	if (k >= weightDim)
		throw std::runtime_error("Invalid weight index");

	if (dp >= PointDim || dq >= PointDim)
		throw std::runtime_error("Invalid derivative direction");

	unsigned int numPoints = X.rows();

	VNLVectorType hessK(numPoints, 0);

	for (unsigned int i = 0; i < numPoints; i++)
	{
		VNLVectorType xi = X.get_row(i);

		hessK[i] = 0;

		for (unsigned int j = 0; j < Y.rows(); j++)
		{
			VNLVectorType yj = Y.get_row(j);
			VNLVectorType wj = W.get_row(j);

			VNLMatrixType H = this->EvaluateKernelHessian(xi, yj);

			hessK[i] += wj[k] * H(dp, dq);
		}
	}

	return hessK;
}


#endif /* _ExactKernel_txx */

