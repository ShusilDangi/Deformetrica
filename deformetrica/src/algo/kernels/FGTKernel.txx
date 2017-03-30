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

#ifndef _FGTKernel_txx
#define _FGTKernel_txx

#include "FGTKernel.h"

#include <exception>
#include <stdexcept>

#include <iostream>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
FGTKernel<TScalar, PointDim>
::FGTKernel() : Superclass()
 {
	this->Init();
	m_FGTObject = 0;
 }



template <class TScalar, unsigned int PointDim>
FGTKernel<TScalar, PointDim>
::FGTKernel(const FGTKernel& o)
 {
	Superclass::m_Sources = o.m_Sources;
	Superclass::m_Weights = o.m_Weights;
	this->SetKernelWidth(o.GetKernelWidth());

	if (o.IsModified())
		this->SetModified();
	else
		this->UnsetModified();

	m_FarRatio = o.m_FarRatio;
	m_NumberOfClusters = o.m_NumberOfClusters;
	m_TruncationOrder = o.m_TruncationOrder;

	if (o.m_FGTObject == 0)
	{
		m_FGTObject = 0;
	}
	else
	{
		delete m_FGTObject;
		m_FGTObject = o.m_FGTObject->Clone();
	}
 }



template <class TScalar, unsigned int PointDim>
FGTKernel<TScalar, PointDim>
::FGTKernel(const VNLMatrixType& X, double h): Superclass(X, h)
 {
	this->Init();
	m_FGTObject = this->BuildFGT(Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
 }



template <class TScalar, unsigned int PointDim>
FGTKernel<TScalar, PointDim>
::FGTKernel(
		const VNLMatrixType& X, const VNLMatrixType& W,
		double h): Superclass(X, W, h)
		{
	this->Init();
	m_FGTObject = this->BuildFGT(Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
		}



template <class TScalar, unsigned int PointDim>
FGTKernel<TScalar, PointDim>*
FGTKernel<TScalar, PointDim>
::Clone() const
 {
	return new FGTKernel(*this);
 }



template <class TScalar, unsigned int PointDim>
FGTKernel<TScalar, PointDim>
::~FGTKernel()
{
	delete m_FGTObject;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
typename FGTKernel<TScalar, PointDim>::VNLMatrixType
FGTKernel<TScalar, PointDim>
::Convolve(const VNLMatrixType& X)
 {
	if (X.columns() != PointDim)
		throw std::runtime_error("Invalid location dimension");

	if (this->IsModified())
	{
		delete m_FGTObject;
		m_FGTObject = this->BuildFGT(
				Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
		this->UnsetModified();
	}

	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	unsigned int comboDim =
			weightDim + weightDim*PointDim + weightDim*PointDim*(PointDim+1)/2;

	VNLMatrixType fgtout(X.rows(), comboDim, 0.0);
	m_FGTObject->Compute_v_i(
			X.data_block(), X.rows(), fgtout.data_block());

	VNLMatrixType V(X.rows(), weightDim, 0.0);
	for (unsigned int i = 0; i < X.rows(); i++)
	{
		for (unsigned int j = 0; j < weightDim; j++)
			V(i, j) = fgtout(i, j);
	}

	return V;
 }



template <class TScalar, unsigned int PointDim>
std::vector< typename FGTKernel<TScalar, PointDim>::VNLMatrixType >
FGTKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X)
 {
	if (X.columns() != PointDim)
		throw std::runtime_error("Invalid location dimension");

	if (this->IsModified())
	{
		delete m_FGTObject;
		m_FGTObject = this->BuildFGT(
				Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
		this->UnsetModified();
	}

	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	unsigned int comboDim =
			weightDim + weightDim*PointDim + weightDim*PointDim*(PointDim+1)/2;

	VNLMatrixType fgtout(X.rows(), comboDim, 0.0);
	m_FGTObject->Compute_v_i(
			X.data_block(), X.rows(), fgtout.data_block());

	std::vector<VNLMatrixType> gradK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLVectorType xi = X.get_row(i);
		VNLVectorType fi = fgtout.get_row(i);

		VNLMatrixType G(weightDim, PointDim, 0.0);

		for (unsigned int r = 0; r < weightDim; r++)
			for (unsigned int c = 0; c < PointDim; c++)
				G(r, c) = fi[r]*xi[c] - fi[weightDim+PointDim*r+c];

		G *= (-2.0 / Superclass::m_KernelWidthSquared);

		gradK.push_back(G);
	}

	return gradK;
 }



template <class TScalar, unsigned int PointDim>
typename FGTKernel<TScalar, PointDim>::VNLVectorType
FGTKernel<TScalar, PointDim>
::ConvolveGradient(const VNLMatrixType& X, unsigned int w, unsigned int dp)
 {
	if (X.columns() != PointDim)
		throw std::runtime_error("Invalid location dimension");

	if (this->IsModified())
	{
		delete m_FGTObject;
		m_FGTObject = this->BuildFGT(
				Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
		this->UnsetModified();
	}

	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	unsigned int comboDim =
			weightDim + weightDim*PointDim + weightDim*PointDim*(PointDim+1)/2;

	VNLMatrixType fgtout(X.rows(), comboDim, 0.0);
	m_FGTObject->Compute_v_i(
			X.data_block(), X.rows(), fgtout.data_block());

	VNLVectorType gradK(X.rows(), 0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLVectorType xi = X.get_row(i);
		VNLVectorType fi = fgtout.get_row(i);

		VNLMatrixType G(weightDim, PointDim, 0.0);

		gradK[i] = fi[w]*xi[dp] - fi[weightDim+PointDim*w+dp];
	}

	gradK *= (-2.0 / Superclass::m_KernelWidthSquared);

	return gradK;
 }



template <class TScalar, unsigned int PointDim>
std::vector< std::vector< typename FGTKernel<TScalar, PointDim>::VNLMatrixType > >
FGTKernel<TScalar, PointDim>
::ConvolveHessian(const VNLMatrixType& X)
 {
	if (X.columns() != PointDim)
		throw std::runtime_error("Invalid location dimension");

	if (this->IsModified())
	{
		delete m_FGTObject;
		m_FGTObject = this->BuildFGT(
				Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
		this->UnsetModified();
	}

	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	unsigned int comboDim =
			weightDim + weightDim*PointDim + weightDim*PointDim*(PointDim+1)/2;

	VNLMatrixType fgtout(X.rows(), comboDim, 0.0);
	m_FGTObject->Compute_v_i(
			X.data_block(), X.rows(), fgtout.data_block());

	std::vector< std::vector<VNLMatrixType> > hessK;

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLVectorType xi = X.get_row(i);
		VNLVectorType fi = fgtout.get_row(i);

		std::vector<VNLMatrixType> H;

		for (unsigned int j = 0; j < weightDim; j++)
		{
			double wconv = fi[j];

			VNLVectorType wyconv(PointDim, 0.0);
			for (unsigned int c = 0; c < PointDim; c++)
				wyconv[c] = fi[weightDim+PointDim*j+c];

			unsigned int wysq_ind =
					weightDim + weightDim*PointDim + j*PointDim*(PointDim+1)/2;

			VNLMatrixType wysqconv(PointDim, PointDim, 0.0);
			for (unsigned int r = 0; r < PointDim; r++)
				for (unsigned int c = r; c < PointDim; c++)
				{
					wysqconv(r, c) = fi[wysq_ind++];
					wysqconv(c, r) = wysqconv(r, c);
				}

			VNLMatrixType& Hj = wysqconv;
			for (unsigned int r = 0; r < PointDim; r++)
				for (unsigned int c = 0; c < PointDim; c++)
				{
					Hj(r, c) += X(i, r)*X(i, c)*wconv
							- wyconv[r]*X(i,c)
							- X(i,r)*wyconv[c];
				}

			Hj *= (4.0 / (Superclass::m_KernelWidthSquared*Superclass::m_KernelWidthSquared));

			for (unsigned int r = 0; r < PointDim; r++)
				Hj(r, r) -= 2.0 / Superclass::m_KernelWidthSquared * wconv;

			H.push_back(Hj);
		}

		hessK.push_back(H);
	}

	return hessK;
 }



template <class TScalar, unsigned int PointDim>
typename FGTKernel<TScalar, PointDim>::VNLVectorType
FGTKernel<TScalar, PointDim>
::ConvolveHessian(const VNLMatrixType& X, unsigned int w,
		unsigned int dp, unsigned int dq)
		{
	if (X.columns() != PointDim)
		throw std::runtime_error("Invalid location dimension");

	if (this->IsModified())
	{
		delete m_FGTObject;
		m_FGTObject = this->BuildFGT(
				Superclass::m_Sources, Superclass::m_Weights, Superclass::m_KernelWidth);
		this->UnsetModified();
	}

	VNLMatrixType& Y = this->GetSources();
	VNLMatrixType& W = this->GetWeights();

	if (Y.rows() != W.rows())
		throw std::runtime_error("Sources and weights count mismatch");

	unsigned int weightDim = W.columns();

	unsigned int comboDim =
			weightDim + weightDim*PointDim + weightDim*PointDim*(PointDim+1)/2;

	VNLMatrixType fgtout(X.rows(), comboDim, 0.0);
	m_FGTObject->Compute_v_i(
			X.data_block(), X.rows(), fgtout.data_block());

	double width4 = Superclass::m_KernelWidthSquared * Superclass::m_KernelWidthSquared;

	VNLVectorType hessK(X.rows(), 0);

	for (unsigned int i = 0; i < X.rows(); i++)
	{
		VNLVectorType xi = X.get_row(i);
		VNLVectorType fi = fgtout.get_row(i);

		double wconv = fi[w];

		double  wyconv_p = fi[weightDim+PointDim*w+dp];
		double  wyconv_q = fi[weightDim+PointDim*w+dq];

		unsigned int wysq_ind =
				weightDim + weightDim*PointDim + w*PointDim*(PointDim+1)/2;

		hessK[i] = fi[wysq_ind + PointDim*dp + dq];

		hessK[i] += X(i,dp)*X(i,dq)*wconv - wyconv_p*X(i,dq) - X(i,dp)*wyconv_q;
		hessK[i] *= (4.0 / width4);

		if (dp == dq)
			hessK[i] -= 2.0 / Superclass::m_KernelWidthSquared * wconv;
	}

	return hessK;
		}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int PointDim>
void
FGTKernel<TScalar, PointDim>
::Init()
 {
	m_FarRatio = 5.0;
	m_TruncationOrder = 8;
	m_NumberOfClusters = 50;
 }



template <class TScalar, unsigned int PointDim>
typename FGTKernel<TScalar, PointDim>::FGCalculatorType*
FGTKernel<TScalar, PointDim>
::BuildFGT(
		const VNLMatrixType& S, const VNLMatrixType& W, double h)
{
	unsigned int numPoints = S.rows();

	unsigned int weightDim = W.columns();

	unsigned int comboDim =
			weightDim + weightDim*PointDim + weightDim*PointDim*(PointDim+1)/2;

	VNLMatrixType C(S.rows(), comboDim, 0.0);

	for (unsigned int i = 0; i < numPoints; i++)
	{
		for (unsigned int j = 0; j < weightDim; j++)
			C(i, j) = W(i, j);

		unsigned int col = weightDim;

		// For gradient term: sum{w * y}
		for (unsigned int r = 0; r < weightDim; r++)
			for (unsigned int c = 0; c < PointDim; c++)
			{
				C(i, col++) = W(i, r) * S(i, c);
			}

		// For Hessian term: sum{w * y^2}
		for (unsigned int p = 0; p < weightDim; p++)
			for (unsigned int r = 0; r < PointDim; r++)
				for (unsigned int c = r; c < PointDim; c++)
				{
					C(i, col++) = W(i, p) * S(i, r) * S(i, c);
				}
	}

	unsigned int numClust = m_NumberOfClusters;
	if (numPoints < numClust)
		numClust = numPoints / 2 + 1;

	TScalar ptRadius = 0;

	FGCalculatorType* fgtobj = new FGCalculatorType(
			PointDim, comboDim,
			C.data_block(), S.data_block(),
			S.rows(), Superclass::m_KernelWidth,
			m_TruncationOrder, numClust, m_FarRatio,
			&ptRadius);

	return fgtobj;
}



#endif /* _FGTKernel_txx */
