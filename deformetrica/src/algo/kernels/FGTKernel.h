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

#ifndef _FGTKernel_h
#define _FGTKernel_h

#include "ExactKernel.h"
#include "FastGauss.h"

/**
 *	\brief 		Kernels with Fast Gauss Transform approximation.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The ExactKernel class inherited from AbstractKernel implements the operations with Gaussian kernels using the Fast Gauss Transforms (i.e. multipole approximation).
 */
template <class TScalar, unsigned int PointDim>
class FGTKernel: public ExactKernel<TScalar, PointDim>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> Superclass;

	/// Vector type.
	typedef typename Superclass::VNLVectorType VNLVectorType;
	/// Matrix type.
	typedef typename Superclass::VNLMatrixType VNLMatrixType;

	/// Fast Gauss type.
	typedef FastGauss<TScalar> FGCalculatorType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	FGTKernel();
	/// Copy constructor.
	FGTKernel(const FGTKernel& o);
	/// See AbstractKernel::AbstractKernel(const VNLMatrixType& X, double h) for details.
	FGTKernel(const VNLMatrixType& X, double h);
	/// See AbstractKernel::AbstractKernel(const VNLMatrixType& X, const VNLMatrixType& W, double h) for details.
	FGTKernel(const VNLMatrixType& X, const VNLMatrixType& W, double h);

	virtual FGTKernel* Clone() const;

	virtual ~FGTKernel();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the number of clusters to \e n.
	void SetNumberOfClusters(unsigned int n) { m_NumberOfClusters = n; }

	/// FGT parameter.
	void SetFarRatio(double f) { m_FarRatio = f; }

	/// Sets the truncation order to \e n.
	void SetTruncationOrder(unsigned int n) { m_TruncationOrder = n; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual VNLMatrixType Convolve(const VNLMatrixType& X);

	virtual std::vector<VNLMatrixType> ConvolveGradient(const VNLMatrixType& X);
	virtual VNLVectorType ConvolveGradient(const VNLMatrixType& X, unsigned int k, unsigned int dp);

	virtual std::vector< std::vector<VNLMatrixType> > ConvolveHessian(const VNLMatrixType & X);
	virtual VNLVectorType ConvolveHessian(const VNLMatrixType& X, unsigned int k, unsigned int dp, unsigned int dq);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void Init();

	FGCalculatorType* BuildFGT(const VNLMatrixType& sources, const VNLMatrixType& weights, double h);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	FGCalculatorType* m_FGTObject;

	/// Number of clusters in the point set.
	unsigned int m_NumberOfClusters;

	/// Truncation order of the Taylor expansion.
	unsigned int m_TruncationOrder;

	/// FGT parameter.
	double m_FarRatio;


}; /* class FGTKernel */


#ifndef MU_MANUAL_INSTANTIATION
#include "FGTKernel.txx"
#endif


#endif /* _FGTKernel_h */
