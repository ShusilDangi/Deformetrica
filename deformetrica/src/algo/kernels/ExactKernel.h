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

#ifndef _ExactKernel_h
#define _ExactKernel_h

#include <cmath>
#include <vector>

#include "AbstractKernel.h"

/**
 *	\brief          Exact kernels.
 *
 *	\copyright	Inria and the University of Utah
 *	\version        Deformetrica 2.1
 *
 *	\details        The ExactKernel class inherited from AbstractKernel implements the operations of
 *                      convolution and evaluation of the kernel using an exact computation (i.e. no approximation).
 */
template <class TScalar, unsigned int PointDim>
class ExactKernel : public AbstractKernel<TScalar, PointDim>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Abstract kernel type.
	typedef AbstractKernel<TScalar, PointDim> Superclass;

	/// Vector type
	typedef typename Superclass::VNLVectorType VNLVectorType;
	/// Matrix type.
	typedef typename Superclass::VNLMatrixType VNLMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	ExactKernel(): Superclass() { }
	/// Copy constructor.
	ExactKernel(const ExactKernel& o);
	/// See AbstractKernel::AbstractKernel(const VNLMatrixType& X, double h).
	ExactKernel(const VNLMatrixType& X, double h): Superclass(X, h) {}
	/// See AbstractKernel::AbstractKernel(const VNLMatrixType& X, const VNLMatrixType& W, double h).
	ExactKernel(const VNLMatrixType& X, const VNLMatrixType& W, double h): Superclass(X, W, h) {}

	virtual ExactKernel* Clone() const;

	virtual ~ExactKernel() { }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	inline TScalar EvaluateKernel(const VNLVectorType& x, const VNLVectorType& y)
	{
		// VNLVectorType dvec = x - y;
		TScalar distsq = (x-y).squared_magnitude();
		return exp(-distsq / Superclass::m_KernelWidthSquared);
	}

	/// Evaluates \f$ K(x_{row_x},y_{row_y}) \f$.
	inline TScalar EvaluateKernel(const VNLMatrixType& X, const VNLMatrixType& Y, size_t row_x, size_t row_y)
	{
		TScalar dist_squared = 0.0;
		TScalar diff;
		for(size_t j = 0; j < X.cols(); j++)
		{
			diff = X(row_x, j) - Y(row_y, j);
			dist_squared += diff*diff;
		}
		return exp(-dist_squared / Superclass::m_KernelWidthSquared);
	}

	// Gradient of kernel(x-y) at x
	virtual VNLVectorType EvaluateKernelGradient(const VNLVectorType& x, const VNLVectorType& y)
	{
		// VNLVectorType dvec = x - y;
		TScalar distsq = (x-y).squared_magnitude();
		TScalar k = exp(-distsq / Superclass::m_KernelWidthSquared);
		// dvec *= -2.0 * k / Superclass::m_KernelWidthSquared;
		// return dvec;
		return ( (x-y) * (-2.0 * k / Superclass::m_KernelWidthSquared) );
	}

	/// Evaluates the gradient of \f$ K(x_{row_x},y_{row_y}) \f$ at \e \f$ x_{row_x} \f$.
	inline VNLVectorType EvaluateKernelGradient(const VNLMatrixType& X, const VNLMatrixType& Y, size_t row_x, size_t row_y)
	{
		VNLVectorType result(PointDim);
		TScalar dist_squared = 0.0;
		TScalar diff;
		for(size_t j = 0; j < X.cols(); j++)
		{
			diff = X(row_x, j) - Y(row_y, j);
			result(j) = diff;
			dist_squared += diff*diff;
		}

		return (result * (-2.0 * exp(-dist_squared / Superclass::m_KernelWidthSquared) / Superclass::m_KernelWidthSquared ));
	}

	// Hessian of kernel(x-y) at x
	virtual VNLMatrixType EvaluateKernelHessian(const VNLVectorType& x, const VNLVectorType& y)
	{
		int dims = x.size();

		VNLMatrixType xy(dims, 1, 0.0);
		TScalar distsq = 0;
		for (int d = 0; d < dims; d++)
		{
			TScalar t = x[d] - y[d];
			xy(d, 0) = t;
			distsq += t*t;
		}

		VNLMatrixType H = xy * xy.transpose();

		TScalar k = exp(-distsq / Superclass::m_KernelWidthSquared);
		H *= 4.0 * k / (Superclass::m_KernelWidthSquared * Superclass::m_KernelWidthSquared);

		for (int i = 0; i < dims; i++)
			H(i, i) -= 2.0 * k / Superclass::m_KernelWidthSquared;

		return H;
	}

	/// Evaluates the hessian of \f$ K(x_{row_x},y_{row_y}) \f$ at \e \f$ x_{row_x} \f$.
	virtual VNLMatrixType EvaluateKernelHessian(const VNLMatrixType& X, const VNLMatrixType& Y, size_t row_x, size_t row_y)
	{
		VNLMatrixType H(PointDim, PointDim, 0.0);
		VNLVectorType x_minus_y(PointDim, 0.0);

		TScalar dist_squared = 0.0;
		TScalar diff;
		for(size_t j = 0; j < X.cols(); j++)
		{
			diff = X(row_x, j) - Y(row_y, j);
			x_minus_y(j) = diff;
			dist_squared += diff*diff;
		}

		for (int j = 0; j < PointDim; j++)
			for (int i = 0; i < PointDim; i++)
				H(i, j) = 4.0 * x_minus_y(i) * x_minus_y(j) / (Superclass::m_KernelWidthSquared * Superclass::m_KernelWidthSquared);

		for (int i = 0; i < PointDim; i++)
			H(i, i) -= 2.0 / Superclass::m_KernelWidthSquared;

		TScalar k_ij = exp(-dist_squared / Superclass::m_KernelWidthSquared);
		return k_ij * H; //(-2.0 * exp(-dist_squared / Superclass::m_KernelWidthSquared) / Superclass::m_KernelWidthSquared ) * result;
	}



	virtual VNLMatrixType Convolve(const VNLMatrixType& X);

	virtual std::vector<VNLMatrixType> ConvolveGradient(const VNLMatrixType& X);
	virtual VNLMatrixType ConvolveGradient(const VNLMatrixType& X, const VNLMatrixType& alpha);

	// Derivative of convolved weight w at direction dp
	virtual VNLVectorType ConvolveGradient(const VNLMatrixType& X, unsigned int k, unsigned int dp);
	virtual VNLMatrixType ConvolveGradient(const VNLMatrixType& X, unsigned int dim);

	virtual std::vector< std::vector<VNLMatrixType> > ConvolveHessian(const VNLMatrixType & X);

	// Second derivative of weight w at directions dp and dq
	virtual VNLVectorType ConvolveHessian(const VNLMatrixType & X, unsigned int k, unsigned int dp, unsigned int dq);
	virtual VNLMatrixType ConvolveHessian(const VNLMatrixType& X, unsigned int row, unsigned int col);



protected:


}; /* class ExactKernel */


#ifndef MU_MANUAL_INSTANTIATION
#include "ExactKernel.txx"
#endif


#endif /* _ExactKernel_h */

