/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BackwardGradientPropagatorUniversal_h
#define _BackwardGradientPropagatorUniversal_h

#include "itkImage.h"
#include "itkVector.h"
#include "Diffeos.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "KernelFactory.h"

#include <vector>

/**
 *	\brief 		Computes adjoint equations.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The BackwardGradientPropagatorUniversal class computes the adjoint equations of the geodesic shooting equations and flow equations, which moves the gradient of the data term back to time t = 0. Its value at time t=0 is used to update initial momenta, initial position of control points, and possibly vertices of template shapes.
 */
template <class TScalar, unsigned int Dimension>
class BackwardGradientPropagatorUniversal
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type.
	typedef vnl_vector<TScalar> VNLVectorType;
	/// Matrix type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// List of matrices type.
	typedef std::vector<VNLMatrixType> VNLMatrixList;

	/// Image type (itk).
	typedef itk::Image<TScalar, Dimension> ImageType;
	/// Image pointer type (itk).
	typedef typename ImageType::Pointer ImagePointer;

	/// Vector type (itk).
	typedef itk::Vector<TScalar, Dimension> VectorType;
	/// Vector image type (itk).
	typedef itk::Image<VectorType, Dimension> VectorImageType;
	/// Vector image pointer type (itk).
	typedef typename VectorImageType::Pointer VectorImagePointer;

	/// Kernel factory type.
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	/// Exact kernel type.
	typedef typename KernelFactoryType::KernelBaseType KernelType;

	/// Deformation type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	BackwardGradientPropagatorUniversal();

	~BackwardGradientPropagatorUniversal();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// y(t): grid points y(1) propagated backward
	/// Trajectories of the shape points.
	void SetPointsTrajectory(VNLMatrixList& PtsT) { m_PointsT = PtsT; }

	// eta(1): data match derivative propagated to end time point
	/// Trajectories of points eta(t), where eta(1) is the derivative of the data term, and eta(t) computed previously by integration of the adjoint equation of the flow equation.
	void SetVectorsTrajectory(const VNLMatrixList& VecT) { m_VectorsT = VecT; }

	// Control point positions and momentas
	/// Set the deformation to \e def.
	void SetDiffeos(DiffeosType* def) { m_Def = def; };

	// Xi(0) block 1
	/// Output variable. Auxiliary variable at time t = \q (of size dimension times the number of control points) 
	VNLMatrixType GetGradientPosAt(long q) const { return m_XiPosT[q]; }

	// Xi(0) block 2
	/// Output variable. Auxiliary variable at time t = \q (of size dimension times the number of momentas) .
	VNLMatrixType GetGradientMomAt(long q) const { return m_XiMomT[q]; }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Computes the auxiliary variables at all time points.
	void Update();


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * \brief		Computes the time derivative of the ajoint equation.
	 *
	 * \details		Computes the time derivative of the ajoint equation. It is used during the integration scheme.
	 *
	 * \param[in]	s		Time index.
	 * \param[out]	dPos	output derivative of the auxiliary variable of size dimension times the number of control points.
	 * \param[out]	dMom	output derivative of the auxiliary variable of size dimension times the number of momentas.
	 */
	void ComputeUpdateAt(unsigned int s, VNLMatrixType& dPos,  VNLMatrixType& dMom);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Entity representing the deformation.
	DiffeosType* m_Def;

	// Derivatives in positions and momentas, each Dim x #control points
	// propagated from time 1 to time 0
	VNLMatrixList m_XiPosT;

	VNLMatrixList m_XiMomT;

	VNLMatrixList m_PointsT;

	VNLMatrixList m_VectorsT;

	/// A kernel (three kernels are used for speed purposes).
	KernelType* m_KernelObj1;
	/// Another kernel (three kernels are used for speed purposes).
	KernelType* m_KernelObj2;
	/// Another kernel (three kernels are used for speed purposes).
	KernelType* m_KernelObj3;


}; /* class BackwardGradientPropagatorUniversal */


#ifndef MU_MANUAL_INSTANTIATION
#include "BackwardGradientPropagatorUniversal.txx"
#endif


#endif /* _BackwardGradientPropagatorUniversal_h */
