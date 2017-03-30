/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _PointCloud_h
#define _PointCloud_h

#include "DeformableObject.h"
#include "Landmark.h"

#include "Diffeos.h"

#include "vnl/vnl_matrix.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <vector>

#include <cstring>
#include <iostream>
#include <sstream>


/**
 *	\brief 		Point clouds.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The PointCloud class inherited from Landmark represents an unstructured unlabelled point set.
 *				This class uses the 0-current representation (aka measures) which does not assume point-to-point correspondence between source and target point clouds to be match.
 */
template <class TScalar, unsigned int Dimension>
class PointCloud : public Landmark<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Landmark type.
	typedef Landmark<TScalar, Dimension> Superclass;

	/// Deformable object type.
	typedef typename Superclass::Superclass DeformableObjectType;

	/// Vector type
	typedef typename Superclass::VNLVectorType VNLVectorType;
	/// Matrix type.
	typedef typename Superclass::VNLMatrixType VNLMatrixType;
	/// List of matrices type.
	typedef typename Superclass::VNLMatrixList VNLMatrixList;

	/// Deformation type.
	typedef typename Superclass::DiffeosType DiffeosType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	PointCloud();
	/// Copy constructor.
	PointCloud(const PointCloud& other);

	virtual PointCloud* Clone() { return new PointCloud(*this); }

	virtual ~PointCloud();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void SetPolyData(vtkPolyData* polyData);

	/// Returns the weights associated to the points.
	inline VNLMatrixType GetPointWeights() const { return m_PointWeights; }

	///	Returns the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e kernelWidth.
	void SetKernelWidth(TScalar h);

	/// Returns the RKHS-norm of itself.
	inline TScalar GetNormSquared() const { return m_NormSquared; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// See DeformableObject::ComputeMatch(DeformableObject* target) for details.
	virtual TScalar ComputeMatch(DeformableObjectType* target);

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual VNLMatrixType ComputeMatchGradient(DeformableObjectType* target);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the weights associated to the points.
	void UpdatePointWeights();

	/// Computes the RKHS-norm of itself.
	TScalar ComputeSelfNorm();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Weights associated to the points (Size : NumberOfPoints x 1).
	VNLMatrixType m_PointWeights;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the point cloud.
	TScalar m_NormSquared;


}; /* class PointCloud */


#ifndef MU_MANUAL_INSTANTIATION
#include "PointCloud.txx"
#endif


#endif /* _PointCloud_h */
