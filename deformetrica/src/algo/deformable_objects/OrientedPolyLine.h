/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _OrientedPolyLine_h
#define _OrientedPolyLine_h

#include "DeformableObject.h"
#include "Landmark.h"

#include "Diffeos.h"

#include "vnl/vnl_matrix.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <vector>

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *	\brief 		Oriented curves.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The OrientedPolyLine class inherited from Landmark represents a set of polygonal lines.
 *				This class uses the current representation of curves, which is sensitive to the orientation (change in orientation changes the sign of the curve in the space of currents).
 */
template <class TScalar, unsigned int Dimension>
class OrientedPolyLine : public Landmark<TScalar, Dimension>
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

	OrientedPolyLine();
	/// Copy constructor.
	OrientedPolyLine(const OrientedPolyLine& other);

	virtual OrientedPolyLine* Clone() { return new OrientedPolyLine(*this); }

	virtual ~OrientedPolyLine();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void SetPolyData(vtkPolyData* polyData);

	/// Returns the centers of the cells.
	inline VNLMatrixType GetCenters() const { return m_Centers; }

	/// Returns the tangents of the cells.
	inline VNLMatrixType GetTangents() const { return m_Tangents; }

	/// Returns the number of cells.
	inline int GetNumberOfCells() const { return m_NumCells; }

	///	Returns the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e h.
	void SetKernelWidth(TScalar h);

	/// Returns the squared RKHS-norm of itself.
	inline TScalar GetNormSquared() const { return m_NormSquared; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// See DeformableObject::ComputeMatch(DeformableObject* target) for details.
	virtual TScalar ComputeMatch(DeformableObjectType* target);
	/// See DeformableObject::ComputeMatch(DeformableObject* target, unsigned int t) for details.
	virtual TScalar ComputeMatch(DeformableObjectType* target, unsigned int t);

	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target) for details.
	virtual VNLMatrixType ComputeMatchGradient(DeformableObjectType* target);
	/// See DeformableObject::ComputeMatchGradient(DeformableObject* target, unsigned int t) for details.
	virtual VNLMatrixType ComputeMatchGradient(DeformableObjectType* target, unsigned int );



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 *	\brief		Computes the centers and the normals from the points.
	 *
	 *	\details	Given a set of polygonal lines , this method computes from the vertices of the curve
	 *				the centers and the tangents of each cell.
	 *
	 *	\param[in]	Pts			The vertices of the curves.
	 *	\param[out]	Centers		The centers of the cells (Size : NumCells x Dimension).
	 *	\param[out]	Tangents	The tangents of the cells (Size : NumCells x Dimension).
	 */
	void ComputeCentersTangents(const VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Tangents);

	/// Computes the RKHS-norm of itself.
	TScalar ComputeSelfNorm();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
	VNLMatrixType m_Centers;

	///	Matrix coordinates of the tangents of the cells  (Size : NumCells x Dimension).
	VNLMatrixType m_Tangents;

	/// Number of cells (i.e. polygonal lines).
	int m_NumCells;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the oriented curve.
	TScalar m_NormSquared;

	/// See Landmark::m_VTKMutex for details.
	itk::SimpleFastMutexLock m_VTKMutex;


}; /* class OrientedPolyLine */


#ifndef MU_MANUAL_INSTANTIATION
#include "OrientedPolyLine.txx"
#endif


#endif /* _OrientedPolyLine_h */
