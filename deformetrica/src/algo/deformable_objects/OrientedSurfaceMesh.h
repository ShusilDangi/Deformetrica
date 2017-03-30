/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _OrientedSurfaceMesh_h
#define _OrientedSurfaceMesh_h

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
 *	\brief 		Oriented surface meshes.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The OrientedSurfaceMesh class inherited from Landmark represents a triangular mesh.
 *				This class uses the current representation of the surface, which is sensitive to the local orientation of the mesh (a flip in ordering of the vertex indices in the connectivity matrix changes the sign of the surface element (i.e. the triangle) in the space of currents).
 */
template <class TScalar, unsigned int Dimension>
class OrientedSurfaceMesh : public Landmark<TScalar, Dimension>
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

	OrientedSurfaceMesh();
	/// Copy constructor.
	OrientedSurfaceMesh(const OrientedSurfaceMesh& other);

	virtual OrientedSurfaceMesh* Clone() { return new OrientedSurfaceMesh(*this); }

	virtual ~OrientedSurfaceMesh();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void SetPolyData(vtkPolyData* polyData) { this->SetPolyData(polyData, false); }
	/// Sets the pointer on a VTK object to \e polyData with possible reorientation.
	void SetPolyData(vtkPolyData* polyData, bool reorient);

	/// Returns the centers of the cells.
	inline VNLMatrixType GetCenters() const { return m_Centers; }

	/// Returns the normals of the cells.
	inline VNLMatrixType GetNormals() const { return m_Normals; }

	/// Returns the number of cells.
	inline int GetNumberOfCells() const { return m_NumCells; }

	///	Returns the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e h.
	void SetKernelWidth(TScalar h);

	/// Returns the RKHS-norm of itself.
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

	/// Possibly reorient normals according to vtk filters.
	/// \warning	Make sure the ordering of point cells is consistent with the direction of the normal.
	void CheckMeshAndNormals(bool reorient);

	/**
	 *	\brief		Computes the centers and the normals from the points.
	 *
	 *	\details	Given a triangular mesh, this method computes from the vertices of the mesh
	 *				the centers and the normals of each face.
	 *
	 *	\param[in]	Pts			The vertices of the mesh.
	 *	\param[out]	Centers		The centers of the cells (Size : NumCells x Dimension).
	 *	\param[out]	Normals		The normals of the cells (Size : NumCells x Dimension).
	 */
	void ComputeCentersNormals(const VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Normals);

	/// Computes the RKHS-norm of itself.
	TScalar ComputeSelfNorm();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
	VNLMatrixType m_Centers;

	///	Matrix coordinates of the normals of the cells  (Size : NumCells x Dimension).
	VNLMatrixType m_Normals;

	/// Number of cells (i.e. triangles) of the mesh.
	int m_NumCells;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the oriented surface.
	TScalar m_NormSquared;

	/// See Landmark::m_VTKMutex for details.
	itk::SimpleFastMutexLock m_VTKMutex;


}; /* class OrientedSurfaceMesh */


#ifndef MU_MANUAL_INSTANTIATION
#include "OrientedSurfaceMesh.txx"
#endif


#endif /* _OrientedSurfaceMesh_h */
