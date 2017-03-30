/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NonOrientedSurfaceMesh_h
#define _NonOrientedSurfaceMesh_h

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
 *	\brief 		Non oriented surface meshes.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The NonOrientedSurfaceMesh class inherited from Landmark represents a triangular mesh.
 *				This class uses the varifold representation of surfaces, which is insensitive to the local orientation of the mesh.
 */
template <class TScalar, unsigned int Dimension>
class NonOrientedSurfaceMesh : public Landmark<TScalar, Dimension>
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

	NonOrientedSurfaceMesh();
	/// Copy constructor.
	NonOrientedSurfaceMesh(const NonOrientedSurfaceMesh& other);

	virtual NonOrientedSurfaceMesh* Clone() { return new NonOrientedSurfaceMesh(*this); }

	virtual ~NonOrientedSurfaceMesh();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void SetPolyData(vtkPolyData* polyData);

	/// Returns the centers of the cells.
	inline VNLMatrixType GetCenters() const { return m_Centers; }

	/// Returns the normals of the cells.
	inline VNLMatrixType GetNormals() const { return m_Normals; }

	/// Returns the matrices of the type \f[ \frac{n_i n_j^T}{\left|n_i\right|^{1/2}\left[n_j\right|^{1/2}} \f], where \f$ n_i \f$ denotes the normals.
	inline VNLMatrixType GetMatrixNormals() const { return m_MatrixNormals; }

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

	/**
	 *	\brief		Computes the centers and the normals from the points.
	 *
	 *	\details	Given a triangular mesh, this method computes from the vertices of the mesh
	 *				the centers and the normals of each face.
	 *
	 *	\param[in]	Pts				The vertices of the mesh.
	 *	\param[out]	Centers			The centers of the cells (Size : NumCells x Dimension).
	 *	\param[out]	Normals			The normals of the cells (Size : NumCells x Dimension).
	 // *	\param[out]	MatrixNormals	The matrices \f[ \frac{n_i n_j^T}{\left|n_i\right|^{1/2}\left|n_j\right|^{1/2}} \f], where \f$ n_i \f$ denotes the normals.
	 */
	void ComputeCentersNormals(VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Normals, VNLMatrixType& MatrixNormals);

	/// Computes the RKHS-norm of itself.
	TScalar ComputeSelfNorm();

	/**
	 *	\brief		computes the product between with symmetric matrices and vectors.
	 *
	 *	\details
	 *
	 *	\param[in]	M	symmetric matrix .
	 *	\param[in]	X	vector .
	 *
	 *	\return		product .
	 */
	VNLVectorType special_product(VNLVectorType M, VNLVectorType X);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
	VNLMatrixType m_Centers;

	///	Matrix coordinates of the normals of the cells  (Size : NumCells x Dimension).
	VNLMatrixType m_Normals;

	/// The matrix \f[ \frac{n_i n_j^T}{\left|n_i\right|^{1/2}\left|n_j\right|^{1/2}} \f], where \f$ n_i \f$ denotes the normals.
	VNLMatrixType m_MatrixNormals;

	/// Number of cells (i.e. triangles) of the mesh.
	int m_NumCells;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	/// Squared RKHS-norm of the oriented surface.
	TScalar m_NormSquared;

	/// See Landmark::m_VTKMutex for details.
	itk::SimpleFastMutexLock m_VTKMutex;

}; /* class NonOrientedSurfaceMesh */


#ifndef MU_MANUAL_INSTANTIATION
#include "NonOrientedSurfaceMesh.txx"
#endif


#endif /* _NonOrientedSurfaceMesh_h */
