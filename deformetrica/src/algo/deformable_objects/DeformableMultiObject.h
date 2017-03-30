/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableMultiObject_h
#define _DeformableMultiObject_h

#include "Diffeos.h"
#include "DeformableObject.h"

#include "vnl/vnl_matrix.h"

#include <vector>

/**
 *	\brief 		Collections of deformable objects (i.e. multi-object)
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The DeformableMultiObject class is used to deal with collections of deformable objects embedded in the current 2D or 3D space.
 *				It extends for such collections the methods of the deformable object class that are designed for a single object at a time.
 */
template <class TScalar, unsigned int Dimension>
class DeformableMultiObject
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef vnl_vector<TScalar> VNLVectorType;
	/// Matrix type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// List of matrices type.
	typedef std::vector<VNLMatrixType> VNLMatrixList;

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	/// List of deformable objects type.
	typedef std::vector<DeformableObjectType*> DeformableObjectList;

	/// Deformation type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableMultiObject();
	/// Copy constructor.
	DeformableMultiObject(const DeformableMultiObject& other);

	/// Makes a copy of the object.
	DeformableMultiObject* Clone() { return new DeformableMultiObject(*this); }

	~DeformableMultiObject();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the list of the deformable objects.
	inline DeformableObjectList GetObjectList() { return m_ObjectList; }
	/// Sets the list of deformable objects \e to list.
	void SetObjectList(DeformableObjectList& list);

	/// Returns the coordinates of the vertices of the deformable objects (for objects of Landmark kind).
	VNLMatrixList GetImageAndPointData();
	/// Sets the coordinates of the vertices of the deformable objects \e Y (for objects of Landmark kind).
	void SetImageAndPointData(VNLMatrixList & Y);

	/// Sets the deformation to \e def.
	inline void SetDiffeos(DiffeosType* def) { m_Def = def; m_Modified = true; }

	/// Returns a vector of each \f$\sigma^2\f$ of the different deformable objects.
	std::vector<TScalar> GetDataSigmaSquared();
	/// Sets each \f$\sigma^2\f$ of the different deformable objects to \e d[i].
	void SetDataSigmaSquared(std::vector<TScalar> d);

	/// Returns true if the parameters of the deformable objects have changed, false otherwise.
	inline bool IsModified() { return m_Modified; }
	/// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
	inline void SetModified() { m_Modified = true; }
	/// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
	inline void SetUnModified() { m_Modified = false; }

	/// Returns the total number of objects.
	inline int GetNumberOfObjects() { return m_NumberOfObjects; }

	/// Returns the number of objects of landmark kind.
	inline int GetNumberOfLandmarkKindObjects() { return m_NumberOfLandmarkKindObjects; } 

	/// Returns the sum of the number of points of each deformable object.
	int GetTotalNumberOfPoints();

	/// See DeformableObject::OutOfBox() for details.
	inline bool OutOfBox() { return m_OutOfBox; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// According to m_OutBox, call the method DeformableObject::Update() for each deformable object.
	void Update();

	/// Returns the list of deformed objects at final time.
	DeformableMultiObject* GetDeformedObject();

	/// Computes DeformableObject::ComputeMatch(DeformableObject* target) for each deformable object.
	std::vector<TScalar> ComputeMatch(DeformableMultiObject* target);
	/// Computes DeformableObject::ComputeMatch(DeformableObject* target, unsigned int t) for each deformable object.
	std::vector<TScalar> ComputeMatch(DeformableMultiObject* target, unsigned int t);


	/// Computes DeformableObject::ComputeMatchGradient(DeformableObject* target) for each deformable object.
	VNLMatrixList ComputeMatchGradient(DeformableMultiObject* target);
	/// Computes DeformableObject::ComputeMatchGradient(DeformableObject* target, unsigned int t) for each deformable object.
	VNLMatrixList ComputeMatchGradient(DeformableMultiObject* target, unsigned int t);

	/// Computes DeformableObject::TransportAlongGeodesic(VNLMatrixType& Eta0, VNLMatrixList& EtaT, VNLMatrixList& YT) for each deformable object.
	void TransportAlongGeodesic(VNLMatrixList& InitialConditions, VNLMatrixList& VectorsT, VNLMatrixList& PointsT);

	/// Calls the method DeformableObject::WriteDeformableObjectAt() for each deformable object.
	void WriteDeformedMultiObjectAt(unsigned int t, std::vector<std::string>& str);
	/// Calls the method DeformableObject::WriteFlow() for each deformable object.
	void WriteFlow(std::vector<std::string>& name, std::vector<std::string>& extension);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates trajectory if the deformable object is of Landmark type.
	void UpdateLandmarkPointsTrajectory();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Collection of m_NumberOfObjects deformable objects.
	DeformableObjectList m_ObjectList;

	/// Total number of objects.
	int m_NumberOfObjects;

	/// Number of objects of landmark kind.
	int m_NumberOfLandmarkKindObjects;

	/// Vector containing at each cell the number of points of the associated deformable object.
	std::vector<int> m_NumberOfPoints;

	/// Entity representing the deformation.
	DiffeosType* m_Def;

	/// See DeformableObject::m_Modified for details.
	bool m_Modified;

	/// See DeformableObject::m_OutOfBox for details.
	bool m_OutOfBox;


}; /* class DeformableMultiObject */


#ifndef MU_MANUAL_INSTANTIATION
#include "DeformableMultiObject.txx"
#endif


#endif /* _DeformableMultiObject_h */
