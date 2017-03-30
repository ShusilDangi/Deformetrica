/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Landmark_h
#define _Landmark_h

#include "DeformableObject.h"

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
 *	\brief 		Landmarks (i.e. labelled point sets)
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The Landmark class inherited from DeformableObject represents a set of labelled points.
 *				This class assumes that the source and the target have the same number of points with a
 *				point-to-point correspondence.
 */
template <class TScalar, unsigned int Dimension>
class Landmark : public DeformableObject<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> Superclass;

	/// Vector type.
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

	Landmark();
	/// Copy constructor.
	Landmark(const Landmark& other);

	virtual Landmark* Clone() { return new Landmark(*this); }

	virtual ~Landmark();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the pointer on a VTK object to \e polyData.
	virtual void SetPolyData(vtkPolyData* polyData);

	/// Returns the coordinates of the points.
	inline VNLMatrixType GetPointCoordinates() const { return m_PointCoordinates; }
	/// Sets the coordinates of the points to \e Y.
	void SetPointCoordinates(const VNLMatrixType& Y);

	/// Returns the number of points of the deformable object.
	inline int GetNumberOfPoints() const { return m_NumberOfPoints; }

	/// Returns the minimum of the coordinates.
	inline VNLVectorType GetMin() const { return m_Min; }

	/// Returns the maximum of the coordinates.
	inline VNLVectorType GetMax() const { return m_Max; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void Update();

	/// See DeformableObject::SetDiffeoAndPointTrajectory() for \b warnings.
	/// \warning	Use with care, it is up to the user to ensure that the point trajectory is consistent
	///				with m_PointCoordinates and with the input diffeos.
	/// this->SetDiffeos(def); this->Update(); should be used whenever possible
	inline void SetDiffeoAndPointTrajectory(DiffeosType* def, VNLMatrixList& traj)
		{ this->SetDiffeos(def); m_PointsT = traj; this->SetUnModified(); }


	virtual VNLMatrixList GetFlow();

	virtual TScalar ComputeMatch(Superclass* target);
	virtual TScalar ComputeMatch(Superclass* target, unsigned int t);

	virtual VNLMatrixType ComputeMatchGradient(Superclass* target);
	virtual VNLMatrixType ComputeMatchGradient(Superclass* target, unsigned int t);

	virtual void TransportAlongGeodesic(VNLMatrixType& Eta0, VNLMatrixList& EtaT, VNLMatrixList& YT);

	virtual void WriteDeformedObjectAt(unsigned int t, std::string filename);

	/// Computes velocity at time \e t.
	virtual VNLMatrixType ComputeVelocityAt(unsigned int t);

protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the point coordinates when Landmark::SetPolyData(vtkPolyData* polyData) is called.
	void UpdatePointCoordinates();

	/// Compute the trajectory of the deformable object (automatically called by Update() method).
	void UpdatePointTrajectory();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Pointer on a geometric structure namely a VTK object.
	vtkSmartPointer<vtkPolyData> m_PolyData;

	///	Matrix coordinates of the points (Size : NumberOfPoints x Dimension).
	VNLMatrixType m_PointCoordinates;

	///	List of different matrices containing positions at different time steps.
	VNLMatrixList m_PointsT;

	///	Number of points of the deformable object.
	int m_NumberOfPoints;

	///	Vector containing the minimum of the coordinates of the landmark.
	VNLVectorType m_Min;
	///	Vector containing the minimum of the coordinates of the landmark.
	VNLVectorType m_Max;

	///	Object used to perform mutex (mutual exclusion) with m_PolyData (important for multithreaded programming).
	itk::SimpleFastMutexLock m_VTKMutex;


}; /* class Landmark */


#ifndef MU_MANUAL_INSTANTIATION
#include "Landmark.txx"
#endif


#endif /* _Landmark_h */
