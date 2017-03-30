/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableMultiObject_txx
#define _DeformableMultiObject_txx

#include "DeformableMultiObject.h"
#include "DeformableObject.h"
#include "Landmark.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#include "KernelFactory.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::DeformableMultiObject()
 {
	m_OutOfBox = false;
	m_Modified = true;
 }



template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::~DeformableMultiObject()
{ 
	for (unsigned int i = 0; i < m_ObjectList.size(); i++)
		delete m_ObjectList[i];
}



template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>
::DeformableMultiObject(const DeformableMultiObject& other)
 {
	// Deep copy so we can safely modify the objects
	m_ObjectList.resize(other.m_ObjectList.size());
	for (unsigned int i = 0; i < other.m_ObjectList.size(); i++)
		m_ObjectList[i] = other.m_ObjectList[i]->Clone();

	m_NumberOfObjects = other.m_NumberOfObjects;
	m_NumberOfLandmarkKindObjects = other.m_NumberOfLandmarkKindObjects;

	m_Def = other.m_Def;

	m_NumberOfPoints = other.m_NumberOfPoints;

	m_OutOfBox = other.m_OutOfBox;
	m_Modified = other.m_Modified;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::SetObjectList(DeformableObjectList& objList)
 {
	// Deep copy so we can safely modify the objects
	m_ObjectList.resize(objList.size());
	for (unsigned int i = 0; i < objList.size(); i++)
		m_ObjectList[i] = objList[i]->Clone();

	m_NumberOfObjects = m_ObjectList.size();

	m_NumberOfLandmarkKindObjects = 0;

	m_NumberOfPoints.resize(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		DeformableObjectType* obj = objList[i];

		if ( obj->IsOfLandmarkKind() )
		{
			m_NumberOfPoints[i] = obj->GetNumberOfPoints(); 
			m_NumberOfLandmarkKindObjects += 1;
		}
	}

	m_Modified = true;
 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::SetDataSigmaSquared(std::vector<TScalar> ds)
 {
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->SetDataSigmaSquared(ds[i]);

 }



template <class TScalar, unsigned int Dimension>
std::vector<TScalar>
DeformableMultiObject<TScalar, Dimension>
::GetDataSigmaSquared()
 {
	std::vector<TScalar> ds(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		ds[i] = m_ObjectList[i]->GetDataSigmaSquared();

	return ds;
 }



template<class TScalar, unsigned int Dimension>
int DeformableMultiObject<TScalar, Dimension>
::GetTotalNumberOfPoints()
 {
	int numberOfTotalPoints = 0;
	for (int i=0; i<m_NumberOfObjects; i++)
	{
		numberOfTotalPoints += m_NumberOfPoints[i];
	}

	return numberOfTotalPoints;
 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::SetImageAndPointData(VNLMatrixList& Y)
 {
	if (Y.size() != m_NumberOfObjects)
		throw std::runtime_error("size of list and number of objects mismatch");

	typedef Landmark<TScalar, Dimension> LandmarkType;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			LandmarkType* obj = dynamic_cast<LandmarkType*>(m_ObjectList[i]);
			obj->SetPointCoordinates(Y[i]);
		}
		else
			throw std::runtime_error("Unknown object type");
	}
 }



template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::VNLMatrixList
DeformableMultiObject<TScalar, Dimension>
::GetImageAndPointData()
 {
	VNLMatrixList Y(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			typedef Landmark<TScalar, Dimension> LandmarkType;
			LandmarkType* obj = dynamic_cast<LandmarkType*>(m_ObjectList[i]);
			Y[i] = obj->GetPointCoordinates();
		}
		else
			throw std::runtime_error("Unknown object type");
	}

	return Y;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::Update()
 {
	if (this->IsModified())
	{
		m_OutOfBox = 0;
		if (m_NumberOfLandmarkKindObjects)
			this->UpdateLandmarkPointsTrajectory(); // Flow + obj[i]->SetPointTrajectory()
	}

	this->SetUnModified();
 }



template <class TScalar, unsigned int Dimension>
DeformableMultiObject<TScalar, Dimension>*
DeformableMultiObject<TScalar, Dimension>
::GetDeformedObject()
 {

	DeformableObjectList defObj(m_NumberOfObjects);

	for (int i = 0; i < m_NumberOfObjects; i++)
		defObj[i] = m_ObjectList[i]->GetDeformedObject();

	DeformableMultiObject* out = new DeformableMultiObject();
	out->SetObjectList(defObj);

	return out;
 }



template <class TScalar, unsigned int Dimension>
typename std::vector<TScalar>
DeformableMultiObject<TScalar, Dimension>
::ComputeMatch(DeformableMultiObject* target)
 {
	if (m_NumberOfObjects != target->GetNumberOfObjects())
	{
		std::cerr << "number of objects mismatched" << std::endl;
	}

	DeformableObjectList targetList = target->GetObjectList();

	std::vector<TScalar> match(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		match[i] = m_ObjectList[i]->ComputeMatch(targetList[i]);

	return match;
 }



template<class TScalar, unsigned int Dimension>
typename std::vector<TScalar> DeformableMultiObject<TScalar, Dimension>
::ComputeMatch(DeformableMultiObject* target, unsigned int t)
 {
	if (m_NumberOfObjects != target->GetNumberOfObjects())
	{
		std::cerr << "number of objects mismatched" << std::endl;
	}

	DeformableObjectList targetList = target->GetObjectList();

	std::vector<TScalar> match(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		match[i] = m_ObjectList[i]->ComputeMatch(targetList[i], t);

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::VNLMatrixList
DeformableMultiObject<TScalar, Dimension>
::ComputeMatchGradient(DeformableMultiObject* target)
 {
	if (m_NumberOfObjects != target->GetNumberOfObjects())
		throw std::runtime_error("number of objects mismatched");

	DeformableObjectList targetList = target->GetObjectList();
	VNLMatrixList gradmatch(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		gradmatch[i] = m_ObjectList[i]->ComputeMatchGradient(targetList[i]);

	return gradmatch;
 }



template<class TScalar, unsigned int Dimension>
typename DeformableMultiObject<TScalar, Dimension>::VNLMatrixList DeformableMultiObject<TScalar, Dimension>
::ComputeMatchGradient(DeformableMultiObject* target, unsigned int t)
 {
	if (m_NumberOfObjects != target->GetNumberOfObjects())
		throw std::runtime_error("number of objects mismatched");

	DeformableObjectList targetList = target->GetObjectList();
	VNLMatrixList gradmatch(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		gradmatch[i] = m_ObjectList[i]->ComputeMatchGradient(targetList[i], t);

	return gradmatch;
 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::TransportAlongGeodesic(VNLMatrixList& InitialConditions, VNLMatrixList& VectorsT, VNLMatrixList& PointsT)
 {
	if (InitialConditions.size() != m_NumberOfObjects)
		throw std::runtime_error("number of initial conditions should equal the number of objects");

	int NumbTimePoints = m_Def->GetNumberOfTimePoints();
	VectorsT.resize(NumbTimePoints);
	PointsT.resize(NumbTimePoints);

	int NumbPts = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
		NumbPts += m_NumberOfPoints[i];

	VNLMatrixType zeroM(NumbPts, Dimension, 0.0);
	for (int t = 0; t < NumbTimePoints; t++)
	{
		VectorsT[t] = zeroM;
		PointsT[t] = zeroM;
	}

	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		VNLMatrixList VecT;
		VNLMatrixList PtsT;
		m_ObjectList[i]->TransportAlongGeodesic(InitialConditions[i], VecT, PtsT);

		for (int t = 0; t < NumbTimePoints; t++)
			for (int r=0; r < m_NumberOfPoints[i]; r++)
			{
				PointsT[t].set_row(counter + r, PtsT[t].get_row(r));
				VectorsT[t].set_row(counter + r,VecT[t].get_row(r));
			}

		counter += m_NumberOfPoints[i];
	}	

 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::WriteDeformedMultiObjectAt(unsigned int t, std::vector<std::string>& outfn)
 {
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->WriteDeformedObjectAt(t, outfn[i]);	
 }



template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::WriteFlow(std::vector<std::string>& outfn, std::vector<std::string>& outext)
 {
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_ObjectList[i]->WriteFlow(outfn[i], outext[i]);
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableMultiObject<TScalar, Dimension>
::UpdateLandmarkPointsTrajectory()
 {
	// Compute the number of points
	int NumbPoints = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
			NumbPoints += m_NumberOfPoints[i];

	typedef Landmark<TScalar, Dimension> LandmarkType;

	// Put the point data into VTK format
	vtkSmartPointer<vtkPoints> X0 = vtkSmartPointer<vtkPoints>::New();
	X0->SetNumberOfPoints(NumbPoints);
	// VNLMatrixType X0(NumbPoints,Dimension);
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			VNLMatrixType Xaux = dynamic_cast<LandmarkType*>(m_ObjectList[i])->GetPointCoordinates();
			for (int r = 0; r < m_NumberOfPoints[i] ; r++)
			{
				// X0.set_row(counter + r, Xaux.get_row(r));
				double p[Dimension];
				for (int dim = 0; dim < Dimension; dim++)
					p[dim] = Xaux(r,dim);

				X0->SetPoint(counter+r, p);
			}

			counter += m_NumberOfPoints[i];
		}

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->SetPoints(X0);

	LandmarkType* allPointData = new LandmarkType();
	allPointData->SetPolyData(pd);
	allPointData->SetDiffeos(m_Def);
	allPointData->Update();

	VNLMatrixList XT = allPointData->GetFlow();
	m_OutOfBox += allPointData->OutOfBox();

	int NumbTimePts = m_Def->GetNumberOfTimePoints();

	counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
		if ( m_ObjectList[i]->IsOfLandmarkKind() )
		{
			VNLMatrixList Xi(NumbTimePts);
			for (int t = 0; t < NumbTimePts; t++)
				Xi[t] = XT[t].get_n_rows(counter, m_NumberOfPoints[i]);

			m_ObjectList[i]->SetDiffeoAndPointTrajectory(m_Def, Xi);
			counter += m_NumberOfPoints[i];
		}

	delete allPointData;
 }


#endif /* _DeformableMultiObject_txx */
