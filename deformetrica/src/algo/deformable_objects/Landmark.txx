/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Landmark_txx
#define _Landmark_txx

#include "Landmark.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkCellData.h"
#include "vtkPointData.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "SimpleTimer.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::Landmark() : Superclass()
 {
	this->SetLandmarkType();
	this->SetModified();
	m_NumberOfPoints = 0;
 }



template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::~Landmark()
{

}



template <class TScalar, unsigned int Dimension>
Landmark<TScalar, Dimension>
::Landmark(const Landmark& other) : Superclass(other)
 {
	this->SetLandmarkType();

	m_PolyData = vtkSmartPointer<vtkPolyData>::New();
	m_PolyData->DeepCopy(other.m_PolyData);

	m_PointCoordinates = other.m_PointCoordinates;
	m_PointsT = other.m_PointsT;

	m_NumberOfPoints = other.m_NumberOfPoints;

	m_Min = other.m_Min;
	m_Max = other.m_Max;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
 {
	m_PolyData = polyData;
	m_NumberOfPoints = m_PolyData->GetNumberOfPoints();

	this->UpdatePointCoordinates();

	VNLVectorType Min(Dimension, 1e20);
	VNLVectorType Max(Dimension, -1e20);

	for (int i = 0; i < m_NumberOfPoints; i++)
	{
		VNLVectorType p = m_PointCoordinates.get_row(i);
		for (unsigned int dim=0; dim<Dimension; dim++)
		{	
			Min[dim] = (p[dim] < Min[dim])?p[dim]:Min[dim];
			Max[dim] = (p[dim] > Max[dim])?p[dim]:Max[dim];
		}
	}

	m_Min = Min;
	m_Max = Max;

	this->SetModified(); 
 }



template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::SetPointCoordinates(const VNLMatrixType& Y)
 {
	if (Y.rows() != m_NumberOfPoints)
		throw std::runtime_error("number of points mismatched");
	if (Y.columns() != Dimension)
		throw std::runtime_error("Dimension mismatched");

	for (unsigned int i = 0; i < m_NumberOfPoints; i++)
	{
		double p[3];
		for (int dim = 0; dim < Dimension; dim++)
			p[dim] = Y(i, dim);

		m_VTKMutex.Lock();
		m_PolyData->GetPoints()->SetPoint(i, p);
		m_VTKMutex.Unlock();
	}

	m_PointCoordinates = Y;

	this->SetModified();
 }




////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::Update()
 {
	if (this->isModified())
		this->UpdatePointTrajectory();

	this->SetUnModified();
 }


template <class TScalar, unsigned int Dimension>
typename Landmark<TScalar, Dimension>::VNLMatrixList
Landmark<TScalar, Dimension>
::GetFlow()
 {
	 return m_PointsT;
 }



template <class TScalar, unsigned int Dimension>
TScalar
Landmark<TScalar, Dimension>
::ComputeMatch(Landmark<TScalar, Dimension>::Superclass* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	Landmark* targetLandmark = dynamic_cast<Landmark*>(target);

	if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
		throw std::runtime_error("Landmark object should have the same number of points");

	unsigned int numTimePts = Superclass::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targPts = targetLandmark->GetPointCoordinates();	
	VNLMatrixType defPts = m_PointsT[numTimePts-1];

	VNLMatrixType D = defPts - targPts;

	TScalar match = D.frobenius_norm();
	match *= match;

	return match;
 }



template<class TScalar, unsigned int Dimension>
TScalar Landmark<TScalar, Dimension>
::ComputeMatch(Landmark<TScalar, Dimension>::Superclass* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	Landmark* targetLandmark = dynamic_cast<Landmark*>(target);

	if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
		throw std::runtime_error("Landmark object should have the same number of points");

	unsigned int numTimePts = Superclass::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targPts = targetLandmark->GetPointCoordinates();
	VNLMatrixType defPts = m_PointsT[t];

	VNLMatrixType D = defPts - targPts;

	TScalar match = D.frobenius_norm();
	match *= match;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename Landmark<TScalar, Dimension>::VNLMatrixType
Landmark<TScalar, Dimension>
::ComputeMatchGradient(Landmark<TScalar, Dimension>::Superclass* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	Landmark* targetLandmark = dynamic_cast<Landmark*>(target);

	if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
		throw std::runtime_error("Landmarks object should have the same number of points");

	unsigned int numTimePts = Superclass::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targPts = targetLandmark->GetPointCoordinates();
	VNLMatrixType defPts = m_PointsT[numTimePts-1];

	VNLMatrixType gradMatch = defPts - targPts;
	gradMatch *= 2.0;

	return gradMatch;
 }



template<class TScalar, unsigned int Dimension>
typename Landmark<TScalar, Dimension>::VNLMatrixType Landmark<TScalar, Dimension>
::ComputeMatchGradient(Landmark<TScalar, Dimension>::Superclass* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	Landmark* targetLandmark = dynamic_cast<Landmark*>(target);

	if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
		throw std::runtime_error("Landmarks object should have the same number of points");

	unsigned int numTimePts = Superclass::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targPts = targetLandmark->GetPointCoordinates();
	VNLMatrixType defPts = m_PointsT[numTimePts - 1];

	VNLMatrixType gradMatch = defPts - targPts;
	gradMatch *= 2.0;

	return gradMatch;
 }

template <class TScalar, unsigned int Dimension>
typename Landmark<TScalar, Dimension>::VNLMatrixType
Landmark<TScalar, Dimension>
::ComputeVelocityAt(unsigned int t)
{
	this->Update();
	
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject();
	kernelObj->SetKernelWidth(Superclass::m_Def->GetKernelWidth());

	VNLMatrixList posT = Superclass::m_Def->GetTrajectoryPositions();
	VNLMatrixList momT = Superclass::m_Def->GetTrajectoryMomentas();

	kernelObj->SetSources(posT[t]);
	kernelObj->SetWeights(momT[t]);
	
	VNLMatrixType V = kernelObj->Convolve(m_PointsT[t]);

	delete kernelObj;
	
	return V;
}


template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::TransportAlongGeodesic(VNLMatrixType& Theta1, VNLMatrixList& ThetaT, VNLMatrixList& XT)
 {
	 XT = m_PointsT;

	long numTimePoints = Superclass::m_Def->GetNumberOfTimePoints();

	// Propagate theta backward
	ThetaT.resize(numTimePoints);
	ThetaT[numTimePoints-1] = Theta1;

	TScalar dt = 1.0 / (numTimePoints - 1);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject();
	kernelObj->SetKernelWidth(Superclass::m_Def->GetKernelWidth());

	VNLMatrixList PositionsT = Superclass::m_Def->GetTrajectoryPositions();
	VNLMatrixList MomentasT = Superclass::m_Def->GetTrajectoryMomentas();

	for (long t = numTimePoints-1; t > 0; t--)
	{
		kernelObj->SetSources(PositionsT[t]);
		kernelObj->SetWeights(MomentasT[t]);

		// Grad mom splatted at CP locations, with convolutions evaluated at x(t-1)
		VNLMatrixType dTheta = kernelObj->ConvolveGradient(m_PointsT[t], ThetaT[t]);

		ThetaT[t-1] = ThetaT[t] + dTheta * dt;

		// Heun's method
		if (Superclass::m_Def->ImprovedEuler())
		{
			kernelObj->SetSources(PositionsT[t-1]);
			kernelObj->SetWeights(MomentasT[t-1]);

			VNLMatrixType dTheta2 = kernelObj->ConvolveGradient(m_PointsT[t-1], ThetaT[t-1]);
			ThetaT[t-1] = ThetaT[t] + (dTheta + dTheta2) * (dt * 0.5);
		}

	}

	delete kernelObj;
 }


template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::WriteDeformedObjectAt(unsigned int t, std::string str)
{
	
	unsigned int numTimePoints = Superclass::m_Def->GetNumberOfTimePoints();
	if (t>numTimePoints)
		throw std::runtime_error("time t out of range");
	
	this->Update();

	VNLMatrixType defPts;
	VNLMatrixType Velocity(m_NumberOfPoints, Dimension);
	if (t==0)
	{
		defPts = m_PointCoordinates;
		Velocity.fill(0.0);
	}
	else
	{
		this->Update();
		defPts = m_PointsT[t];
		Velocity = this->ComputeVelocityAt(t);
	}
		
	m_VTKMutex.Lock();
	vtkSmartPointer<vtkPolyData> outData = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkDoubleArray> velocityVecs = vtkSmartPointer<vtkDoubleArray>::New();
	velocityVecs->SetNumberOfComponents(3);
	velocityVecs->SetNumberOfTuples(m_NumberOfPoints);
	velocityVecs->SetName("velocity");
	
	outData->DeepCopy(m_PolyData);
	
	m_VTKMutex.Unlock();

	for (unsigned int i = 0; i < m_NumberOfPoints; i++)
	{
		double p[3];
		double v[3];
		// In 2D, the z-coordinate is set to 0;
		p[2] = 0.0;
		v[2] = 0.0;

		for (int dim = 0; dim < Dimension; dim++)
		{
			p[dim] = defPts(i, dim);
			v[dim] = Velocity(i, dim);
		}

		m_VTKMutex.Lock();
		outData->GetPoints()->SetPoint(i, p);
		velocityVecs->SetTuple(i, v);
		m_VTKMutex.Unlock();
	}
	
	m_VTKMutex.Lock();
	outData->GetPointData()->SetVectors(velocityVecs);
	m_VTKMutex.Unlock();

	vtkSmartPointer<vtkPolyDataWriter> writer =
			vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(str.c_str());

#if (VTK_MAJOR_VERSION == 5)
	writer->SetInput(outData);
#else
	writer->SetInputData(outData);
#endif

	writer->Update();
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::UpdatePointCoordinates()
 {
	m_PointCoordinates.set_size(m_NumberOfPoints, Dimension);

	for (unsigned int i = 0; i < m_NumberOfPoints; i++)
	{
		// Here we used 3 since 2D points still have a z-coordinate that is equal to 0 in vtkPolyData.
		// This coordinate is removed in m_PointCoordinates
		double p[3];
		m_VTKMutex.Lock();
		m_PolyData->GetPoint(i, p);
		m_VTKMutex.Unlock();

		for (int dim = 0; dim < Dimension; dim++)
			m_PointCoordinates(i, dim) = p[dim];
	}

 }



template <class TScalar, unsigned int Dimension>
void
Landmark<TScalar, Dimension>
::UpdatePointTrajectory()
 {

	int NumbTimePoints = Superclass::m_Def->GetNumberOfTimePoints();
	TScalar dt = 1.0 / (NumbTimePoints - 1);
	
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject();
	kernelObj->SetKernelWidth(Superclass::m_Def->GetKernelWidth());

	VNLMatrixList posT = Superclass::m_Def->GetTrajectoryPositions();
	VNLMatrixList momT = Superclass::m_Def->GetTrajectoryMomentas();

	m_PointsT.resize(NumbTimePoints);
	for (unsigned int t = 0; t < NumbTimePoints ; t++)
		m_PointsT[t] = m_PointCoordinates;


	if (momT[0].frobenius_norm() < 1e-20)
		return;


	for (unsigned int t = 0; t < NumbTimePoints-1 ; t++)
	{
		kernelObj->SetSources(posT[t]);
		kernelObj->SetWeights(momT[t]);

		VNLMatrixType dX = kernelObj->Convolve(m_PointsT[t]);

		m_PointsT[t + 1] = m_PointsT[t] + (dX * dt);

		if (Superclass::m_Def->ImprovedEuler())
		{
			kernelObj->SetSources(posT[t + 1]);
			kernelObj->SetWeights(momT[t + 1]);

			VNLMatrixType dX2 = kernelObj->Convolve(m_PointsT[t + 1]);

			m_PointsT[t + 1] = m_PointsT[t] + (dX + dX2) * (dt * 0.5);
		}

		Superclass::m_OutOfBox = Superclass::m_Def->CheckBoundingBox(m_PointsT,t+1);
		if (Superclass::m_OutOfBox)
		{
			std::cout << "Landmark deformation: out of box at time t = " << t+1 << std::endl;
			break;
		}
				
	}

	delete kernelObj;

 }



#endif /* _Landmark_txx */
