/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _OrientedPolyLine_txx
#define _OrientedPolyLine_txx

#include "OrientedPolyLine.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include <cstring>
#include <iostream>
#include <sstream>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::OrientedPolyLine() : Superclass()
 {
	this->SetOrientedPolyLineType();
	this->SetModified();
	m_KernelWidth = 0;
 }



template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::~OrientedPolyLine()
{

}



template <class TScalar, unsigned int Dimension>
OrientedPolyLine<TScalar, Dimension>
::OrientedPolyLine(const OrientedPolyLine& other) : Superclass(other)
 {
	this->SetOrientedPolyLineType();

	m_Centers = other.m_Centers;
	m_Tangents = other.m_Tangents;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_NormSquared = other.m_NormSquared;

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
 {
	Superclass::SetPolyData(polyData);

	m_NumCells = polyData->GetNumberOfCells();

	VNLMatrixType Pts = this->GetPointCoordinates();	
	this->ComputeCentersTangents(Pts,m_Centers,m_Tangents);

	// overload m_Min and m_Max using the centers instead of the points (useful when using the results of matching pursuit as target)
	VNLVectorType Min(Dimension, 1e20);
	VNLVectorType Max(Dimension, -1e20);

	for (int i = 0; i < m_NumCells; i++)
	{
		VNLVectorType p = m_Centers.get_row(i);
		for (unsigned int dim=0; dim<Dimension; dim++)
		{	
			Min[dim] = (p[dim] < Min[dim])?p[dim]:Min[dim];
			Max[dim] = (p[dim] > Max[dim])?p[dim]:Max[dim];
		}
	}

	Superclass::m_Min = Min;
	Superclass::m_Max = Max;

	this->SetModified();

	if (m_KernelWidth != 0)
		m_NormSquared = this->ComputeSelfNorm();

 }



template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::SetKernelWidth(TScalar h)
 {
	m_KernelWidth = h;

	if (Superclass::m_NumberOfPoints != 0)
		m_NormSquared = this->ComputeSelfNorm();
 }




////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
TScalar
OrientedPolyLine<TScalar, Dimension>
::ComputeMatch(OrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	OrientedPolyLine* targetOrientedPolyLine = dynamic_cast<OrientedPolyLine*>(target);

	if (m_KernelWidth != targetOrientedPolyLine->GetKernelWidth())
	{
		std::cerr << "DEBUG kw = "  << m_KernelWidth << " target kw = " << targetOrientedPolyLine->GetKernelWidth() << std::endl;
		throw std::runtime_error("Kernel width of curve currents mismatched");
	}

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetOrientedPolyLine->GetCenters();
	VNLMatrixType targTangents = targetOrientedPolyLine->GetTangents();

	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];
	VNLMatrixType defCenters, defTangents;
	this->ComputeCentersTangents(defPts, defCenters, defTangents);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defTangents);

	VNLMatrixType SdotS = kernelObject->Convolve(defCenters);

	TScalar match = targetOrientedPolyLine->GetNormSquared();

	// std::cout << "targetNormSquared = " << match << std::endl;

	for (int i = 0; i < m_NumCells; i ++)
		match += dot_product(SdotS.get_row(i), defTangents.get_row(i));

	// TScalar SourceNormSquared = match -  targetOrientedPolyLine->GetNormSquared();
	// std::cout << "targetNormSquared = " << SourceNormSquared << std::endl;

	VNLMatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetOrientedPolyLine->GetNumberOfCells(); i++)
		match -= 2.0 * dot_product( SdotT.get_row(i),  targTangents.get_row(i) );

	// std::cout << "Inner Product = " << match - SourceNormSquared - targetOrientedPolyLine->GetNormSquared()<< std::endl;
	// 
	// std::cout << "Total = " << match << std::endl;

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
TScalar
OrientedPolyLine<TScalar, Dimension>
::ComputeMatch(DeformableObjectType* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	OrientedPolyLine* targetOrientedPolyLine = dynamic_cast<OrientedPolyLine*>(target);

	if (m_KernelWidth != targetOrientedPolyLine->GetKernelWidth())
	{
		std::cerr << "DEBUG kw = "  << m_KernelWidth << " target kw = " << targetOrientedPolyLine->GetKernelWidth() << std::endl;
		throw std::runtime_error("Kernel width of curve currents mismatched");
	}

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetOrientedPolyLine->GetCenters();
	VNLMatrixType targTangents = targetOrientedPolyLine->GetTangents();

	// Make sure t is inside the time interval
	if (t >= numTimePts)
	{
		std::cerr << "time for matching is outside of time interval" << std::endl;
	}


	VNLMatrixType defPts = Superclass::m_PointsT[t];
	VNLMatrixType defCenters, defTangents;
	this->ComputeCentersTangents(defPts, defCenters, defTangents);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defTangents);

	VNLMatrixType SdotS = kernelObject->Convolve(defCenters);

	TScalar match = targetOrientedPolyLine->GetNormSquared();

	// std::cout << "targetNormSquared = " << match << std::endl;

	for (int i = 0; i < m_NumCells; i ++)
		match += dot_product(SdotS.get_row(i), defTangents.get_row(i));

	// TScalar SourceNormSquared = match -  targetOrientedPolyLine->GetNormSquared();
	// std::cout << "targetNormSquared = " << SourceNormSquared << std::endl;

	VNLMatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetOrientedPolyLine->GetNumberOfCells(); i++)
		match -= 2.0 * dot_product( SdotT.get_row(i),  targTangents.get_row(i) );

	// std::cout << "Inner Product = " << match - SourceNormSquared - targetOrientedPolyLine->GetNormSquared()<< std::endl;
	//
	// std::cout << "Total = " << match << std::endl;

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename OrientedPolyLine<TScalar, Dimension>::VNLMatrixType
OrientedPolyLine<TScalar, Dimension>
::ComputeMatchGradient(OrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	OrientedPolyLine* targetOrientedPolyLine = dynamic_cast<OrientedPolyLine*>(target);

	if (m_KernelWidth != targetOrientedPolyLine->GetKernelWidth())
		throw std::runtime_error("Kernel width of curve currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetOrientedPolyLine->GetCenters();
	VNLMatrixType targTangents = targetOrientedPolyLine->GetTangents();

	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];
	VNLMatrixType defCenters, defTangents;
	this->ComputeCentersTangents(defPts, defCenters, defTangents);


	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defTangents);
	VNLMatrixType KtauS = kernelObject->Convolve(defCenters);
	VNLMatrixType gradKtauS = kernelObject->ConvolveGradient(defCenters, defTangents);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targTangents);
	VNLMatrixType KtauT = kernelObject->Convolve(defCenters);
	VNLMatrixType gradKtauT = kernelObject->ConvolveGradient(defCenters, defTangents);

	VNLMatrixType gradKtau = gradKtauS - gradKtauT;

	VNLMatrixType gradmatch(this->GetNumberOfPoints(),Dimension,0.0);
	for (int f = 0; f < m_NumCells; f++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(f, ptIds);
		m_VTKMutex.Unlock();

		int indM = ptIds->GetId(0);
		int indP = ptIds->GetId(1);

		VNLVectorType Ktau = KtauS.get_row(f) - KtauT.get_row(f);
		Ktau *= 2.0;
		gradmatch.set_row(indM, gradmatch.get_row(indM) - Ktau);
		gradmatch.set_row(indP, gradmatch.get_row(indP) + Ktau);

		gradmatch.set_row(indM, gradmatch.get_row(indM) + gradKtau.get_row(f));
		gradmatch.set_row(indP, gradmatch.get_row(indP) + gradKtau.get_row(f));

	}

	delete kernelObject;

	return gradmatch;
 }



template <class TScalar, unsigned int Dimension>
typename OrientedPolyLine<TScalar, Dimension>::VNLMatrixType
OrientedPolyLine<TScalar, Dimension>
::ComputeMatchGradient(OrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	OrientedPolyLine* targetOrientedPolyLine = dynamic_cast<OrientedPolyLine*>(target);

	if (m_KernelWidth != targetOrientedPolyLine->GetKernelWidth())
		throw std::runtime_error("Kernel width of curve currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetOrientedPolyLine->GetCenters();
	VNLMatrixType targTangents = targetOrientedPolyLine->GetTangents();

	// Make sure t is inside the time interval
	if (t >= numTimePts)
	{
		std::cerr << "time for matching is outside of time interval" << std::endl;
	}

	VNLMatrixType defPts = Superclass::m_PointsT[t];
	VNLMatrixType defCenters, defTangents;
	this->ComputeCentersTangents(defPts, defCenters, defTangents);


	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defTangents);
	VNLMatrixType KtauS = kernelObject->Convolve(defCenters);
	VNLMatrixList gradKtauS = kernelObject->ConvolveGradient(defCenters);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targTangents);
	VNLMatrixType KtauT = kernelObject->Convolve(defCenters);
	VNLMatrixList gradKtauT = kernelObject->ConvolveGradient(defCenters);

	VNLMatrixType gradmatch(this->GetNumberOfPoints(),Dimension,0.0);

	for (int f = 0; f < m_NumCells; f++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(f, ptIds);
		m_VTKMutex.Unlock();

		int indM = ptIds->GetId(0);
		int indP = ptIds->GetId(1);

		VNLVectorType Ktau = KtauS.get_row(f) - KtauT.get_row(f);
		Ktau *= 2.0;
		gradmatch.set_row(indM, gradmatch.get_row(indM) - Ktau);
		gradmatch.set_row(indP, gradmatch.get_row(indP) + Ktau);

		VNLMatrixType delta = gradKtauS[f] - gradKtauT[f];
		VNLVectorType gradKtau = delta.transpose() * defTangents.get_row(f);
		gradmatch.set_row(indM, gradmatch.get_row(indM) + gradKtau);
		gradmatch.set_row(indP, gradmatch.get_row(indP) + gradKtau);

	}

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedPolyLine<TScalar, Dimension>
::ComputeCentersTangents(const VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Tangents)
 {
	Centers.set_size(m_NumCells, Dimension);
	Tangents.set_size(m_NumCells, Dimension);

	// Superclass::m_PolyData->BuildCells();

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(i, ptIds);
		m_VTKMutex.Unlock();

		if (ptIds->GetNumberOfIds() != 2)
			throw std::runtime_error("Not a polygonal line!");

		VNLVectorType p0 = Pts.get_row(ptIds->GetId(0));
		VNLVectorType p1 = Pts.get_row(ptIds->GetId(1));
		Centers.set_row(i, (p0 + p1) / 2.0 );
		Tangents.set_row(i, p1 - p0);

	}
 }



template <class TScalar, unsigned int Dimension>
TScalar
OrientedPolyLine<TScalar, Dimension>
::ComputeSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_Tangents);

	VNLMatrixType selfKW = kernelObject->Convolve(m_Centers);

	TScalar norm2 = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
		norm2 += dot_product(selfKW.get_row(i), m_Tangents.get_row(i));

	delete kernelObject;

	return norm2;
 }




#endif /* _OrientedPolyLine_txx */
