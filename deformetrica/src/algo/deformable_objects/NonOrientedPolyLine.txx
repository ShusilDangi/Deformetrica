/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NonOrientedPolyLine_txx
#define _NonOrientedPolyLine_txx

#include "NonOrientedPolyLine.h"

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
NonOrientedPolyLine<TScalar, Dimension>
::NonOrientedPolyLine() : Superclass()
 {
	this->SetNonOrientedPolyLineType();
	this->SetModified();
	m_KernelWidth = 0;
 }



template <class TScalar, unsigned int Dimension>
NonOrientedPolyLine<TScalar, Dimension>
::NonOrientedPolyLine(const NonOrientedPolyLine& other) : Superclass(other)
 {
	this->SetNonOrientedPolyLineType();

	m_Centers = other.m_Centers;
	m_Tangents = other.m_Tangents;
	m_MatrixTangents = other.m_MatrixTangents;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_NormSquared = other.m_NormSquared;

 }



template <class TScalar, unsigned int Dimension>
NonOrientedPolyLine<TScalar, Dimension>
::~NonOrientedPolyLine()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
NonOrientedPolyLine<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
 {
	Superclass::SetPolyData(polyData);

	m_NumCells = polyData->GetNumberOfCells();

	VNLMatrixType Pts = this->GetPointCoordinates();
	this->ComputeCentersTangents(Pts,m_Centers, m_Tangents, m_MatrixTangents);

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
NonOrientedPolyLine<TScalar, Dimension>
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
NonOrientedPolyLine<TScalar, Dimension>
::ComputeMatch(NonOrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	NonOrientedPolyLine* targetNonOrientedPolyLine = dynamic_cast<NonOrientedPolyLine*>(target);

	if (m_KernelWidth != targetNonOrientedPolyLine->GetKernelWidth())
	{
		std::cerr << "DEBUG kw = "  << m_KernelWidth << " target kw = " << targetNonOrientedPolyLine->GetKernelWidth() << std::endl;
		throw std::runtime_error("Kernel width of curve currents mismatched");
	}

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetNonOrientedPolyLine->GetCenters();
	VNLMatrixType targTangents = targetNonOrientedPolyLine->GetTangents();

	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];
	VNLMatrixType defCenters, defTangents, defMatrixTangents;
	this->ComputeCentersTangents(defPts, defCenters, defTangents, defMatrixTangents);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defMatrixTangents);

	VNLMatrixType SdotS = kernelObject->Convolve(defCenters);

	TScalar match = targetNonOrientedPolyLine->GetNormSquared();

	for (int i = 0; i < m_NumCells; i++)
	{
		VNLVectorType Mi = special_product(SdotS.get_row(i), defTangents.get_row(i));
		match += dot_product(defTangents.get_row(i), Mi);
	}

	VNLMatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetNonOrientedPolyLine->GetNumberOfCells(); i++)
	{
		VNLVectorType Mi = special_product(SdotT.get_row(i), targTangents.get_row(i));
		match -= 2.0 * dot_product(targTangents.get_row(i), Mi);
	}

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedPolyLine<TScalar, Dimension>::VNLMatrixType
NonOrientedPolyLine<TScalar, Dimension>
::ComputeMatchGradient(NonOrientedPolyLine<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	NonOrientedPolyLine* targetNonOrientedPolyLine = dynamic_cast<NonOrientedPolyLine*>(target);

	if (m_KernelWidth != targetNonOrientedPolyLine->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetNonOrientedPolyLine->GetCenters();
	VNLMatrixType targTangents = targetNonOrientedPolyLine->GetTangents();
	VNLMatrixType targMatrixTangents = targetNonOrientedPolyLine->GetMatrixTangents();

	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];
	VNLMatrixType defCenters, defTangents, defMatrixTangents;
	this->ComputeCentersTangents(defPts, defCenters, defTangents, defMatrixTangents);


	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defMatrixTangents);
	VNLMatrixType KtauS = kernelObject->Convolve(defCenters);
	VNLMatrixList gradKtauS = kernelObject->ConvolveGradient(defCenters);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targMatrixTangents);
	VNLMatrixType KtauT = kernelObject->Convolve(defCenters);
	VNLMatrixList gradKtauT = kernelObject->ConvolveGradient(defCenters);

	VNLMatrixType gradmatch(this->GetNumberOfPoints(),Dimension,0.0);

	for (int f = 0; f < m_NumCells; f++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(f, ptIds);
		m_VTKMutex.Unlock();

		int ind0 = ptIds->GetId(0);
		int ind1 = ptIds->GetId(1);

		VNLVectorType defN = defTangents.get_row(f);
		TScalar defN_mag2 = defN.squared_magnitude();

		if (defN_mag2 > 1e-20)
		{
			VNLVectorType Ktau = special_product(KtauS.get_row(f) - KtauT.get_row(f), defN);
			Ktau *= 4.0;

			VNLVectorType aux = special_product(defMatrixTangents.get_row(f), Ktau);
			Ktau = (Ktau - aux / (2*defN_mag2) ) / sqrt(defN_mag2);

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) - Ktau );
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + Ktau );

			VNLMatrixType delta = gradKtauS[f] - gradKtauT[f];
			VNLVectorType gradKtau(Dimension);
			for (int p = 0; p < Dimension; p++)
			{
				VNLVectorType Mf = special_product(delta.get_column(p), defN);
				gradKtau(p) = dot_product(Mf, defN);
			}

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau);
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau);
		}

	}

	delete kernelObject;

	return gradmatch;
 }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
NonOrientedPolyLine<TScalar, Dimension>
::ComputeCentersTangents(VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Tangents, VNLMatrixType& MatrixTangents)
 {
	Centers.set_size(m_NumCells, Dimension);
	Tangents.set_size(m_NumCells, Dimension);
	MatrixTangents.set_size(m_NumCells, Dimension*(Dimension+1)/2);


	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(i, ptIds);
		m_VTKMutex.Unlock();

		VNLVectorType p0 = Pts.get_row(ptIds->GetId(0));
		VNLVectorType p1 = Pts.get_row(ptIds->GetId(1));
		Centers.set_row(i, (p0 + p1) / 2.0 );

		VNLVectorType Ti = p1-p0;
		Ti /= sqrt( Ti.magnitude() + 1e-20 ); // divided by norm^(1/2)
		Tangents.set_row(i, Ti);

		VNLVectorType aux(Dimension*(Dimension+1)/2);
		int index = 0;
		for (int p = 0; p < Dimension; p++)
			for (int q = p; q < Dimension; q++)
				aux(index++) = Ti(p)*Ti(q);

		MatrixTangents.set_row(i, aux);
	}
 }



template <class TScalar, unsigned int Dimension>
TScalar
NonOrientedPolyLine<TScalar, Dimension>
::ComputeSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixTangents);

	VNLMatrixType selfKW = kernelObject->Convolve(m_Centers);

	TScalar norm2 = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		VNLVectorType Mi = special_product(selfKW.get_row(i), m_Tangents.get_row(i));
		norm2 += dot_product(m_Tangents.get_row(i), Mi);
	}

	delete kernelObject;

	return norm2;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedPolyLine<TScalar, Dimension>::VNLVectorType
NonOrientedPolyLine<TScalar, Dimension>
::special_product(VNLVectorType M, VNLVectorType X)
 {
	VNLMatrixType Ms(Dimension,Dimension);
	int index = 0;
	for (int p = 0; p < Dimension; p++)
		for (int q = p; q < Dimension; q++)
		{
			Ms(p,q) = M(index++);
			Ms(q,p) = Ms(p,q);
		}

	return (Ms * X);
 }


#endif /* _NonOrientedPolyLine_txx */
