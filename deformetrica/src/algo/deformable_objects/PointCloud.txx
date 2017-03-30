/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _PointCloud_txx
#define _PointCloud_txx

#include "PointCloud.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

#include <cstring>
#include <iostream>
#include <sstream>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::PointCloud() : Superclass()
{
	this->SetPointCloudType();
	this->SetModified();
	m_KernelWidth = 0;
}



template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::~PointCloud()
{

}



template <class TScalar, unsigned int Dimension>
PointCloud<TScalar, Dimension>
::PointCloud(const PointCloud& other) : Superclass(other)
 {
	this->SetPointCloudType();

	m_PointWeights = other.m_PointWeights;

	m_KernelWidth = other.m_KernelWidth;
	m_NormSquared = other.m_NormSquared;

 }




////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
PointCloud<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
 {
	Superclass::SetPolyData(polyData);
	this->UpdatePointWeights();
	this->SetModified();

	if (m_KernelWidth != 0)
		m_NormSquared = this->ComputeSelfNorm();

 }



template <class TScalar, unsigned int Dimension>
void
PointCloud<TScalar, Dimension>
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
PointCloud<TScalar, Dimension>
::ComputeMatch(PointCloud<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	PointCloud* targetPointCloud = dynamic_cast<PointCloud*>(target);

	if (m_KernelWidth != targetPointCloud->GetKernelWidth())
		throw std::runtime_error("Kernel width of point clouds mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targPts = targetPointCloud->GetPointCoordinates();
	VNLMatrixType targWts = targetPointCloud->GetPointWeights();
	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defPts);
	kernelObject->SetWeights(m_PointWeights);

	VNLMatrixType SdotS = kernelObject->Convolve(defPts);

	TScalar match = targetPointCloud->GetNormSquared();

	for (int i = 0; i < this->GetNumberOfPoints(); i ++)
		match += SdotS(i,0) * m_PointWeights(i,0);

	VNLMatrixType SdotT = kernelObject->Convolve(targPts);
	for (int i = 0; i < targetPointCloud->GetNumberOfPoints(); i++)
		match -= SdotT(i,0) * (2.0 * targWts(i,0));

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename PointCloud<TScalar, Dimension>::VNLMatrixType
PointCloud<TScalar, Dimension>
::ComputeMatchGradient(PointCloud<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	PointCloud* targetPointCloud = dynamic_cast<PointCloud*>(target);

	if (m_KernelWidth != targetPointCloud->GetKernelWidth())
		throw std::runtime_error("Kernel width of point clouds mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targPts = targetPointCloud->GetPointCoordinates();
	VNLMatrixType targWts = targetPointCloud->GetPointWeights();
	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defPts);
	kernelObject->SetWeights(m_PointWeights);
	VNLMatrixList grad_SdotS = kernelObject->ConvolveGradient(defPts);

	kernelObject->SetSources(targPts);
	kernelObject->SetWeights(targWts);
	VNLMatrixList grad_SdotT = kernelObject->ConvolveGradient(defPts);

	int numPts = this->GetNumberOfPoints();
	VNLMatrixType gradmatch(numPts, Dimension, 0.0);

	for (int i = 0; i < numPts; i ++)
		gradmatch.set_row(i, m_PointWeights(i,0) * (grad_SdotS[i].get_row(0) - grad_SdotT[i].get_row(0)));
	gradmatch *= 2.0;

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
PointCloud<TScalar, Dimension>
::UpdatePointWeights()
 {
	m_PointWeights.set_size(Superclass::m_NumberOfPoints,1);
	m_PointWeights.fill(1.0);

	vtkSmartPointer<vtkPointData> pd = Superclass::m_PolyData->GetPointData();
	int numCmp = pd->GetNumberOfComponents();
	if (numCmp==0)
		std::cout << "Warning: No weights detected: use unit weight for each point" << std::endl;
	else if (numCmp==1)
	{
		vtkSmartPointer<vtkDataArray> scalars = pd->GetScalars();
		int numT = scalars->GetNumberOfTuples();
		if (numT != Superclass::m_NumberOfPoints)
			std::cout << "number of scalars and number of points mismatched! Defaulting to unit weights" << std::endl;
		else
			for (int i = 0; i < numT; i++)
				m_PointWeights.set_row(i, scalars->GetComponent(i,0)); // 0 because of scalar data
	}

 }



template <class TScalar, unsigned int Dimension>
TScalar
PointCloud<TScalar, Dimension>
::ComputeSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(this->GetPointCoordinates());
	kernelObject->SetWeights(m_PointWeights);

	VNLMatrixType selfKW = kernelObject->Convolve(this->GetPointCoordinates());

	TScalar norm2 = 0;

	for (unsigned int i = 0; i < this->GetNumberOfPoints(); i++)
		norm2 += selfKW(i, 0) * m_PointWeights(i,0);

	delete kernelObject;

	return norm2;
 }


#endif /* _PointCloud_txx */
