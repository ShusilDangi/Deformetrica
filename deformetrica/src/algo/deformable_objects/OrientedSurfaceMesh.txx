/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _OrientedSurfaceMesh_txx
#define _OrientedSurfaceMesh_txx

#include "OrientedSurfaceMesh.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkTriangleFilter.h"
#include "myvtkPolyDataNormals.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include "vnl/vnl_cross.h"
#include "vnl/vnl_matrix.h"

#include <cstring>
#include <iostream>
#include <sstream>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
OrientedSurfaceMesh<TScalar, Dimension>
::OrientedSurfaceMesh() : Superclass()
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetOrientedSurfaceMeshType();
	this->SetModified();
	m_KernelWidth = 0;
 }


template <class TScalar, unsigned int Dimension>
OrientedSurfaceMesh<TScalar, Dimension>
::~OrientedSurfaceMesh()
{

}

template <class TScalar, unsigned int Dimension>
OrientedSurfaceMesh<TScalar, Dimension>
::OrientedSurfaceMesh(const OrientedSurfaceMesh& other) : Superclass(other)
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetOrientedSurfaceMeshType();

	m_Centers = other.m_Centers;
	m_Normals = other.m_Normals;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_NormSquared = other.m_NormSquared;

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedSurfaceMesh<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData, bool reorient)
 {
	Superclass::SetPolyData(polyData);

	this->CheckMeshAndNormals(reorient);

	VNLMatrixType Pts = this->GetPointCoordinates();
	this->ComputeCentersNormals(Pts,m_Centers,m_Normals);

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
OrientedSurfaceMesh<TScalar, Dimension>
::SetKernelWidth(TScalar h)
 {
	m_KernelWidth = h;

	if (Superclass::m_NumberOfPoints != 0)
		m_NormSquared = this->ComputeSelfNorm();
 }




////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int Dimension>
TScalar OrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatch(OrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target)
 {
	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	return this->ComputeMatch(target, numTimePts-1);
 }



template<class TScalar, unsigned int Dimension>
TScalar OrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatch(OrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	OrientedSurfaceMesh* targetOrientedSurfaceMesh = dynamic_cast<OrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();

	// Make sure t is inside the time interval
	if (t >= numTimePts)
	{
		std::cerr << "time for matching is outside of time interval" << std::endl;
	}

	VNLMatrixType targCenters = targetOrientedSurfaceMesh->GetCenters();
	VNLMatrixType targNormals = targetOrientedSurfaceMesh->GetNormals();

	// Get the points at the time point of interest
	VNLMatrixType defPts = Superclass::m_PointsT[t];
	VNLMatrixType defCenters, defNormals;
	this->ComputeCentersNormals(defPts, defCenters, defNormals);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defNormals);

	VNLMatrixType SdotS = kernelObject->Convolve(defCenters);

	TScalar match = targetOrientedSurfaceMesh->GetNormSquared();

	for (int i = 0; i < m_NumCells; i++)
		match += dot_product(SdotS.get_row(i), defNormals.get_row(i));

	VNLMatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetOrientedSurfaceMesh->GetNumberOfCells(); i++)
		match -= 2.0 * dot_product(SdotT.get_row(i), targNormals.get_row(i));

	delete kernelObject;

	return match;
 }



template<class TScalar, unsigned int Dimension>
typename OrientedSurfaceMesh<TScalar, Dimension>::VNLMatrixType OrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatchGradient(OrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target)
 {
	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	return this->ComputeMatchGradient(target, numTimePts-1);
 }



template<class TScalar, unsigned int Dimension>
typename OrientedSurfaceMesh<TScalar, Dimension>::VNLMatrixType OrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatchGradient(OrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	OrientedSurfaceMesh* targetOrientedSurfaceMesh = dynamic_cast<OrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();

	// Make sure t is inside the time interval
	if (t >= numTimePts)
	{
		std::cerr << "time for matching is outside of time interval" << std::endl;
	}

	VNLMatrixType targCenters = targetOrientedSurfaceMesh->GetCenters();
	VNLMatrixType targNormals = targetOrientedSurfaceMesh->GetNormals();

	VNLMatrixType defPts = Superclass::m_PointsT[t];
	VNLMatrixType defCenters, defNormals;
	this->ComputeCentersNormals(defPts, defCenters, defNormals);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defNormals);
	VNLMatrixType KtauS = kernelObject->Convolve(defCenters);
	VNLMatrixType gradKtauS = kernelObject->ConvolveGradient(defCenters, defNormals);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targNormals);
	VNLMatrixType KtauT = kernelObject->Convolve(defCenters);
	VNLMatrixType gradKtauT = kernelObject->ConvolveGradient(defCenters, defNormals);

	VNLMatrixType gradKtau = (gradKtauS - gradKtauT) * 2.0 / 3.0;

	VNLMatrixType gradmatch(this->GetNumberOfPoints(), Dimension, 0.0);
	for (int f = 0; f < m_NumCells; f++)
	{
		vtkSmartPointer < vtkIdList > ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(f, ptIds);
		m_VTKMutex.Unlock();

		int ind0 = ptIds->GetId(0);
		int ind1 = ptIds->GetId(1);
		int ind2 = ptIds->GetId(2);

		VNLVectorType e0 = defPts.get_row(ind2) - defPts.get_row(ind1);
		VNLVectorType e1 = defPts.get_row(ind0) - defPts.get_row(ind2);
		VNLVectorType e2 = defPts.get_row(ind1) - defPts.get_row(ind0);

		VNLVectorType Ktau = KtauS.get_row(f) - KtauT.get_row(f);
		gradmatch.set_row(ind0, gradmatch.get_row(ind0) + vnl_cross_3d(e0, Ktau));
		gradmatch.set_row(ind1, gradmatch.get_row(ind1) + vnl_cross_3d(e1, Ktau));
		gradmatch.set_row(ind2, gradmatch.get_row(ind2) + vnl_cross_3d(e2, Ktau));

		gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau.get_row(f));
		gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau.get_row(f));
		gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau.get_row(f));

	}

	delete kernelObject;

	return gradmatch;
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
OrientedSurfaceMesh<TScalar, Dimension>
::CheckMeshAndNormals(bool reorient)
 {

	vtkSmartPointer<vtkTriangleFilter> trif =
			vtkSmartPointer<vtkTriangleFilter>::New();
	
#if (VTK_MAJOR_VERSION == 5)
	trif->SetInput(Superclass::m_PolyData);
#else
	trif->SetInputData(Superclass::m_PolyData);
#endif
	
	trif->PassVertsOff();
	trif->PassLinesOff();
	trif->Update();

	vtkSmartPointer<myvtkPolyDataNormals> normalf =
			vtkSmartPointer<myvtkPolyDataNormals>::New();

#if (VTK_MAJOR_VERSION == 5)
	normalf->SetInput(trif->GetOutput());
#else
	normalf->SetInputData(trif->GetOutput());
#endif

	normalf->ComputePointNormalsOff();
	normalf->ComputeCellNormalsOn();
	normalf->SplittingOff();

	if (reorient)
	{
		normalf->ConsistencyOn();
		normalf->AutoOrientNormalsOn(); // Should have closed surface
		normalf->FlipNormalsOn();
	}
	else
	{
		normalf->ConsistencyOff();
		normalf->AutoOrientNormalsOff();
	}
	normalf->Update();

	Superclass::m_PolyData = normalf->GetOutput();
	Superclass::m_PolyData->BuildLinks();

	m_NumCells = Superclass::m_PolyData->GetNumberOfCells();

 }



template <class TScalar, unsigned int Dimension>
void
OrientedSurfaceMesh<TScalar, Dimension>
::ComputeCentersNormals(const VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Normals)
 {
	Centers.set_size(m_NumCells, Dimension);
	Normals.set_size(m_NumCells, Dimension);

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

		m_VTKMutex.Lock();
		Superclass::m_PolyData->GetCellPoints(i, ptIds);
		m_VTKMutex.Unlock();

		if (ptIds->GetNumberOfIds() != 3)
			throw std::runtime_error("Not a triangle cell!");

		VNLVectorType p0 = Pts.get_row(ptIds->GetId(0));
		VNLVectorType p1 = Pts.get_row(ptIds->GetId(1));
		VNLVectorType p2 = Pts.get_row(ptIds->GetId(2));

		Centers.set_row(i, (p0 + p1 + p2) / 3.0 );
		Normals.set_row(i, vnl_cross_3d(p2 - p0, p1 - p0) / 2);

	}
 }



template <class TScalar, unsigned int Dimension>
TScalar
OrientedSurfaceMesh<TScalar, Dimension>
::ComputeSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_Normals);

	VNLMatrixType selfKW = kernelObject->Convolve(m_Centers);

	TScalar norm2 = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
		norm2 += dot_product(selfKW.get_row(i), m_Normals.get_row(i));

	delete kernelObject;

	return norm2;
 }


#endif /* _OrientedSurfaceMesh_txx */
