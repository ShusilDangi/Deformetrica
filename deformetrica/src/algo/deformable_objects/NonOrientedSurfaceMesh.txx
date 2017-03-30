/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _NonOrientedSurfaceMesh_txx
#define _NonOrientedSurfaceMesh_txx

#include "NonOrientedSurfaceMesh.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkTriangleFilter.h"
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
NonOrientedSurfaceMesh<TScalar, Dimension>
::NonOrientedSurfaceMesh() : Superclass()
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetNonOrientedSurfaceMeshType();
	this->SetModified();
	m_KernelWidth = 0;
 }



template <class TScalar, unsigned int Dimension>
NonOrientedSurfaceMesh<TScalar, Dimension>
::~NonOrientedSurfaceMesh()
{

}



template <class TScalar, unsigned int Dimension>
NonOrientedSurfaceMesh<TScalar, Dimension>
::NonOrientedSurfaceMesh(const NonOrientedSurfaceMesh& other) : Superclass(other)
 {
	if (Dimension != 3)
		throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

	this->SetNonOrientedSurfaceMeshType();

	m_Centers = other.m_Centers;
	m_Normals = other.m_Normals;
	m_MatrixNormals = other.m_MatrixNormals;
	m_NumCells = other.m_NumCells;

	m_KernelWidth = other.m_KernelWidth;
	m_NormSquared = other.m_NormSquared;

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
NonOrientedSurfaceMesh<TScalar, Dimension>
::SetPolyData(vtkPolyData* polyData)
 {
	Superclass::SetPolyData(polyData);

	m_NumCells = Superclass::m_PolyData->GetNumberOfCells();

	VNLMatrixType Pts = this->GetPointCoordinates();
	this->ComputeCentersNormals(Pts,m_Centers, m_Normals, m_MatrixNormals);

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
NonOrientedSurfaceMesh<TScalar, Dimension>
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
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatch(NonOrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	NonOrientedSurfaceMesh* targetNonOrientedSurfaceMesh = dynamic_cast<NonOrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetNonOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface unsigned currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetNonOrientedSurfaceMesh->GetCenters();
	VNLMatrixType targNormals = targetNonOrientedSurfaceMesh->GetNormals();

	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];
	VNLMatrixType defCenters, defNormals, defMatrixNormals;
	this->ComputeCentersNormals(defPts, defCenters, defNormals, defMatrixNormals);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defMatrixNormals);

	VNLMatrixType SdotS = kernelObject->Convolve(defCenters);

	TScalar match = targetNonOrientedSurfaceMesh->GetNormSquared();

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		VNLVectorType Mi = special_product(SdotS.get_row(i), defNormals.get_row(i));
		match += dot_product(defNormals.get_row(i), Mi);		
	}

	VNLMatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetNonOrientedSurfaceMesh->GetNumberOfCells(); i++)
	{
		VNLVectorType Mi = special_product(SdotT.get_row(i), targNormals.get_row(i));				
		match -= 2.0 * dot_product(targNormals.get_row(i), Mi);		
	}

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
TScalar
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatch(NonOrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable object types mismatched");

	NonOrientedSurfaceMesh* targetNonOrientedSurfaceMesh = dynamic_cast<NonOrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetNonOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface unsigned currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();

	// Make sure t is inside the time interval
	if (t >= numTimePts)
	{
		std::cerr << "time for matching is outside of time interval" << std::endl;
	}

	VNLMatrixType targCenters = targetNonOrientedSurfaceMesh->GetCenters();
	VNLMatrixType targNormals = targetNonOrientedSurfaceMesh->GetNormals();

	VNLMatrixType defPts = Superclass::m_PointsT[t];
	VNLMatrixType defCenters, defNormals, defMatrixNormals;
	this->ComputeCentersNormals(defPts, defCenters, defNormals, defMatrixNormals);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defMatrixNormals);

	VNLMatrixType SdotS = kernelObject->Convolve(defCenters);

	TScalar match = targetNonOrientedSurfaceMesh->GetNormSquared();

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		VNLVectorType Mi = special_product(SdotS.get_row(i), defNormals.get_row(i));
		match += dot_product(defNormals.get_row(i), Mi);		
	}

	VNLMatrixType SdotT = kernelObject->Convolve(targCenters);
	for (int i = 0; i < targetNonOrientedSurfaceMesh->GetNumberOfCells(); i++)
	{
		VNLVectorType Mi = special_product(SdotT.get_row(i), targNormals.get_row(i));				
		match -= 2.0 * dot_product(targNormals.get_row(i), Mi);		
	}

	delete kernelObject;

	return match;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedSurfaceMesh<TScalar, Dimension>::VNLMatrixType
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatchGradient(NonOrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	NonOrientedSurfaceMesh* targetNonOrientedSurfaceMesh = dynamic_cast<NonOrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetNonOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();
	VNLMatrixType targCenters = targetNonOrientedSurfaceMesh->GetCenters();
	VNLMatrixType targNormals = targetNonOrientedSurfaceMesh->GetNormals();
	VNLMatrixType targMatrixNormals = targetNonOrientedSurfaceMesh->GetMatrixNormals();

	VNLMatrixType defPts = Superclass::m_PointsT[numTimePts-1];
	VNLMatrixType defCenters, defNormals, defMatrixNormals;
	this->ComputeCentersNormals(defPts, defCenters, defNormals, defMatrixNormals);


	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defMatrixNormals);
	VNLMatrixType KtauS = kernelObject->Convolve(defCenters);
	VNLMatrixList gradKtauS = kernelObject->ConvolveGradient(defCenters);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targMatrixNormals);
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
		int ind2 = ptIds->GetId(2);

		VNLVectorType e0 = defPts.get_row(ind2) - defPts.get_row(ind1);
		VNLVectorType e1 = defPts.get_row(ind0) - defPts.get_row(ind2);
		VNLVectorType e2 = defPts.get_row(ind1) - defPts.get_row(ind0);

		VNLVectorType defN = defNormals.get_row(f);
		TScalar defN_mag2 = defN.squared_magnitude();

		if (defN_mag2 > 1e-20)
		{
			VNLVectorType Ktau = special_product(KtauS.get_row(f) - KtauT.get_row(f), defN);
			Ktau *= 2.0;

			VNLVectorType aux = special_product(defMatrixNormals.get_row(f), Ktau);
			Ktau = (Ktau - aux / (2*defN_mag2) ) / sqrt(defN_mag2);

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + vnl_cross_3d(e0, Ktau) );
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + vnl_cross_3d(e1, Ktau) );
			gradmatch.set_row(ind2, gradmatch.get_row(ind2) + vnl_cross_3d(e2, Ktau) );

			VNLMatrixType delta = gradKtauS[f] - gradKtauT[f];
			VNLVectorType gradKtau(Dimension);
			for (int p = 0; p < Dimension; p++)
			{
				VNLVectorType Mf = special_product(delta.get_column(p), defN);
				gradKtau(p) = dot_product(Mf, defN);
			}
			gradKtau *= 2.0/3.0;

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau);
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau);
			gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau);		
		}

	}

	delete kernelObject;

	return gradmatch;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedSurfaceMesh<TScalar, Dimension>::VNLMatrixType
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeMatchGradient(NonOrientedSurfaceMesh<TScalar, Dimension>::DeformableObjectType* target, unsigned int t)
 {
	if (this->GetType() != target->GetType())
		throw std::runtime_error("Deformable objects types mismatched");

	NonOrientedSurfaceMesh* targetNonOrientedSurfaceMesh = dynamic_cast<NonOrientedSurfaceMesh*>(target);

	if (m_KernelWidth != targetNonOrientedSurfaceMesh->GetKernelWidth())
		throw std::runtime_error("Kernel width of surface currents mismatched");

	unsigned int numTimePts = DeformableObjectType::m_Def->GetNumberOfTimePoints();

	// Make sure t is inside the time interval
	if (t >= numTimePts)
	{
		std::cerr << "time for matching is outside of time interval" << std::endl;
	}

	VNLMatrixType targCenters = targetNonOrientedSurfaceMesh->GetCenters();
	VNLMatrixType targNormals = targetNonOrientedSurfaceMesh->GetNormals();
	VNLMatrixType targMatrixNormals = targetNonOrientedSurfaceMesh->GetMatrixNormals();

	VNLMatrixType defPts = Superclass::m_PointsT[t];
	VNLMatrixType defCenters, defNormals, defMatrixNormals;
	this->ComputeCentersNormals(defPts, defCenters, defNormals, defMatrixNormals);


	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);

	kernelObject->SetSources(defCenters);
	kernelObject->SetWeights(defMatrixNormals);
	VNLMatrixType KtauS = kernelObject->Convolve(defCenters);
	VNLMatrixList gradKtauS = kernelObject->ConvolveGradient(defCenters);

	kernelObject->SetSources(targCenters);
	kernelObject->SetWeights(targMatrixNormals);
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
		int ind2 = ptIds->GetId(2);

		VNLVectorType e0 = defPts.get_row(ind2) - defPts.get_row(ind1);
		VNLVectorType e1 = defPts.get_row(ind0) - defPts.get_row(ind2);
		VNLVectorType e2 = defPts.get_row(ind1) - defPts.get_row(ind0);

		VNLVectorType defN = defNormals.get_row(f);
		TScalar defN_mag2 = defN.squared_magnitude();

		if (defN_mag2 > 1e-20)
		{
			VNLVectorType Ktau = special_product(KtauS.get_row(f) - KtauT.get_row(f), defN);
			Ktau *= 2.0;

			VNLVectorType aux = special_product(defMatrixNormals.get_row(f), Ktau);
			Ktau = (Ktau - aux / (2*defN_mag2) ) / sqrt(defN_mag2);

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + vnl_cross_3d(e0, Ktau) );
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + vnl_cross_3d(e1, Ktau) );
			gradmatch.set_row(ind2, gradmatch.get_row(ind2) + vnl_cross_3d(e2, Ktau) );

			VNLMatrixType delta = gradKtauS[f] - gradKtauT[f];
			VNLVectorType gradKtau(Dimension);
			for (int p = 0; p < Dimension; p++)
			{
				VNLVectorType Mf = special_product(delta.get_column(p), defN);
				gradKtau(p) = dot_product(Mf, defN);
			}
			gradKtau *= 2.0/3.0;

			gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau);
			gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau);
			gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau);		
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
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeCentersNormals(VNLMatrixType& Pts, VNLMatrixType& Centers, VNLMatrixType& Normals, VNLMatrixType& MatrixNormals)
 {
	Centers.set_size(m_NumCells, Dimension);
	Normals.set_size(m_NumCells, Dimension);
	MatrixNormals.set_size(m_NumCells, Dimension*(Dimension+1)/2);

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

		VNLVectorType Ni = vnl_cross_3d(p2 - p0, p1 - p0) / 2;
		Ni /= sqrt( Ni.magnitude() + 1e-20 ); // divided by norm^(1/2)
		Normals.set_row(i, Ni);

		VNLVectorType aux(Dimension*(Dimension+1)/2);
		int index = 0;
		for (int p = 0; p < Dimension; p++)
			for (int q = p; q < Dimension; q++)
				aux(index++) = Ni(p)*Ni(q);

		MatrixNormals.set_row(i, aux);
	}
 }



template <class TScalar, unsigned int Dimension>
TScalar
NonOrientedSurfaceMesh<TScalar, Dimension>
::ComputeSelfNorm()
 {
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObject = kFactory->CreateKernelObject();
	kernelObject->SetKernelWidth(m_KernelWidth);
	kernelObject->SetSources(m_Centers);
	kernelObject->SetWeights(m_MatrixNormals);

	VNLMatrixType selfKW = kernelObject->Convolve(m_Centers);

	TScalar norm2 = 0;

	for (unsigned int i = 0; i < m_NumCells; i++)
	{
		VNLVectorType Mi = special_product(selfKW.get_row(i), m_Normals.get_row(i));
		norm2 += dot_product(m_Normals.get_row(i), Mi);
	}

	delete kernelObject;

	return norm2;
 }



template <class TScalar, unsigned int Dimension>
typename NonOrientedSurfaceMesh<TScalar, Dimension>::VNLVectorType
NonOrientedSurfaceMesh<TScalar, Dimension>
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



#endif
