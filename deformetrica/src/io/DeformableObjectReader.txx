/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectReader_txx
#define _DeformableObjectReader_txx

#include "Landmark.h"
#include "PointCloud.h"
#include "OrientedPolyLine.h"
#include "NonOrientedPolyLine.h"
#include "OrientedSurfaceMesh.h"
#include "NonOrientedSurfaceMesh.h"

#include "KernelFactory.h"

#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "itksys/SystemTools.hxx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeformableObjectReader<TScalar, Dimension>
::DeformableObjectReader()
 {
	m_Bbox.set_size(Dimension,2);
	m_Min.set_size(Dimension);
	m_Max.set_size(Dimension);
	m_IsTemplate = false;
 }



template <class TScalar, unsigned int Dimension>
DeformableObjectReader<TScalar, Dimension>
::~DeformableObjectReader()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
DeformableObjectReader<TScalar, Dimension>
::Update()
 {

	typedef Landmark<TScalar, Dimension> LandmarkType;
	typedef PointCloud<TScalar, Dimension> PointCloudType;
	typedef OrientedPolyLine<TScalar, Dimension> OrientedPolyLineType;
	typedef NonOrientedPolyLine<TScalar, Dimension> NonOrientedPolyLineType;
	typedef OrientedSurfaceMesh<TScalar, Dimension> OrientedSurfaceMeshType;
	typedef NonOrientedSurfaceMesh<TScalar, Dimension> NonOrientedSurfaceMeshType;

	const char* ObjectType = m_ParamObject->GetDeformableObjectType().c_str();

	if (itksys::SystemTools::Strucmp(ObjectType,"Landmark") == 0)
	{
		LandmarkType* objectLandmark = new LandmarkType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectLandmark->SetPolyData(PolyData);

			m_Min = objectLandmark->GetMin();
			m_Max = objectLandmark->GetMax();
		}
		m_Object = objectLandmark;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"PointCloud") == 0)
	{
		PointCloudType* objectPointCloud = new PointCloudType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectPointCloud->SetPolyData(PolyData);

			m_Min = objectPointCloud->GetMin();
			m_Max = objectPointCloud->GetMax();

			m_Bbox.set_column(0,m_Min);
			m_Bbox.set_column(1,m_Max);
			typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
			KernelFactoryType* kfac = KernelFactoryType::Instantiate();
			kfac->SetDataDomain(m_Bbox);

			objectPointCloud->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectPointCloud;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"OrientedPolyLine") == 0)
	{
		OrientedPolyLineType* objectOrientedPolyLine = new OrientedPolyLineType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectOrientedPolyLine->SetPolyData(PolyData);

			m_Min = objectOrientedPolyLine->GetMin();
			m_Max = objectOrientedPolyLine->GetMax();

			m_Bbox.set_column(0,m_Min);
			m_Bbox.set_column(1,m_Max);
			typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
			KernelFactoryType* kfac = KernelFactoryType::Instantiate();
			kfac->SetDataDomain(m_Bbox);

			objectOrientedPolyLine->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectOrientedPolyLine;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"NonOrientedPolyLine") == 0)
	{
		NonOrientedPolyLineType* objectNonOrientedPolyLine = new NonOrientedPolyLineType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectNonOrientedPolyLine->SetPolyData(PolyData);

			m_Min = objectNonOrientedPolyLine->GetMin();
			m_Max = objectNonOrientedPolyLine->GetMax();

			m_Bbox.set_column(0,m_Min);
			m_Bbox.set_column(1,m_Max);
			typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
			KernelFactoryType* kfac = KernelFactoryType::Instantiate();
			kfac->SetDataDomain(m_Bbox);

			objectNonOrientedPolyLine->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectNonOrientedPolyLine;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"OrientedSurfaceMesh") == 0)
	{
		bool ReOrientSurface = m_ParamObject->ReOrient();

		OrientedSurfaceMeshType* objectOrientedSurfaceMesh = new OrientedSurfaceMeshType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectOrientedSurfaceMesh->SetPolyData(PolyData, ReOrientSurface);

			m_Min = objectOrientedSurfaceMesh->GetMin();
			m_Max = objectOrientedSurfaceMesh->GetMax();

			m_Bbox.set_column(0,m_Min);
			m_Bbox.set_column(1,m_Max);
			typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
			KernelFactoryType* kfac = KernelFactoryType::Instantiate();
			kfac->SetDataDomain(m_Bbox);

			objectOrientedSurfaceMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectOrientedSurfaceMesh;
	}

	else if (itksys::SystemTools::Strucmp(ObjectType,"NonOrientedSurfaceMesh") == 0)
	{
		NonOrientedSurfaceMeshType* objectNonOrientedSurfaceMesh = new NonOrientedSurfaceMeshType();
		{
			vtkSmartPointer<vtkPolyDataReader> reader =
					vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(m_FileName);
			reader->Update();
			vtkSmartPointer<vtkPolyData> PolyData = reader->GetOutput();

			objectNonOrientedSurfaceMesh->SetPolyData(PolyData);

			m_Min = objectNonOrientedSurfaceMesh->GetMin();
			m_Max = objectNonOrientedSurfaceMesh->GetMax();

			m_Bbox.set_column(0,m_Min);
			m_Bbox.set_column(1,m_Max);
			typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
			KernelFactoryType* kfac = KernelFactoryType::Instantiate();
			kfac->SetDataDomain(m_Bbox);

			objectNonOrientedSurfaceMesh->SetKernelWidth(m_ParamObject->GetKernelWidth());
		}
		m_Object = objectNonOrientedSurfaceMesh;
	}

	else
	{
		std::cerr << "Unknown object type: " << m_ParamObject->GetDeformableObjectType() <<
				"\nCurrently available types are: Landmark, PointCloud, OrientedPolyLine, "
				"NonOrientedPolyLine, OrientedSurfaceMesh, NonOrientedSurfaceMesh" << std::endl;
		return;
	}

	if (m_IsTemplate)
	{
		m_Object->SetDataSigma(m_ParamObject->GetDataSigma());
	}

	m_Bbox.set_column(0, m_Min);
	m_Bbox.set_column(1, m_Max);

 }

#endif
