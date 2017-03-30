/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _deform_txx
#define _deform_txx

#include "KernelFactory.h"
#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "DeformableObject.h"

#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"

#include "DeformableObjectReader.h"

#if ITK_VERSION_MAJOR >= 4
#include "itkFFTWGlobalConfiguration.h"
#endif

#include "itksys/SystemTools.hxx"

#include <cstring>
#include <iostream>
#include <sstream>



template <unsigned int Dimension>
void deform(SparseDiffeoParameters* paramDiffeos,
		bool useReverseFlow,
		vnl_matrix<float>& CP0,
		vnl_matrix<float>& MOM0,
		int numObjects,
		std::vector<DeformableObjectParameters::Pointer> paramObjectsList,
		const std::vector<char*> objectfnList)
{
	typedef float TScalar;

	std::cout << "Shoot and Flow\n===" << std::endl;
	std::cout << "Number of objects to flow: " << numObjects << std::endl;
	std::cout << "\n===\n" << std::endl;
	std::cout << "Deformations Parameters:" << std::endl;
	paramDiffeos->PrintSelf(std::cout);
	std::cout << "\n===\n" << std::endl;
	for (int i = 0; i < numObjects; i++)
	{

		std::cout << "Shoot and Flow\n===" << std::endl;
		std::cout << "Number of objects to flow: " << numObjects << std::endl;
		std::cout << "\n===\n" << std::endl;
		std::cout << "Deformations Parameters:" << std::endl;
		paramDiffeos->PrintSelf(std::cout);
		std::cout << "\n===\n" << std::endl;
		for (int i = 0; i < numObjects; i++)
		{
			std::cout << "Object " << i << ": " << std::endl;
			std::cout << "Source file: " << objectfnList[i] << std::endl;
			paramObjectsList[i]->PrintSelf(std::cout);
			std::cout << "\n===\n" << std::endl;
		}
		std::cout << "number of CPs: " << CP0.rows() << std::endl;

		if (CP0.rows() != MOM0.rows())
			throw std::runtime_error("Number of CPs and Momentas mismatched");

		vnl_vector<TScalar> minData(Dimension, 1e20);
		vnl_vector<TScalar> maxData(Dimension, -1e20);

		for (int k = 0; k < CP0.rows(); k++)
		{
			for (int d = 0; d < Dimension; d++)
			{
				minData[d] = (CP0(k, d) < minData[d]) ? CP0(k, d) : minData[d];
				maxData[d] = (CP0(k, d) > maxData[d]) ? CP0(k, d) : maxData[d];
			}
		}
#if ITK_VERSION_MAJOR >= 4
		itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

		typedef Diffeos<TScalar, Dimension> DiffeosType;
		typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
		typedef std::vector<DeformableObjectType*> DeformableObjectList;
		typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;

		// create the deformation object
		DiffeosType* def = new DiffeosType();
		def->SetKernelWidth(paramDiffeos->GetKernelWidth());
		def->SetNumberOfTimePoints(paramDiffeos->GetNumberOfTimePoints());
		def->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor()); // to define the bounding box
		def->UseImprovedEuler();

		typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
		KernelFactoryType* kfac = KernelFactoryType::Instantiate();
		if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "fgt") == 0)
		{
			kfac->UseFGTKernel();
		}
		else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0)
		{
			kfac->UseP3MKernel();
		}

#ifdef USE_CUDA
		else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") == 0)
		{
			kfac->UseExactKernel();
		}
		else
		{
			if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "cudaexact") != 0)
				std::cerr << "Unknown kernel type: defaulting to CudaExact" << std::endl;
			kfac->UseCUDAExactKernel();
		}
#else
		else
		{
			if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "exact") != 0)
				std::cerr << "Unknown kernel type: defaulting to Exact" << std::endl;
			kfac->UseExactKernel();
		}
#endif
		kfac->SetWorkingSpacingRatio(paramDiffeos->GetP3MWorkingSpacingRatio());
		kfac->SetPaddingFactor(paramDiffeos->GetP3MPaddingFactor());


		// create source object
		vnl_matrix<TScalar> bbox;
		DeformableObjectList objectList(numObjects);
		for (int i = 0; i < numObjects; i++)
		{
			typedef DeformableObjectReader<TScalar, Dimension> ObjectReaderType;
			ObjectReaderType* reader = new ObjectReaderType();
			reader->SetObjectParameters(paramObjectsList[i]);

			reader->SetFileName(objectfnList[i]);
			reader->SetTemplateType();
			reader->Update();

			objectList[i] = reader->GetOutput();
			bbox = reader->GetBoundingBox();
			for (unsigned int d = 0; d < Dimension; d++)
			{
				minData[d] = (bbox(d, 0) < minData[d]) ? bbox(d, 0) : minData[d];
				maxData[d] = (bbox(d, 1) > maxData[d]) ? bbox(d, 1) : maxData[d];
			}
		}
		DeformableMultiObjectType* object = new DeformableMultiObjectType();
		object->SetObjectList(objectList);

		// set bounding box
		vnl_matrix<TScalar> DataDomain(Dimension, 2);
		DataDomain.set_column(0, minData);
		DataDomain.set_column(1, maxData);
		std::cout << "Working domain: origin =  [" << minData << "] length = [" << maxData - minData << "]" << std::endl;

		def->SetDataDomain(DataDomain);
		kfac->SetDataDomain(DataDomain);

		// shoot and flow
		std::cout << "Shoot..." << std::endl;
		def->SetStartPositions(CP0);
		def->SetStartMomentas(MOM0);
		def->Update();

		if (def->OutOfBox())
			throw std::runtime_error("out of box");

		if (useReverseFlow)
		{
			def->ReverseFlow();
		}


		std::cout << "... and Flow" << std::endl;
		object->SetDiffeos(def);
		object->Update();

		//if (object->OutOfBox())
		//	throw std::runtime_error("out of box");

		std::cout << "Write output files" << std::endl;

		std::vector<std::string> outname(numObjects);
		std::vector<std::string> outext(numObjects);
		std::vector<std::string> outfn(numObjects);
		for (int i = 0; i < numObjects; i++)
		{
			std::string sourcefn;
			sourcefn.assign(objectfnList[i]);
			int length = sourcefn.size();
			int index = length - 1;
			while ((sourcefn[index] != '.') && (index >= 0))
				index--;

			outname[i] = sourcefn.substr(0, index);
			outext[i] = sourcefn.substr(index, length - index);

			std::ostringstream oss;
			oss << outname[i] << "_flow";
			outname[i] = oss.str();
		}

		object->WriteFlow(outname, outext);

		delete object;

	}

#endif

}
