/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _estimateAtlas_txx
#define _estimateAtlas_txx

#include "DeformetricaConfig.h"

#include "itkVersion.h"
#include "vtkVersion.h"

#include "KernelFactory.h"
#include "SparseDiffeoAtlasEstimator.h"
#include "Diffeos.h"

#include "DeformableMultiObject.h"
#include "DeformableObject.h"

#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"

#include "DeformableObjectReader.h"

#include "SimpleTimer.h"

#if ITK_VERSION_MAJOR >= 4
#include "itkFFTWGlobalConfiguration.h"
#endif

#include "itksys/SystemTools.hxx"

#include "vnl/vnl_vector.h"

#include <cstring>
#include <iostream>
#include <sstream>

template <unsigned int Dimension>
void estimateAtlas(
		SparseDiffeoParameters* paramDiffeos,
		int numSubjects,
		int numObjects,
		std::vector<DeformableObjectParameters::Pointer>& paramObjectsList,
		const std::vector<char*>& templatefnList,
		const std::vector< std::vector<char*> >& subjectfnList)
{
#ifdef USE_DOUBLE_PRECISION
	typedef double TScalar;
	std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
	typedef float TScalar;
	std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif

	std::cout << "Sparse diffeomorphic atlas estimation\n===" << std::endl;
	std::cout << "ITK version " << itk::Version::GetITKVersion() << std::endl;
	std::cout << "VTK version " << vtkVersion::GetVTKVersion() << std::endl;
	std::cout << "Deformation Parameters:" << std::endl;
	paramDiffeos->PrintSelf(std::cout);
	std::cout << "\n===\n" << std::endl;
	for (int i = 0; i < numObjects; i++)
	{
		std::cout << "Object " << i << ": " << std::endl;
		std::cout << "template file: " << templatefnList[i] << std::endl;
		paramObjectsList[i]->PrintSelf(std::cout);
		std::cout << "\n===\n" << std::endl;
	}


#if ITK_VERSION_MAJOR >= 4
	itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

	typedef Diffeos<TScalar, Dimension> DiffeosType;

	vnl_vector<TScalar> minData(Dimension,1e20);
	vnl_vector<TScalar> maxData(Dimension,-1e20);

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


	// create template and target objects
	DeformableObjectList templateObjectList(numObjects);
	std::vector<DeformableObjectList> targetObjectList(numSubjects);
	for (int s = 0; s < numSubjects; s++)
	{
		DeformableObjectList Aux(numObjects);
		targetObjectList[s] = Aux;
	}

	// scan the list of objects
	typedef DeformableObjectReader<TScalar,Dimension> ObjectReaderType;
	vnl_matrix<TScalar> bbox;
	for (int i = 0; i < numObjects; i++)
	{
		for (int s = 0; s < numSubjects; s++)
		{
			ObjectReaderType* reader = new ObjectReaderType();
			reader->SetObjectParameters(paramObjectsList[i]);

			reader->SetFileName(subjectfnList[i][s]);
			reader->Update();

			targetObjectList[s][i] = reader->GetOutput();
			bbox = reader->GetBoundingBox();
			for (unsigned int d=0; d<Dimension; d++)
			{
				minData[d] = (bbox(d,0) < minData[d])?bbox(d,0):minData[d];
				maxData[d] = (bbox(d,1) > maxData[d])?bbox(d,1):maxData[d];
			}
		}

		ObjectReaderType* reader = new ObjectReaderType();
		reader->SetObjectParameters(paramObjectsList[i]);

		reader->SetFileName(templatefnList[i]);
		reader->SetTemplateType();
		reader->Update();

		templateObjectList[i] = reader->GetOutput();
		bbox = reader->GetBoundingBox();
		for (unsigned int d=0; d<Dimension; d++)
		{
			minData[d] = (bbox(d,0) < minData[d])?bbox(d,0):minData[d];
			maxData[d] = (bbox(d,1) > maxData[d])?bbox(d,1):maxData[d];
		}
	}


	DeformableMultiObjectType* templateObj = new DeformableMultiObjectType();
	templateObj->SetObjectList(templateObjectList);

	std::vector< DeformableMultiObjectType* > target(numSubjects);
	for (unsigned int s = 0; s < numSubjects; s++)
	{
		target[s] = new DeformableMultiObjectType();
		target[s]->SetObjectList(targetObjectList[s]);
	}


	vnl_matrix<TScalar> DataDomain(Dimension,2);
	DataDomain.set_column(0,minData);
	DataDomain.set_column(1,maxData);

	std::cout << "Working domain: origin =  [" << minData << "] length = [" << maxData - minData << "]" << std::endl; 

	// Set up DataDomain for kernel factory and deformation
	def->SetDataDomain(DataDomain);
	kfac->SetDataDomain(DataDomain);

	double InitialCPSpacing = paramDiffeos->GetInitialCPSpacing();
	if (InitialCPSpacing < 1e-20)
	{
		InitialCPSpacing = paramDiffeos->GetKernelWidth();
		std::cout << "InitialCPSpacing set to " << InitialCPSpacing << std::endl;
	}

	// create timer
	SimpleTimer timer;

	// create the AtlasEstimator!
	typedef SparseDiffeoAtlasEstimator<TScalar, Dimension> AtlasEstimator;
	AtlasEstimator* atlasbuilder = new AtlasEstimator();
	atlasbuilder->SetTargetList(target);
	atlasbuilder->SetTemplate(templateObj);
	atlasbuilder->SetDiffeos(def);
	atlasbuilder->SetInitialCPSpacing(InitialCPSpacing);
	paramDiffeos->FreezeCP()?atlasbuilder->SetFreezeCP():atlasbuilder->UnsetFreezeCP();
	atlasbuilder->SetInitialCPPosition(paramDiffeos->GetInitialCPPosition_fn()); // null string means no InitialCPPosition given (will be initialized with a regular lattice)
	atlasbuilder->SetInitialMomentas(paramDiffeos->GetInitialMomenta_fn()); // null string means no InitialMomenta given (will be initialized with zero)
	atlasbuilder->SetSparsityPrior(paramDiffeos->GetSparsityPrior());
	atlasbuilder->SetRegularityWeights(paramDiffeos->GetRegularityWeights());
	if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientdescent") == 0)
	{
		atlasbuilder->SetGradientDescent();
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "ista") == 0)
	{
		atlasbuilder->SetISTA();
	}
	else
	{
		if ( (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "f_ista") != 0) &&
				(itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fista") != 0) )
			std::cerr << "Unknown optimization method type: defaulting to FISTA" << std::endl;
			atlasbuilder->SetFISTA();
	}
	atlasbuilder->SetSmoothingKernelWidth(paramDiffeos->GetSmoothingKernelWidthRatio() * paramDiffeos->GetKernelWidth());
	atlasbuilder->SetMaxIterations(paramDiffeos->GetMaxIterations());
	atlasbuilder->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
	atlasbuilder->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
	atlasbuilder->SetInitialStepMultiplier(paramDiffeos->GetInitialStepMultiplier());
	atlasbuilder->SetNumberOfThreads(paramDiffeos->GetNumberOfThreads());

	atlasbuilder->Update();

	timer.Stop();

	std::cout << "Atlas estimation took "
			<< timer.GetElapsedHours() << " hours, "
			<< timer.GetElapsedMinutes() << " minutes, "
			<< timer.GetElapsedSeconds() << " seconds"
			<< std::endl;


	std::cout << "Write output files" << std::endl;

	std::vector<std::string> outname(numObjects);
	std::vector<std::string> outext(numObjects);
	std::vector<std::string> outfn(numObjects);
	for (int i = 0; i < numObjects; i++)
	{
		std::string sourcefn;
		sourcefn.assign(templatefnList[i]);
		int length = sourcefn.size();
		int index = length - 1;
		while( (sourcefn[index] != '.') && (index>=0) )
			index--;

		outname[i] = sourcefn.substr(0,index);
		outext[i] = sourcefn.substr(index,length-index);

		std::ostringstream oss;
		oss << outname[i] << "_template" << outext[i] << std::ends;
		outfn[i] = oss.str();
	}
	atlasbuilder->WriteFlow(outname, outext);

	DeformableMultiObjectType* templ = atlasbuilder->GetTemplate();
	templ->WriteDeformedMultiObjectAt(0, outfn);

	delete atlasbuilder;

}

#endif
