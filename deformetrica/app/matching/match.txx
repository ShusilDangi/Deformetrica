/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _match_txx
#define _match_txx

#include "DeformetricaConfig.h"

#include "KernelFactory.h"
#include "SparseDiffeoMatcher.h"
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
void match(
		SparseDiffeoParameters* paramDiffeos, int numObjects, std::vector<DeformableObjectParameters::Pointer> paramObjectsList,  const std::vector<char*> sourcefnList, const std::vector<char*> targetfnList)
{
#ifdef USE_DOUBLE_PRECISION
	typedef double TScalar;
	std::cout << "(Computations are in double precision)" << std::endl << std::endl;
#else
	typedef float TScalar;
	std::cout << "(Computations are in simple precision)" << std::endl << std::endl;
#endif
	std::cout << "Sparse diffeomorphic matching\n===" << std::endl;
	std::cout << "Number of objects detected: " << numObjects << std::endl;
	std::cout << "\n===\n" << std::endl;
	std::cout << "Deformations Parameters:" << std::endl;
	paramDiffeos->PrintSelf(std::cout);
	std::cout << "\n===\n" << std::endl;
	for (int i = 0; i < numObjects; i++)
	{
		std::cout << "Object " << i << ": " << std::endl;
		std::cout << "Source file: " << sourcefnList[i] << std::endl;
		std::cout << "Target file: " << targetfnList[i] << std::endl;
		paramObjectsList[i]->PrintSelf(std::cout);
		std::cout << "\n===\n" << std::endl;
	}

#if ITK_VERSION_MAJOR >= 4
	itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_ESTIMATE);
#endif

	typedef itk::Image<TScalar, Dimension> ImageType;
	typedef typename ImageType::Pointer ImageTypePointer;

	typedef Diffeos<TScalar, Dimension> DiffeosType;

	double minI = 0;
	double maxI = 255;

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
	if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "p3m") == 0)
	{
		kfac->UseP3MKernel();
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetKernelType().c_str(), "fgt") == 0)
	{
		kfac->UseFGTKernel();
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



	// create source and target
	vnl_matrix<TScalar> bbox;

	DeformableObjectList sourceObjectList(numObjects);
	DeformableObjectList targetObjectList(numObjects);
	for (int i = 0; i < numObjects; i++)
	{
		typedef DeformableObjectReader<TScalar,Dimension> ObjectReaderType;
		ObjectReaderType* reader = new ObjectReaderType();
		reader->SetObjectParameters(paramObjectsList[i]);

		reader->SetFileName(targetfnList[i]);
		reader->Update();

		targetObjectList[i] = reader->GetOutput();
		bbox = reader->GetBoundingBox();
		for (unsigned int d=0; d<Dimension; d++)
		{
			minData[d] = (bbox(d,0) < minData[d])?bbox(d,0):minData[d];
			maxData[d] = (bbox(d,1) > maxData[d])?bbox(d,1):maxData[d];
		}

		reader->SetFileName(sourcefnList[i]);
		reader->SetTemplateType();
		reader->Update();

		sourceObjectList[i] = reader->GetOutput();
		bbox = reader->GetBoundingBox();
		for (unsigned int d=0; d<Dimension; d++)
		{
			minData[d] = (bbox(d,0) < minData[d])?bbox(d,0):minData[d];
			maxData[d] = (bbox(d,1) > maxData[d])?bbox(d,1):maxData[d];
		}
	}


	DeformableMultiObjectType* source = new DeformableMultiObjectType();
	DeformableMultiObjectType* target = new DeformableMultiObjectType();
	source->SetObjectList(sourceObjectList);
	target->SetObjectList(targetObjectList);

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
	
	// create the matcher!
	typedef SparseDiffeoMatcher<TScalar, Dimension> SparseMatcher;
	SparseMatcher* matcher = new SparseMatcher();
	matcher->SetTarget(target);
	matcher->SetSource(source);
	matcher->SetDiffeos(def);
	matcher->SetInitialCPSpacing(InitialCPSpacing);
	paramDiffeos->FreezeCP()?matcher->SetFreezeCP():matcher->UnsetFreezeCP();
	matcher->SetInitialCPPosition(paramDiffeos->GetInitialCPPosition_fn()); // null string means no InitialCPPosition given (will be initialized with a regular lattice)
	matcher->SetInitialMomentas(paramDiffeos->GetInitialMomenta_fn()); // null string means no InitialMomenta given (will be initialized with zero)
	matcher->SetSparsityPrior(paramDiffeos->GetSparsityPrior());
	if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "gradientdescent") == 0)
	{
		matcher->SetGradientDescent();
	}
	else if (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "ista") == 0)
	{
		matcher->SetISTA();
	}
	else
	{
		if ( (itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "f_ista") != 0) &&
				(itksys::SystemTools::Strucmp(paramDiffeos->GetOptimizationMethodType().c_str(), "fista") != 0) )
			std::cerr << "Unknown optimization method type: defaulting to FISTA" << std::endl;
			matcher->SetFISTA();
	}
	matcher->SetMaxIterations(paramDiffeos->GetMaxIterations());
	matcher->SetMaxLineSearchIterations(paramDiffeos->GetMaxLineSearchIterations());
	matcher->SetAdaptiveTolerance(paramDiffeos->GetAdaptiveTolerance());
	matcher->SetInitialStepMultiplier(paramDiffeos->GetInitialStepMultiplier());
	matcher->Update();

	timer.Stop();

	std::cout << "Matching took "
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
		sourcefn.assign(sourcefnList[i]);
		int length = sourcefn.size();
		int index = length - 1;
		while( (sourcefn[index] != '.') && (index>=0) )
			index--;

		outname[i] = sourcefn.substr(0,index);
		outext[i] = sourcefn.substr(index,length-index);

		std::ostringstream oss;
		oss << outname[i] << "_final" << outext[i] << std::ends;
		outfn[i] = oss.str();
	}

	matcher->WriteFlow(outname, outext);


	delete matcher;


}

#endif
