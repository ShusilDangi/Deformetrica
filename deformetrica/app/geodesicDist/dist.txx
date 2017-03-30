/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _dist_txx
#define _dist_txx

#include "KernelFactory.h"

#include "SparseDiffeoParametersXMLFile.h"

#include <cstring>
#include <iostream>
#include <sstream>

typedef vnl_matrix<float> VNLMatrixType;
typedef float TScalar;


template <unsigned int Dimension>
float dist(SparseDiffeoParameters* paramDiffeos,
		vnl_matrix<float>& CP0,
		vnl_matrix<float>& MOM0,
		const std::vector<char*> objectfnList)
{
	// std::cout << "Compute Geodesic Distance\n===" << std::endl;
	// std::cout << "Control Points File: " << objectfnList[1] << std::endl;
	// std::cout << "Moments File: " << objectfnList[2] << std::endl;
	// std::cout << "\n===\n" << std::endl;
	// std::cout << "Deformations Parameters:" << std::endl;
	// paramDiffeos->PrintSelf(std::cout);
	// std::cout << "\n===\n" << std::endl;
	// std::cout << "Source file: " << objectfnList[0] << std::endl;
	// std::cout << "\n===\n" << std::endl;
	// std::cout << "number of CPs: " << CP0.rows() << std::endl;

	if (CP0.rows() != MOM0.rows())
		throw std::runtime_error("Number of CPs and Momentas mismatched");

	// compute the geodesic distance: the squared norm of the initial velocity speed
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
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

	KernelType* momKernelObj  = kfac->CreateKernelObject();
	momKernelObj->SetKernelWidth(paramDiffeos->GetKernelWidth());
	std::cout << "Kernel Width: " << paramDiffeos->GetKernelWidth() <<std::endl;
	momKernelObj->SetSources(CP0);
	momKernelObj->SetWeights(MOM0);
	VNLMatrixType kMom = momKernelObj->Convolve(CP0);

	TScalar dist = 0.0;
	for (unsigned int i = 0; i < CP0.rows(); i++)
	{
		dist += dot_product(kMom.get_row(i), MOM0.get_row(i));
		std::cout << kMom.get_row(i) << std::endl;
	}

	delete momKernelObj;

	return dist;

}
#endif
