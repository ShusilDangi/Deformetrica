/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformetricaConfig.h"

#include "KernelFactory.h"
#include "SparseDiffeoMatcher.h"

#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "deform.txx"
#include "readMatrixDLM.txx"

#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
	if (argc < 7 || ((argc - 5) % 2 != 0))
	{
		std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml Direction CP.txt MOM.txt paramsObject1.xml object1 paramsObject2.xml object2 ... " << std::endl;
		std::cerr << "where Direction = -1 for using the inverse flow, +1 otherwise" << std::endl;
		return -1;
	}

	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

	int numObjects = (argc - 5) / 2;

	// read general parameters for diffeomorphic matching
	SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
	if (paramDiffeos.IsNull())
	{
		std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
		return -1;
	}

	int dir = atoi(argv[2]);
	if (dir == 1)
		std::cout << "Use the forward flow of diffeomorphisms" << std::endl;
	else if (dir == -1)
		std::cout << "Use the backward flow of diffeomorphisms" << std::endl;
	else
	{
		std::cout << "Unknown value for direction: " << dir << " is neither +1 or -1" << std::endl;
		return -1;
	}

	bool useInverseFlow = (dir == -1);

	// read initial CPs and Momentas
	typedef vnl_matrix<float> VNLMatrixType;
	VNLMatrixType CP0 = readMatrixDLM<float>(argv[3]);
	VNLMatrixType MOM0 = readMatrixDLM<float>(argv[4]);

	// read the list of objects: params, source and target
	std::vector<char*> objectfn(numObjects);
	std::vector<DeformableObjectParameters::Pointer> paramObjects(numObjects);

	for (unsigned int i = 0; i < numObjects; i++)
	{
		paramObjects[i] = readDeformableObjectParametersXML(argv[5 + 2 * i]);
		if (paramObjects[i].IsNull())
		{
			std::cerr << "Failed creating XML object, bad input file " << argv[5 + 2 * i] << "?" << std::endl;
			return -1;
		}
		objectfn[i] = argv[6 + 2 * i];
	}

	try
	{
		deform <2> (paramDiffeos, useInverseFlow, CP0, MOM0, numObjects, paramObjects, objectfn);
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
		return -1;
	}
	catch (std::exception& e)
	{
		std::cerr << "Exception: " << e.what() << std::endl;
		return -1;
	}
	catch (std::string& s)
	{
		std::cerr << "Exception: " << s << std::endl;
		return -1;
	}
	catch (...)
	{
		std::cerr << "Unknown exception" << std::endl;
		return -1;
	}

	return 0;
}
