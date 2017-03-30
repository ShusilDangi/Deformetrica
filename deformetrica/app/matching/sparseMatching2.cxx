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

#include "match.txx"

#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
	if (argc < 5 || ( (argc-2) % 3 != 0 ) )
	{
		std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml paramsObject1.xml source1 target1 paramsObject2.xml source2 target2 ... " << std::endl;
		return -1;
	}

	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

	int numObjects = (argc - 2) / 3;

	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

	// read general parameters for diffeomorphic matching
	SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
	if (paramDiffeos.IsNull())
	{
		std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
		return -1;
	}

	// read the list of objects: params, source and target
	std::vector<char*> sourcefn;
	std::vector<char*> targetfn;
	std::vector<DeformableObjectParameters::Pointer> paramObjects;

	sourcefn.resize(numObjects);
	targetfn.resize(numObjects);
	paramObjects.resize(numObjects);

	for (unsigned int i = 0; i < numObjects; i++)
	{
		paramObjects[i] = readDeformableObjectParametersXML(argv[2 + 3*i]);
		if (paramObjects[i].IsNull())
		{
			std::cerr << "Failed creating XML object, bad input file " << argv[2 + 3*i] << "?" << std::endl;
			return -1;
		}
		sourcefn[i] = argv[3 + 3*i];
		targetfn[i] = argv[4 + 3*i];
	}


	try
	{
		match<2>(paramDiffeos, numObjects, paramObjects, sourcefn, targetfn);
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
