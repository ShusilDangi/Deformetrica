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

#include "estimateAtlas.txx"

#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
	if (argc < 6 )
	{
		std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml NumberOfObjects paramsObject1.xml InitialTemplate1 Subject1 Subject2 Subject3 ... paramsObject2.xml InitialTemplate2 Subject1 Subject2 Subject3 ... " << std::endl;
		return -1;
	}

	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

	int numObjects = atoi(argv[2]);
	std::cout << "Number of Objects = " << numObjects << std::endl;

	if ( (argc - 3) % numObjects != 0)
	{
		std::cerr << "Number of files mismatched with the number of objects" << std::endl;
		return -1;
	}

	int numSubjects = (argc - 3) / numObjects - 2;
	std::cout << "Number of Subjects = " << numSubjects << std::endl;

	if (numSubjects < 2)
	{
		std::cerr << "Atlas building requires at least 2 subjects" << std::endl;
		return -1;
	}


	// read general parameters for diffeomorphic matching
	SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
	if (paramDiffeos.IsNull())
	{
		std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
		return -1;
	}

	// read the list of objects: params, source and target
	std::vector<char*> templatefn;
	std::vector< std::vector< char* > > subjectfn;
	std::vector<DeformableObjectParameters::Pointer> paramObjects;

	templatefn.resize(numObjects);
	paramObjects.resize(numObjects);
	subjectfn.resize(numObjects);
	for (int i = 0; i < numObjects; i++)
		subjectfn[i].resize(numSubjects);

	unsigned int indx = 3;
	for (unsigned int i = 0; i < numObjects; i++)
	{
		paramObjects[i] = readDeformableObjectParametersXML(argv[indx++]);
		if (paramObjects[i].IsNull())
		{
			std::cerr << "Failed creating XML object, bad input file " << argv[indx] << "?" << std::endl;
			return -1;
		}
		templatefn[i] = argv[indx++];
		for (unsigned s = 0; s < numSubjects; s++)
			subjectfn[i][s] = argv[indx++];
	}

	try
	{
		estimateAtlas<2>(paramDiffeos, numSubjects, numObjects, paramObjects, templatefn, subjectfn);
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
