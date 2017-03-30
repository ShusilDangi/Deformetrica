/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*    																				   *
*	 Author: Shusil Dangi, Rochester Institute of Technology / Kitware                 *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformetricaConfig.h"

#include "dist.txx"
#include "readMatrixDLM.txx"

#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
	if (argc < 4 )
	{
		std::cerr << "Usage: " << argv[0] << " paramsDiffeos.xml CP.txt MOM.txt " << std::endl;
		return -1;
	}

	std::cout << "Deformetrica " << DEFORMETRICA_VERSION_MAJOR << "." << DEFORMETRICA_VERSION_MINOR << std::endl;

	// read general parameters for diffeomorphic matching
	SparseDiffeoParameters::Pointer paramDiffeos = readSparseDiffeoParametersXML(argv[1]);
	if (paramDiffeos.IsNull())
	{
		std::cerr << "Failed creating XML object, bad input file " << argv[1] << "?" << std::endl;
		return -1;
	}

	// read initial CPs and Momentas
	VNLMatrixType CP0 = readMatrixDLM<float>(argv[2]);
	std::vector< vnl_matrix<TScalar> > MOM0 = readMultipleMatrixDLM<float>(argv[3]);

	//read the list of objects: params, cp, mom files
	std::vector<char*> objectfn(argc-1);
	for (unsigned int i = 1; i < argc; i++)
	{
		objectfn[i-1] = argv[i];
	}

	try
	{
		float geoDist;
		for(unsigned int i = 0; i<MOM0.size(); i++)
		{
			geoDist = dist<3>(paramDiffeos, CP0, MOM0[i], objectfn);
			std::cout<< "Geodesic Distance: " << geoDist << std::endl;
		}
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
