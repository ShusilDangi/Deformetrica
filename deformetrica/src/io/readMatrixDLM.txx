/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _readMatrixDLM_txx
#define _readMatrixDLM_txx

#include "vnl/vnl_matrix.h"

#include <fstream>
#include <sstream>
#include <string>

template <class TScalar>
vnl_matrix<TScalar> readMatrixDLM(const char* fn)
{

	std::cout << "Reading " << fn << std::endl;

	std::ifstream infile(fn);

	if (!infile.is_open())
	{
		std::cout << "error opening file" << std::endl; // to re-route error stream to std::cout
		throw std::runtime_error("error opening file");
	}

	unsigned int numRows = 0;

	bool isfirst = true;
	std::string firstline;

	std::string line;
	while (!infile.eof())
	{
		getline(infile, line);

		if (isfirst)
		{
			firstline = line;
			isfirst = false;
		}

		numRows++;
	}

	unsigned int numCols = 0;

	std::stringstream ss(firstline);
	while (!ss.eof())
	{
		float f;
		ss >> f;
		numCols++;
	}

	//infile.clear();
	//infile.seekg(0);
	infile.close();
	infile.open(fn);

	--numRows;

	std::cout << "Reading " << numRows << " x " << numCols << " matrix" << std::endl;
	vnl_matrix<TScalar> M(numRows, numCols, 0);

	TScalar a = 0;
	for (unsigned int i = 0; i < numRows; i++)
	{
		for (unsigned int j = 0; j < numCols; j++)
		{
			infile >> a;
			M(i, j) = a;
		}
	}

	infile.close();

	return M;
}


template <class TScalar>
std::vector< vnl_matrix<TScalar> > readMultipleMatrixDLM(const char* fn)
{

	std::cout << "Reading " << fn << std::endl;

	std::ifstream infile(fn);

	if (!infile.is_open())
	{
		std::cout << "error opening file" << std::endl; // to re-route error stream to std::cout
		throw std::runtime_error("error opening file");
	}

	std::string firstline;
	getline(infile, firstline);

	std::stringstream ss(firstline);
	unsigned int N, numRows, numCols;
	ss >> N;
	ss >> numRows;
	ss >> numCols;
	std::cout << "Reading " << N << " times a " << numRows << " x " << numCols << " matrix" << std::endl;

	getline(infile, firstline);

	std::vector< vnl_matrix<TScalar> > M(N);
	vnl_matrix<TScalar> Mn(numRows, numCols, 0);

	for (unsigned int n = 0; n < N; n++)
	{
		Mn.fill(0);
		TScalar a = 0;

		for (unsigned int i = 0; i < numRows; i++)
		{
			for (unsigned int j = 0; j < numCols; j++)
			{
				infile >> a;
				Mn(i, j) = a;
			}
		}
		M[n] = Mn;
		getline(infile, firstline);
	}

	infile.close();

	return M;
}


#endif /* _readMatrixDLM_txx */
