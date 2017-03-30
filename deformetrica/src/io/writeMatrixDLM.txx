/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _writeMatrixDLM_txx
#define _writeMatrixDLM_txx

#include "vnl/vnl_matrix.h"

#include <fstream>



template <class TScalar>
void writeMatrixDLM(const char* fn, const vnl_matrix<TScalar>& M)
{
	unsigned int numRows = M.rows();
	unsigned int numCols = M.columns();

	std::ofstream outfile(fn);

	for (unsigned int i = 0; i < numRows; i++)
	{
		for (unsigned int j = 0; j < (numCols-1); j++)
			outfile << M(i, j) << " ";
		outfile << M(i, numCols-1) << std::endl;
	}

	outfile.close();
}



template <class TScalar>
void writeMultipleMatrixDLM(const char* fn, const std::vector< vnl_matrix<TScalar> >& M)
{
	std::ofstream outfile(fn);

	unsigned int N = M.size();
	if ( N==0 )
		outfile.close();

	unsigned int numRows = M[0].rows();
	unsigned int numCols = M[0].columns();

	outfile << N << " " << numRows << " " << numCols << std::endl << std::endl;

	for (unsigned int n = 0; n < N; n++)
	{
		vnl_matrix<TScalar> Mn = M[n];
		for (unsigned int i = 0; i < numRows; i++)
		{
			for (unsigned int j = 0; j < (numCols-1); j++)
				outfile << Mn(i, j) << " ";
			outfile << Mn(i, numCols-1) << std::endl;
		}
		outfile << std::endl;
	}

	outfile.close();
}


#endif /* _writeMatrixDLM_txx */
