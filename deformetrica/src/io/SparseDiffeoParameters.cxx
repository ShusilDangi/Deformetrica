/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoParameters_cxx
#define _SparseDiffeoParameters_cxx

#include "SparseDiffeoParameters.h"

#include "itksys/SystemTools.hxx"
#include "vnl/vnl_matrix.h"

SparseDiffeoParameters
::SparseDiffeoParameters()
{
	m_KernelWidth = 0.0;

	m_KernelType = "p3m";
	m_P3MWorkingSpacingRatio = 0.2; // 1/5 gives a relative approximation error of about 5%, use 0.3 for increased speed
	m_P3MPaddingFactor = 3.0; // enlarge grids by 3 x kernelwidth to avoid side effects (FFTs have circular boundary conditions). It is also used to define a bounding box

	m_NumberOfTimePoints = 10;
	m_InitialCPSpacing = 0.0;

	m_SparsityPrior = 0.0;

	m_OptimizationMethodType = "F_ISTA";

	m_CPsAtShapePoints = false;
	m_FreezeCP = false;
	m_InitialCPPosition_fn = "";
	m_InitialMomenta_fn = "";

	m_SmoothingKernelWidthRatio = 1.0;

	m_MaxIterations = 100;
	m_MaxLineSearchIterations = 10;

	m_StepExpand = 2.0;
	m_StepShrink = 0.5;

	m_AdaptiveTolerance = 1e-4;
	m_InitialStepMultiplier = 0.1;

	m_NumberOfThreads = 1;
}

SparseDiffeoParameters
::~SparseDiffeoParameters()
{

}

bool
SparseDiffeoParameters
::CheckValues()
{
	if (m_KernelWidth < 1e-20)
		return false;

	if (m_NumberOfTimePoints < 2)
		return false;

	if (m_MaxIterations < 1)
		return false;
	if (m_MaxLineSearchIterations < 1)
		return false;

	if (m_AdaptiveTolerance <= 0.0)
		return false;
	if (m_InitialStepMultiplier < 1e-20)
		return false;

	if (m_InitialCPSpacing < 1e-20)
		return false;

	return true;
}

void
SparseDiffeoParameters
::PrintSelf(std::ostream& os)
{
	os << "Kernel width = " << m_KernelWidth << std::endl;
	os << std::endl;
	os << "Kernel type = " << m_KernelType << std::endl;
	os << "P3M working spacing ratio = " << m_P3MWorkingSpacingRatio << std::endl;
	os << "P3M padding factor = " << m_P3MPaddingFactor << std::endl;
	os << std::endl;  
	os << "Number of time points = " <<  m_NumberOfTimePoints << std::endl;
	os << "Initial CP spacing = " << m_InitialCPSpacing << std::endl;
	os << "Freeze CP = " << (m_FreezeCP?"On":"Off") << std::endl;
	os << "Initial set of CP loaded from " << m_InitialCPPosition_fn << std::endl;
	os << "Initial momenta loaded from " << m_InitialMomenta_fn << std::endl;
	os << std::endl;
	os << "Sparsity prior = " << m_SparsityPrior << std::endl;
	os << std::endl;
	os << "SmoothingKernelWidthRatio = " << m_SmoothingKernelWidthRatio << std::endl;
	os << "Optimization method: " << m_OptimizationMethodType << std::endl;
	os << "Max descent iterations = " << m_MaxIterations << std::endl;
	os << "Max line search iterations = " << m_MaxLineSearchIterations << std::endl;
	os << "Step expand = " << m_StepExpand << std::endl;
	os << "Step shrink = " << m_StepShrink << std::endl;
	os << "Adaptive tolerance = " << m_AdaptiveTolerance << std::endl;
	os << "Initial step multiplier = " << m_InitialStepMultiplier << std::endl;
	os << "Number of threads = " << m_NumberOfThreads << std::endl;
	os << std::endl;
}


#endif /* _SparseDiffeoParameters_cxx */
