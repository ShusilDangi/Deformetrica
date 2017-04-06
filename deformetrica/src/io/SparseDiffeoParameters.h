/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoParameters_h
#define _SparseDiffeoParameters_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>


class SparseDiffeoParameters : public itk::Object
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef SparseDiffeoParameters Self;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	itkNewMacro(Self);

	// Make sure all values are OK
	virtual bool CheckValues();

	virtual void PrintSelf(std::ostream& os);

	// Compulsory parameters
	itkGetMacro(KernelWidth, double);
	itkSetMacro(KernelWidth, double);

	// Optional parameters...
	// ... for the diffeos
	itkGetMacro(KernelType, std::string);
	itkSetMacro(KernelType, std::string);

	itkGetMacro(NumberOfTimePoints, unsigned int);
	itkSetMacro(NumberOfTimePoints, unsigned int);

	itkGetMacro(InitialCPSpacing, double);
	itkSetMacro(InitialCPSpacing, double);

	itkGetMacro(P3MWorkingSpacingRatio, double);
	itkSetMacro(P3MWorkingSpacingRatio, double);

	itkGetMacro(P3MPaddingFactor, double);
	itkSetMacro(P3MPaddingFactor, double);


	// ... for the objective function
	itkGetMacro(SparsityPrior, double);
	itkSetMacro(SparsityPrior, double);

	// ... for the optimization method
	itkGetMacro(OptimizationMethodType, std::string);
	itkSetMacro(OptimizationMethodType, std::string);

	itkGetMacro(InitialCPPosition_fn, std::string);
	itkSetMacro(InitialCPPosition_fn, std::string);

	itkGetMacro(InitialMomenta_fn, std::string);
	itkSetMacro(InitialMomenta_fn, std::string);

	inline void SetCPsAtShapePoints(bool yesNo) {m_CPsAtShapePoints = yesNo; }
	inline bool GetCpsAtShapePoints() {return m_CPsAtShapePoints; }

	inline void SetFreezeCP() { m_FreezeCP = true; }
	inline void UnsetFreezeCP() { m_FreezeCP = false; }
	inline bool FreezeCP() { return m_FreezeCP; }

 	itkGetMacro(SmoothingKernelWidthRatio, double);
	itkSetMacro(SmoothingKernelWidthRatio, double);

	itkGetMacro(MaxIterations, unsigned int);
	itkSetMacro(MaxIterations, unsigned int);

	itkGetMacro(MaxLineSearchIterations, unsigned int);
	itkSetMacro(MaxLineSearchIterations, unsigned int);

	itkGetMacro(StepExpand, double);
	itkSetMacro(StepExpand, double);

	itkGetMacro(StepShrink, double);
	itkSetMacro(StepShrink, double);

	itkGetMacro(AdaptiveTolerance, double);
	itkSetMacro(AdaptiveTolerance, double);

	itkGetMacro(InitialStepMultiplier, double);
	itkSetMacro(InitialStepMultiplier, double);

	itkGetMacro(NumberOfThreads, unsigned int);
	itkSetMacro(NumberOfThreads, unsigned int);

	// Required if different weights need to be assigned to different subjects
	itkGetMacro(Weights, std::string);
	itkSetMacro(Weights, std::string);


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	SparseDiffeoParameters();

	~SparseDiffeoParameters();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	double m_KernelWidth;

	std::string m_KernelType;
	unsigned int m_NumberOfTimePoints;
	double m_InitialCPSpacing;
	double m_P3MWorkingSpacingRatio;
	double m_P3MPaddingFactor;

	double m_SparsityPrior;

	std::string m_OptimizationMethodType;

	bool m_CPsAtShapePoints;
	bool m_FreezeCP;
	std::string m_InitialCPPosition_fn;
	std::string m_InitialMomenta_fn;

	double m_SmoothingKernelWidthRatio;

	unsigned int m_MaxIterations;
	unsigned int m_MaxLineSearchIterations;

	double m_StepExpand;
	double m_StepShrink;
	double m_AdaptiveTolerance;
	double m_InitialStepMultiplier;

	unsigned int m_NumberOfThreads;

	std::string m_Weights;


}; /* class SparseDiffeoParameters */

#endif /* _SparseDiffeoParameters_h */
