/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectParameters_h
#define _DeformableObjectParameters_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>


class DeformableObjectParameters : public itk::Object
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef DeformableObjectParameters Self;
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
	itkGetMacro(DeformableObjectType, std::string);
	itkSetMacro(DeformableObjectType, std::string);

	itkGetMacro(DataSigma, double);
	itkSetMacro(DataSigma, double);

	itkGetMacro(SmoothnessPrior, double);
	itkSetMacro(SmoothnessPrior, double);

	// compulsory only if DeformableObjectType is a current
	itkGetMacro(KernelWidth, double);
	itkSetMacro(KernelWidth, double);

	itkGetMacro(KernelType, std::string);
	itkSetMacro(KernelType, std::string);

	itkGetMacro(P3MWorkingSpacingRatio, double);
	itkSetMacro(P3MWorkingSpacingRatio, double);

	itkGetMacro(P3MPaddingFactor, double);
	itkSetMacro(P3MPaddingFactor, double);


	inline bool ReOrient() { return m_reOrient; }
	inline void SetReOrient() { m_reOrient = true; }
	inline void UnsetReOrient() { m_reOrient = false; }

protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObjectParameters();

	~DeformableObjectParameters();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string m_DeformableObjectType;

	double m_DataSigma;
	double m_SmoothnessPrior;

	std::string m_KernelType;

	double m_KernelWidth;
	double m_P3MWorkingSpacingRatio;
	double m_P3MPaddingFactor;

	bool m_reOrient;


}; /* class DeformableObjectParameters */


#endif /* _DeformableObjectParameters_h */
