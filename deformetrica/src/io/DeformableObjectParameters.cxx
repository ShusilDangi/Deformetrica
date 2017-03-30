/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "DeformableObjectParameters.h"

#include "itksys/SystemTools.hxx"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

DeformableObjectParameters
::DeformableObjectParameters()
{
	m_DeformableObjectType = "unknown";

	m_DataSigma = 0.0;

	m_KernelType = "p3m";
	m_P3MWorkingSpacingRatio = 0.2; // 1/5 gives a relative approximation error of about 5%, use 0.3 for increased speed
	m_P3MPaddingFactor = 3.0; // enlarge grids by 3 x kernelwidth to avoid side effects (FFTs have circular boundary conditions). It is also used to define a bounding box

	m_KernelWidth = 0.0;

	m_reOrient = false;
}



DeformableObjectParameters
::~DeformableObjectParameters()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

bool
DeformableObjectParameters
::CheckValues()
{
	if (m_DataSigma < 1e-20)
		return false;

	if (itksys::SystemTools::Strucmp(m_DeformableObjectType.c_str(),"Landmarks") != 0)
	{
		if (m_KernelWidth < 1e-20)
			return false;
	}

	return true;
}



void
DeformableObjectParameters
::PrintSelf(std::ostream& os)
{
	os << "Deformable object type: " << m_DeformableObjectType << std::endl;
	os << "Data sigma = " << m_DataSigma << std::endl;
	os << std::endl;
	os << "Kernel type = " << m_KernelType << std::endl;
	os << "Kernel width = " << m_KernelWidth << std::endl;
	os << std::endl;
	os << "P3M working spacing ratio = " << m_P3MWorkingSpacingRatio << std::endl;
	os << "P3M padding factor = " << m_P3MPaddingFactor << std::endl;
	os << "Reorient normals: " << (m_reOrient?"On":"Off") << std::endl;
	os << std::endl;
}
