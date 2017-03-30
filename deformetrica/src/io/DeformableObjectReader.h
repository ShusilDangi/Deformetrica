/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectReader_h
#define _DeformableObjectReader_h

#include "DeformableObject.h"
#include "DeformableObjectParametersXMLFile.h"

#include "vnl/vnl_matrix.h"

#include <cstring>
#include <iostream>
#include <sstream>

template <class TScalar, unsigned int Dimension>
class DeformableObjectReader
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// Matrix type.
	typedef vnl_vector<TScalar> VNLVectorType;

	/// Deformable object type.
	typedef DeformableObject<TScalar,Dimension> DeformableObjectType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObjectReader();

	~DeformableObjectReader();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the file name to \e fn.
	void SetFileName(char* fn) { m_FileName = fn; }


	void SetObjectParameters(DeformableObjectParameters::Pointer param) { m_ParamObject = param; }


	void SetTemplateType() { m_IsTemplate = true; }

	/// Returns the ouput i.e. the deformable object.
	DeformableObjectType* GetOutput() { return m_Object; }

	/// Returns the bounding box.
	VNLMatrixType GetBoundingBox() { return m_Bbox; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void Update();



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Name of the file where the deformable object will be read.
	char* m_FileName;	

	DeformableObjectParameters::Pointer m_ParamObject;

	/// Our deformable object which will be extracted from the file.
	DeformableObjectType* m_Object;


	VNLMatrixType m_Bbox;

	bool m_IsTemplate;

	VNLVectorType m_Min;
	VNLVectorType m_Max;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "DeformableObjectReader.txx"
#endif


#endif
