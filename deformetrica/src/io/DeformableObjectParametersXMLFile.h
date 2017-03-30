/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectParametersXMLFile_h
#define _DeformableObjectParametersXMLFile_h

#include "itkXMLFile.h"

#include "DeformableObjectParameters.h"

class DeformableObjectParametersXMLFileReader: public itk::XMLReader<DeformableObjectParameters>
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef DeformableObjectParametersXMLFileReader Self;
	typedef itk::XMLReader<DeformableObjectParameters> Superclass;
	typedef itk::SmartPointer<Self> Pointer;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// RTTI
	itkTypeMacro(Self, Superclass)

	itkNewMacro(Self);

	virtual int CanReadFile(const char* name);



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObjectParametersXMLFileReader();

	~DeformableObjectParametersXMLFileReader();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void StartElement(const char * name, const char **atts);
	virtual void EndElement(const char *name);
	virtual void CharacterDataHandler(const char *inData, int inLength);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObjectParameters::Pointer m_PObject;
	std::string m_CurrentString;


}; /* class DeformableObjectParametersXMLFileReader */





class DeformableObjectParametersXMLFileWriter: public itk::XMLWriterBase<DeformableObjectParameters>
{

public:


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef DeformableObjectParametersXMLFileWriter Self;
	typedef itk::XMLWriterBase<DeformableObjectParameters> Superclass;
	typedef itk::SmartPointer<Self> Pointer;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// RTTI
	itkTypeMacro(Self, Superclass)

	itkNewMacro(Self);

	virtual int CanWriteFile(const char* name);

	virtual int WriteFile();



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObjectParametersXMLFileWriter();

	~DeformableObjectParametersXMLFileWriter();



}; /* class DeformableObjectParametersXMLFileWriter */



// Convenience functions
DeformableObjectParameters::Pointer readDeformableObjectParametersXML(const char* fn);
bool writeDeformableObjectParametersXML(const char* fn, DeformableObjectParameters* p);



#endif /* _DeformableObjectParametersXMLFile_h */
