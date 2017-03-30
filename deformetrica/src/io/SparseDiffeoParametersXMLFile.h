/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoParametersXMLFile_h
#define _SparseDiffeoParametersXMLFile_h

#include "itkXMLFile.h"

#include "SparseDiffeoParameters.h"

class SparseDiffeoParametersXMLFileReader: public itk::XMLReader<SparseDiffeoParameters>
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef SparseDiffeoParametersXMLFileReader Self;
	typedef itk::XMLReader<SparseDiffeoParameters> Superclass;
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

	SparseDiffeoParametersXMLFileReader();

	~SparseDiffeoParametersXMLFileReader();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void StartElement(const char * name, const char **atts);
	virtual void EndElement(const char *name);
	virtual void CharacterDataHandler(const char *inData, int inLength);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	SparseDiffeoParameters::Pointer m_PObject;
	std::string m_CurrentString;


}; /* class SparseDiffeoParametersXMLFileReader */





class SparseDiffeoParametersXMLFileWriter: public itk::XMLWriterBase<SparseDiffeoParameters>
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	typedef SparseDiffeoParametersXMLFileWriter Self;
	typedef itk::XMLWriterBase<SparseDiffeoParameters> Superclass;
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

	SparseDiffeoParametersXMLFileWriter();

	~SparseDiffeoParametersXMLFileWriter();



}; /* class SparseDiffeoParametersXMLFileWriter */





// Convenience functions
SparseDiffeoParameters::Pointer readSparseDiffeoParametersXML(const char* fn);
bool writeSparseDiffeoParametersXML(const char* fn, SparseDiffeoParameters* p);




#endif /* _SparseDiffeoParametersXMLFile_h */
