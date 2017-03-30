/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObjectParametersXMLFile_cxx
#define _DeformableObjectParametersXMLFile_cxx

#include "DeformableObjectParametersXMLFile.h"

#include "itksys/SystemTools.hxx"

#include <fstream>
#include <sstream>
#include <string>

#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

DeformableObjectParametersXMLFileReader
::DeformableObjectParametersXMLFileReader()
{
	m_PObject = 0;
}



DeformableObjectParametersXMLFileReader
::~DeformableObjectParametersXMLFileReader()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

int
DeformableObjectParametersXMLFileReader
::CanReadFile(const char* name)
{
	if(!itksys::SystemTools::FileExists(name) ||
			itksys::SystemTools::FileIsDirectory(name) ||
			itksys::SystemTools::FileLength(name) == 0)
		return 0;
	return 1;
}



void
DeformableObjectParametersXMLFileReader
::StartElement(const char* name, const char** atts)
{
	m_CurrentString = "";

	if(itksys::SystemTools::Strucmp(name,"DEFORMABLE-OBJECT-PARAMETERS") == 0)
	{
		m_PObject = DeformableObjectParameters::New();
	}
}



void
DeformableObjectParametersXMLFileReader
::EndElement(const char* name)
{
	if(itksys::SystemTools::Strucmp(name,"DEFORMABLE-OBJECT-PARAMETERS") == 0)
	{
		m_OutputObject = &(*m_PObject);
	}
	if(itksys::SystemTools::Strucmp(name,"DEFORMABLE-OBJECT-TYPE") == 0)
	{
		m_PObject->SetDeformableObjectType(m_CurrentString);
	}
	else if(itksys::SystemTools::Strucmp(name,"KERNEL-TYPE") == 0)
	{
		m_PObject->SetKernelType(m_CurrentString);
	}
	else if(itksys::SystemTools::Strucmp(name,"P3M-WORKING-SPACING-RATIO") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetP3MWorkingSpacingRatio(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"P3M-PADDING-FACTOR") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetP3MPaddingFactor(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"DATA-SIGMA") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetDataSigma(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"SMOOTHNESS-PRIOR") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetSmoothnessPrior(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"KERNEL-WIDTH") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetKernelWidth(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"REORIENT-NORMALS") == 0)
	{
		if (itksys::SystemTools::Strucmp(m_CurrentString.c_str(),"On") == 0)
			m_PObject->SetReOrient();
		else if (itksys::SystemTools::Strucmp(m_CurrentString.c_str(),"Off") == 0)
			m_PObject->UnsetReOrient();
		else
			std::cout << "WARNING: Cannot read: " << m_CurrentString.c_str() << ", should be On or Off" << std::endl;
	}

}



void
DeformableObjectParametersXMLFileReader
::CharacterDataHandler(const char* inData, int inLength)
{
	for (int i = 0; i < inLength; i++)
		m_CurrentString += inData[i];
}





////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

DeformableObjectParametersXMLFileWriter
::DeformableObjectParametersXMLFileWriter()
{

}



DeformableObjectParametersXMLFileWriter
::~DeformableObjectParametersXMLFileWriter()
{

}





////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

int
DeformableObjectParametersXMLFileWriter
::CanWriteFile(const char* name)
{
	return true;
}



// Support function
template <typename T>
static void
WriteField(DeformableObjectParametersXMLFileWriter* writer, const char* attname,
		T value, std::ofstream& output)
{
	writer->WriteStartElement(attname, output);
	output << value;
	writer->WriteEndElement(attname, output);
	output << std::endl;
}



int
DeformableObjectParametersXMLFileWriter
::WriteFile()
{
	if (m_InputObject == 0)
		itkExceptionMacro(<< "No object to write");

	if (!m_InputObject->CheckValues())
		itkExceptionMacro(<< "Invalid values");

	if (m_Filename.length() == 0)
		itkExceptionMacro(<< "No file name specified");

	std::ofstream output(m_Filename.c_str());
	if (output.fail())
		itkExceptionMacro(<< "Can not open " << m_Filename);

	// Header
	WriteStartElement("?xml version=\"1.0\"?",output);
	output << std::endl;
	WriteStartElement("!DOCTYPE DEFORMABLE-OBJECT-PARAMETERS",output);
	output << std::endl;

	WriteStartElement("DEFORMABLE-OBJECT-PARAMETERS", output);
	output << std::endl;

	DeformableObjectParameters::Pointer p = m_InputObject;

	WriteField<std::string>(this, "DEFORMABLE-OBJECT-TYPE", p->GetDeformableObjectType(), output);
	WriteField<double>(this, "DATA-SIGMA", p->GetDataSigma(), output);
	WriteField<double>(this, "SMOOTHNESS-PRIOR", p->GetSmoothnessPrior(), output);

	WriteField<std::string>(this, "KERNEL-TYPE", p->GetKernelType(), output);
	WriteField<double>(this, "KERNEL-WIDTH", p->GetKernelWidth(), output);

	WriteField<double>(this, "P3M-WORKING-SPACING-RATIO", p->GetP3MWorkingSpacingRatio(), output);
	WriteField<double>(this, "P3M-PADDING-FACTOR", p->GetP3MPaddingFactor(), output);

	WriteField<const char*>(this, "REORIENT-NORMALS", p->ReOrient()?"On":"Off", output);

	// Finish
	WriteEndElement("DEFORMABLE-OBJECT-PARAMETERS", output);
	output << std::endl;
	output.close();

	return 0;
}





// Definition of some convenience functions
DeformableObjectParameters::Pointer
readDeformableObjectParametersXML(const char* fn)
{
	DeformableObjectParametersXMLFileReader::Pointer reader =
			DeformableObjectParametersXMLFileReader::New();
	try
	{
		reader->SetFilename(fn);
		reader->GenerateOutputInformation();
	}
	catch (...)
	{
		return 0;
	}
	return reader->GetOutputObject();
}



bool
writeDeformableObjectParametersXML(const char* fn, DeformableObjectParameters* p)
{
	if (p == 0)
		return false;
	if (!p->CheckValues())
		return false;

	DeformableObjectParametersXMLFileWriter::Pointer writer =
			DeformableObjectParametersXMLFileWriter::New();

	std::string outfn = fn;

	/*
	// Enforce XML file extension
	std::string ext = mu::get_ext(fn);
	if (ext.compare("xml") != 0)
		outfn += std::string(".xml");
	*/

	writer->SetFilename(outfn.c_str());
	writer->SetObject(p);
	writer->WriteFile();

	return true;
}


#endif /* _DeformableObjectParametersXMLFile_h */
