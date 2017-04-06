/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoParametersXMLFile_cxx
#define _SparseDiffeoParametersXMLFile_cxx

#include "SparseDiffeoParametersXMLFile.h"

#include "itksys/SystemTools.hxx"

#include <fstream>
#include <sstream>
#include <string>

#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

SparseDiffeoParametersXMLFileReader
::SparseDiffeoParametersXMLFileReader()
{
	m_PObject = 0;
}



SparseDiffeoParametersXMLFileReader
::~SparseDiffeoParametersXMLFileReader()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

int
SparseDiffeoParametersXMLFileReader
::CanReadFile(const char* name)
{
	if(!itksys::SystemTools::FileExists(name) ||
			itksys::SystemTools::FileIsDirectory(name) ||
			itksys::SystemTools::FileLength(name) == 0)
		return 0;
	return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
SparseDiffeoParametersXMLFileReader
::StartElement(const char* name, const char** atts)
{
	m_CurrentString = "";

	if(itksys::SystemTools::Strucmp(name,"SPARSE-DIFFEO-PARAMETERS") == 0)
	{
		m_PObject = SparseDiffeoParameters::New();
	}
}


void
SparseDiffeoParametersXMLFileReader
::EndElement(const char* name)
{
	if(itksys::SystemTools::Strucmp(name,"SPARSE-DIFFEO-PARAMETERS") == 0)
	{
		m_OutputObject = &(*m_PObject);
	}
	else if(itksys::SystemTools::Strucmp(name,"NUMBER-OF-TIMEPOINTS") == 0)
	{
		long n = atol(m_CurrentString.c_str());
		m_PObject->SetNumberOfTimePoints(n);
	}
	else if(itksys::SystemTools::Strucmp(name,"OPTIMIZATION-METHOD-TYPE") == 0)
	{
		m_PObject->SetOptimizationMethodType(m_CurrentString);
	}
	else if(itksys::SystemTools::Strucmp(name,"MAX-ITERATIONS") == 0)
	{
		long n = atol(m_CurrentString.c_str());
		m_PObject->SetMaxIterations(n);
	}
	else if(itksys::SystemTools::Strucmp(name,"MAX-LINE-SEARCH-ITERATIONS") == 0)
	{
		long n = atol(m_CurrentString.c_str());
		m_PObject->SetMaxLineSearchIterations(n);
	}
	else if(itksys::SystemTools::Strucmp(name,"INITIAL-CP-SPACING") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetInitialCPSpacing(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"CP-AT-SHAPE-POINTS") == 0)
	{
		if(itksys::SystemTools::Strucmp(m_CurrentString.c_str(),"ON") == 0)
			m_PObject->SetCPsAtShapePoints(true);
		else
			m_PObject->SetCPsAtShapePoints(false);
	}
	else if(itksys::SystemTools::Strucmp(name,"FREEZE-CP") == 0)
	{
		if(itksys::SystemTools::Strucmp(m_CurrentString.c_str(),"ON") == 0)
			m_PObject->SetFreezeCP();
		else if (itksys::SystemTools::Strucmp(m_CurrentString.c_str(),"OFF") == 0)
			m_PObject->UnsetFreezeCP();
		else
			std::cout << "WARNING: Cannot read: " << m_CurrentString.c_str() << ", should be On or Off" << std::endl;
	}
	else if(itksys::SystemTools::Strucmp(name,"INITIAL-CP-POSITION") == 0)
	{
		m_PObject->SetInitialCPPosition_fn(m_CurrentString);
	}
	else if(itksys::SystemTools::Strucmp(name,"INITIAL-MOMENTA") == 0)
	{
		m_PObject->SetInitialMomenta_fn(m_CurrentString);
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
	else if(itksys::SystemTools::Strucmp(name,"SPARSITY-PRIOR") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetSparsityPrior(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"KERNEL-WIDTH") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetKernelWidth(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"SMOOTHING-KERNEL-WIDTH-RATIO") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetSmoothingKernelWidthRatio(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"STEP-EXPAND") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetStepExpand(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"STEP-SHRINK") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetStepShrink(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"ADAPTIVE-TOLERANCE") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetAdaptiveTolerance(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"INITIAL-STEP-MULTIPLIER") == 0)
	{
		double d = atof(m_CurrentString.c_str());
		m_PObject->SetInitialStepMultiplier(d);
	}
	else if(itksys::SystemTools::Strucmp(name,"NUMBER-OF-THREADS") == 0)
	{
		long n = atol(m_CurrentString.c_str());
		m_PObject->SetNumberOfThreads(n);
	}
	else if(itksys::SystemTools::Strucmp(name,"WEIGHTS") == 0)
	{
		m_PObject->SetWeights(m_CurrentString);
	}

}


void
SparseDiffeoParametersXMLFileReader
::CharacterDataHandler(const char* inData, int inLength)
{
	for (int i = 0; i < inLength; i++)
		m_CurrentString += inData[i];
}





////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

SparseDiffeoParametersXMLFileWriter
::SparseDiffeoParametersXMLFileWriter()
{

}



SparseDiffeoParametersXMLFileWriter
::~SparseDiffeoParametersXMLFileWriter()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

int
SparseDiffeoParametersXMLFileWriter
::CanWriteFile(const char* name)
{
	return true;
}



// Support function
template <typename T>
static void
WriteField(SparseDiffeoParametersXMLFileWriter* writer, const char* attname,
		T value, std::ofstream& output)
{
	writer->WriteStartElement(attname, output);
	output << value;
	writer->WriteEndElement(attname, output);
	output << std::endl;
}



int
SparseDiffeoParametersXMLFileWriter
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
	WriteStartElement("!DOCTYPE SPARSE-DIFFEO-PARAMETERS",output);
	output << std::endl;

	WriteStartElement("SPARSE-DIFFEO-PARAMETERS", output);
	output << std::endl;

	SparseDiffeoParameters::Pointer p = m_InputObject;

	WriteField<double>(this, "KERNEL-WIDTH", p->GetKernelWidth(), output);

	WriteField<std::string>(this, "KERNEL-TYPE", p->GetKernelType(), output);

	WriteField<unsigned int>(this, "NUMBER-OF-TIMEPOINTS", p->GetNumberOfTimePoints(), output);

	WriteField<double>(this, "INITIAL-CP-SPACING", p->GetInitialCPSpacing(), output);
	WriteField<std::string>(this, "INITIAL-CP-POSITION", p->GetInitialCPPosition_fn(), output);
	WriteField<std::string>(this, "INITIAL-MOMENTA", p->GetInitialMomenta_fn(), output);
	WriteField<const char*>(this, "FREEZE-CP", p->FreezeCP()?"On":"Off", output);

	WriteField<double>(this, "P3M-WORKING-SPACING-RATIO", p->GetP3MWorkingSpacingRatio(), output);
	WriteField<double>(this, "P3M-PADDING-FACTOR", p->GetP3MPaddingFactor(), output);

	WriteField<double>(this, "SPARSITY-PRIOR", p->GetSparsityPrior(), output);

	WriteField<std::string>(this, "OPTIMIZATION-METHOD-TYPE", p->GetOptimizationMethodType(), output);

	WriteField<unsigned int>(this, "MAX-ITERATIONS", p->GetMaxIterations(), output);
	WriteField<unsigned int>(this, "MAX-LINE-SEARCH-ITERATIONS", p->GetMaxLineSearchIterations(), output);

	WriteField<double>(this, "STEP-EXPAND", p->GetStepExpand(), output);
	WriteField<double>(this, "STEP-SHRINK", p->GetStepShrink(), output);
	WriteField<double>(this, "ADAPTIVE-TOLERANCE", p->GetAdaptiveTolerance(), output);
	WriteField<double>(this, "INITIAL-STEP-MULTIPLIER", p->GetInitialStepMultiplier(), output);

	WriteField<double>(this, "SMOOTHING-KERNEL-WIDTH-RATIO", p->GetSmoothingKernelWidthRatio(), output);
	WriteField<unsigned int>(this, "NUMBER-OF-THREADS", p->GetNumberOfThreads(), output);

	WriteField<std::string>(this, "REGULARITY-WEIGHTS", p->GetWeights(), output);

	// Finish
	WriteEndElement("SPARSE-DIFFEO-PARAMETERS", output);
	output << std::endl;
	output.close();

	return 0;
}





// Definition of some convenience functions
SparseDiffeoParameters::Pointer
readSparseDiffeoParametersXML(const char* fn)
{
	SparseDiffeoParametersXMLFileReader::Pointer reader =
			SparseDiffeoParametersXMLFileReader::New();
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
writeSparseDiffeoParametersXML(const char* fn, SparseDiffeoParameters* p)
{
	if (p == 0)
		return false;
	if (!p->CheckValues())
		return false;

	SparseDiffeoParametersXMLFileWriter::Pointer writer =
			SparseDiffeoParametersXMLFileWriter::New();

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



#endif /* _SparseDiffeoParametersXMLFile_cxx */
