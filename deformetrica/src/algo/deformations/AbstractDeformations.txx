/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AbstractDeformations_txx
#define _AbstractDeformations_txx

#include "KernelFactory.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
AbstractDeformations<TScalar, Dimension>
::AbstractDeformations() :
 m_NumberOfTimePoints(15), m_KernelWidth(1.0), m_UseImprovedEuler(true), m_Modified(true), m_Type(null)
 {
 }



template <class TScalar, unsigned int Dimension>
AbstractDeformations<TScalar, Dimension>
::~AbstractDeformations()
{}



template <class TScalar, unsigned int Dimension>
AbstractDeformations<TScalar, Dimension>
::AbstractDeformations(const AbstractDeformations& other)
 {
	m_NumberOfTimePoints = other.m_NumberOfTimePoints;

	m_StartPositions = other.m_StartPositions;
	m_StartMomentas = other.m_StartMomentas;

	m_DataDomain = other.m_DataDomain;

	m_KernelWidth = other.m_KernelWidth;
	m_UseImprovedEuler = other.m_UseImprovedEuler;
	m_Modified = true;

	m_PaddingFactor = other.m_PaddingFactor;

	m_Type = other.m_Type;

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
bool
AbstractDeformations<TScalar, Dimension>
::CheckBoundingBox(VNLMatrixList& X, int t)
{
	VNLMatrixType Xt = X[t];
	int outOfBox = 0;

	for (int d=0; d<Dimension; d++)
	{
		outOfBox += (Xt.get_column(d).min_value() < m_BoundingBox(d,0));
		outOfBox += (Xt.get_column(d).max_value() > m_BoundingBox(d,1));
	}

	return (outOfBox > 0);
}



template<class TScalar, unsigned int Dimension>
void
AbstractDeformations<TScalar, Dimension>
::Update()
 {
	if (this->IsModified())
	{
		this->InitBoundingBox();
		this->Shoot();

		this->UnsetModified();
	}
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
AbstractDeformations<TScalar, Dimension>
::InitBoundingBox()
 {
	TScalar h = m_KernelWidth;
	TScalar paddingFactor = m_PaddingFactor;

	TScalar offset = paddingFactor*h/2;

	VNLMatrixType boundingBox = m_DataDomain;
	boundingBox.set_column(0, boundingBox.get_column(0) - offset);
	boundingBox.set_column(1, boundingBox.get_column(1) + offset);

	m_BoundingBox = boundingBox;
	m_OutOfBox = false;
 }



#endif /* _AbstractDeformations_txx */
