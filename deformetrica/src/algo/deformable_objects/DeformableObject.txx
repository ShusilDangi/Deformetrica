/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObject_txx
#define _DeformableObject_txx

#include "DeformableObject.h"

#include "KernelFactory.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
DeformableObject<TScalar, Dimension>
::DeformableObject()
 {
	m_OutOfBox = false;
	m_Modified = false;
	m_Type = null;
 }



template <class TScalar, unsigned int Dimension>
DeformableObject<TScalar, Dimension>
::~DeformableObject()
{ 

}



template <class TScalar, unsigned int Dimension>
DeformableObject<TScalar, Dimension>
::DeformableObject(const DeformableObject& other)
 {
	m_Type = other.m_Type;

	m_DataSigmaSquared = other.m_DataSigmaSquared;

	m_Def = other.m_Def;

	m_OutOfBox = other.m_OutOfBox;
	m_Modified = other.m_Modified;
 }


#endif /* _DeformableObject_txx */
