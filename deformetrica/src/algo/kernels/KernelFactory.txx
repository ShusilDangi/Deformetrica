/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _KernelFactory_txx
#define _KernelFactory_txx

#include "ExactKernel.h"
#ifdef USE_CUDA
	#include "CUDAExactKernel.h"
#endif
#include "FGTKernel.h"
#include "P3MKernel.h"

#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>*
KernelFactory<TScalar, PointDim>::m_SingletonInstance = 0;



template<class TScalar, unsigned int PointDim>
itk::SimpleFastMutexLock
KernelFactory<TScalar, PointDim>::m_Mutex;



template<class TScalar, unsigned int PointDim>
bool
KernelFactory<TScalar, PointDim>::m_IsInstantiated = false;



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>
::KernelFactory()
 {
	m_WhichKernel = Exact;
	m_DataDomain = 0;
	m_WorkingSpacingRatio = 0.0;

 }



template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>
::~KernelFactory()
{
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
void
KernelFactory<TScalar, PointDim>
::SetDataDomain(VNLMatrixType DD)
 {
	m_Mutex.Lock();
	m_DataDomain = DD;
	m_Mutex.Unlock();
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int PointDim>
KernelFactory<TScalar, PointDim>*
KernelFactory<TScalar, PointDim>
::Instantiate()
 {
	if (!m_IsInstantiated)
	{
		m_Mutex.Lock();
		if (m_SingletonInstance == 0)
			m_SingletonInstance = new KernelFactory<TScalar, PointDim>();
		m_Mutex.Unlock();
		m_IsInstantiated = true;
	}

	return m_SingletonInstance;
 }



template<class TScalar, unsigned int PointDim>
typename KernelFactory<TScalar, PointDim>::KernelBaseType*
KernelFactory<TScalar, PointDim>
::CreateKernelObject()
 {
	switch (m_WhichKernel)
	{
	case Exact:
		return new ExactKernel<TScalar, PointDim>();
#ifdef USE_CUDA
	case CUDAExact:
		return new CUDAExactKernel<TScalar, PointDim>();
#endif
	case FGT:
		return new FGTKernel<TScalar, PointDim>();
	case P3M:
		typedef P3MKernel<TScalar, PointDim> P3MKernelType;
		P3MKernelType* obj = new P3MKernelType();
		if (m_DataDomain.size())
			obj->SetDataDomain(this->GetDataDomain());

		obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
		obj->SetPaddingFactor(this->GetPaddingFactor());
		return obj;
	}
	return 0;
 }



template<class TScalar, unsigned int PointDim>
typename KernelFactory<TScalar, PointDim>::KernelBaseType*
KernelFactory<TScalar, PointDim>
::CreateKernelObject(
		const VNLMatrixType& X,
		const VNLMatrixType& W,
		TScalar h)
{
	switch (m_WhichKernel)
	{
	case Exact:
		return new ExactKernel<TScalar, PointDim>(X, W, h);
#ifdef USE_CUDA
		return new CUDAExactKernel<TScalar, PointDim>(X, W, h);
#endif
	case FGT:
		return new FGTKernel<TScalar, PointDim>(X, W, h);
	case P3M:
		typedef P3MKernel<TScalar, PointDim> P3MKernelType;
		P3MKernelType* obj = new P3MKernelType(X, W, h);
		if (m_DataDomain.size())
			obj->SetDataDomain(this->GetDataDomain());
		obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
		obj->SetPaddingFactor(this->GetPaddingFactor());
		return obj;
	}
	return 0;
}



#endif /* _KernelFactory_txx */
