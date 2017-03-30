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

#ifndef _KernelFactory_h
#define _KernelFactory_h

#include "itkSimpleFastMutexLock.h"

#include "ExactKernel.h"
#include "FGTKernel.h"


/**
 *	\brief      A kernel factory.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    The KernelFactory class enables to instantiate objects whose type is derived from an
 *                  abstract type (it is the principle of the factory method pattern). On top of that, this
 *                  class implements the singleton pattern, which means that only one kernel can be used during code execution.
 */
template <class TScalar, unsigned int PointDim>
class KernelFactory
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Possible type of kernels.
	typedef enum
	{
		Exact,      /*!< Kernel with exact computation (see ExactKernel). */
#ifdef USE_CUDA
		CUDAExact,  /*!< Kernel with exact CUDA computation (see CUDAExactKernel). */
#endif
		FGT,        /*!< Kernel with Fast Gauss Transform computation (see FGTKernel). */
		P3M         /*!< Kernel with linearly spaced grid computation (see P3MKernel). */
	} KernelEnumType;

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> KernelBaseType;

	/// Vector type
	typedef typename KernelBaseType::VNLVectorType VNLVectorType;
	/// Matrix type.
	typedef typename KernelBaseType::VNLMatrixType VNLMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the type of the kernel.
	KernelEnumType GetKernelCode() const { return m_WhichKernel; }
	/// Sets the type of the kernel to ExactKernel.
	void UseExactKernel() { m_Mutex.Lock(); m_WhichKernel = Exact; m_Mutex.Unlock();}
#ifdef USE_CUDA
	/// Sets the type of the kernel to CUDAExactKernel.
	void UseCUDAExactKernel() { m_Mutex.Lock(); m_WhichKernel = CUDAExact; m_Mutex.Unlock();}
#endif
	/// Sets the type of the kernel to FGTKernel.
	void UseFGTKernel() { m_Mutex.Lock(); m_WhichKernel = FGT; m_Mutex.Unlock();}
	/// Sets the type of the kernel to P3MKernel.
	void UseP3MKernel() { m_Mutex.Lock(); m_WhichKernel = P3M; m_Mutex.Unlock();}

	/// Returns the data domain.
	VNLMatrixType GetDataDomain() const { return m_DataDomain; }
	/// Sets the data domain to \e DD.
	void SetDataDomain(VNLMatrixType DD);

	/// Get the spacing of the lattice for grid optimization (p3m).
	TScalar GetWorkingSpacingRatio() const { return m_WorkingSpacingRatio; }
	/// Set the spacing of the lattice for grid optimization (p3m).
	void SetWorkingSpacingRatio(TScalar d) { m_Mutex.Lock(); m_WorkingSpacingRatio = d; m_Mutex.Unlock(); }

	/// Returns the padding factor.
	TScalar GetPaddingFactor() const { return m_PaddingFactor; }
	/// Sets the padding factor to \e d.
	void SetPaddingFactor(TScalar d) { m_Mutex.Lock(); m_PaddingFactor = d; m_Mutex.Unlock(); }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Instantiates an object of KernelFactory type with the Singleton strategy.
	static KernelFactory<TScalar, PointDim>* Instantiate();

	/// Returns the instance of the object, NULL in case of error.
	KernelBaseType* CreateKernelObject();

	/// Returns the instance of the object, NULL in case of error.
	KernelBaseType* CreateKernelObject(
			const VNLMatrixType& X,
			const VNLMatrixType& W,
			TScalar h);

// NOTE: make sure to assign the target image when using P3M
//   void SetWorkingImage(ImageType* img);
// ImagePointer GetWorkingImage() const { return m_WorkingImage; }



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	KernelFactory();

	~KernelFactory();

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Type of the kernel.
	KernelEnumType m_WhichKernel;

	// Information of image domain for heuristically determining parameters
	// ImagePointer m_WorkingImage;

	/// See AbstractDeformations::m_DataDomain for details.
	VNLMatrixType m_DataDomain;

	/// See P3MKernel::m_WorkingSpacingRatio for details.
	TScalar m_WorkingSpacingRatio;

	/// See P3MKernel::m_PaddingFactor for details.
	TScalar m_PaddingFactor;



private:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Object used to perform mutex (important for multithreaded programming).
	static itk::SimpleFastMutexLock m_Mutex;

	/// Boolean which enables to instantiate once a kernel factory.
	static bool m_IsInstantiated;

	/// The unique kernel factory.
	static KernelFactory* m_SingletonInstance;


}; /* class KernelFactory */


#ifndef MU_MANUAL_INSTANTIATION
#include "KernelFactory.txx"
#endif


#endif /* _KernelFactory_h */
