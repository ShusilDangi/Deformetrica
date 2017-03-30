/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Diffeos_txx
#define _Diffeos_txx

#include "Diffeos.h"

#include "KernelFactory.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
Diffeos<TScalar, Dimension>
::Diffeos() : Superclass()
 {
	this->SetDiffeosType();
 }



template <class TScalar, unsigned int Dimension>
Diffeos<TScalar, Dimension>
::Diffeos(const Diffeos& other) : Superclass(other)
 {
	this->SetDiffeosType();
 }



template <class TScalar, unsigned int Dimension>
Diffeos<TScalar, Dimension>
::~Diffeos() {}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::ReverseFlow()
 {
	int midPoint = static_cast<int>( 0.5*(this->GetNumberOfTimePoints() + 1) );

	for (int t = 0; t < midPoint; t++)
	{
		VNLMatrixType auxP = Superclass::m_PositionsT[t];
		VNLMatrixType auxM = - Superclass::m_MomentasT[t];

		Superclass::m_PositionsT[t] = Superclass::m_PositionsT[this->GetNumberOfTimePoints() - t - 1];
		Superclass::m_MomentasT[t] = - Superclass::m_MomentasT[this->GetNumberOfTimePoints() - t- 1];

		Superclass::m_PositionsT[this->GetNumberOfTimePoints() - t- 1] = auxP;
		Superclass::m_MomentasT[this->GetNumberOfTimePoints() - t - 1] = auxM;
	}

 }



template <class TScalar, unsigned int Dimension>
void
Diffeos<TScalar, Dimension>
::Shoot()
 {
	unsigned int numCP = Superclass::m_StartPositions.rows();
	unsigned int numTimePoints = this->GetNumberOfTimePoints();

	VNLMatrixList& outPos = Superclass::m_PositionsT;
	VNLMatrixList& outMoms = Superclass::m_MomentasT;

	outPos.resize(numTimePoints);
	outMoms.resize(numTimePoints);
	for (unsigned int t = 0; t < numTimePoints; t++)
	{
		outPos[t] = Superclass::m_StartPositions;
		outMoms[t].set_size(numCP, Dimension);
		outMoms[t].fill(0);
	}
	outPos[0] = Superclass::m_StartPositions;
	outMoms[0] = Superclass::m_StartMomentas;

	// Special case: nearly zero momentas yield no motion
	if (outMoms[0].frobenius_norm() < 1e-20)
		return;

	TScalar timeStep = 1.0 / (numTimePoints-1);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	KernelType* kernelObj = kFactory->CreateKernelObject();
	kernelObj->SetKernelWidth(this->GetKernelWidth());

	for (unsigned int t = 0; t < (numTimePoints-1); t++)
	{
		kernelObj->SetSources(outPos[t]);
		kernelObj->SetWeights(outMoms[t]);

		typename KernelType::VNLMatrixType dPos = kernelObj->Convolve(outPos[t]);
		typename KernelType::VNLMatrixType dMom = kernelObj->ConvolveGradient(outPos[t], outMoms[t]);

		outPos[t+1] = outPos[t] + dPos * timeStep;
		outMoms[t+1] = outMoms[t] - dMom * timeStep;

//		// Heun's method
//		if (m_UseImprovedEuler)
//		{
//			kernelObj->SetSources(outPos[t+1]);
//			kernelObj->SetWeights(outMoms[t+1]);
//
//			kGradMom = kernelObj->ConvolveGradient(outPos[t+1]);
//
//			typename KernelType::VNLMatrixType dPos2 =
//					kernelObj->Convolve(outPos[t+1]);
//
//			typename KernelType::VNLMatrixType dMom2(numCP, Dimension, 0);
//			for (unsigned int i = 0; i < numCP; i++)
//				dMom2.set_row(i,
//						(kGradMom[i].transpose() * outMoms[t+1].get_row(i)) );
//
//			outPos[t+1] = outPos[t] + (dPos + dPos2) * (timeStep * 0.5);
//			outMoms[t+1] = outMoms[t] - (dMom + dMom2) * (timeStep * 0.5);
//		}

		Superclass::m_OutOfBox = this->CheckBoundingBox(outPos,t+1);
		if (this->OutOfBox())
		{
			std::cerr << "Geodesic shooting: out of box at time t = " << t+1 << std::endl;
			break;
		}
	}

	delete kernelObj;
 }


#endif /* _Diffeos_txx */
