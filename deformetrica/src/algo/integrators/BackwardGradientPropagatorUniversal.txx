/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _BackwardGradientPropagatorUniversal_txx
#define _BackwardGradientPropagatorUniversal_txx

#include "SimpleTimer.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
BackwardGradientPropagatorUniversal<TScalar, Dimension>
::BackwardGradientPropagatorUniversal()
 {
	m_KernelObj1 = 0;
	m_KernelObj2 = 0;
	m_KernelObj3 = 0;
 }



template <class TScalar, unsigned int Dimension>
BackwardGradientPropagatorUniversal<TScalar, Dimension>
::~BackwardGradientPropagatorUniversal()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
BackwardGradientPropagatorUniversal<TScalar, Dimension>
::Update()
 {
	long numTimePoints = m_Def->GetNumberOfTimePoints();

	VNLMatrixList PositionsT = m_Def->GetTrajectoryPositions();
	long numCP = PositionsT[0].rows();

	// Initialization: Xi[T-1] = 0
	VNLMatrixType zeroM(numCP, Dimension, 0);

	m_XiPosT.resize(numTimePoints);
	m_XiMomT.resize(numTimePoints);
	for (long t = 0; t < numTimePoints; t++)
	{
		m_XiPosT[t] = zeroM;
		m_XiMomT[t] = zeroM;
	}

	// Set up kernel objects
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	TScalar h = m_Def->GetKernelWidth();
	m_KernelObj1 = kFactory->CreateKernelObject();
	m_KernelObj1->SetKernelWidth(h);
	m_KernelObj2 = kFactory->CreateKernelObject();
	m_KernelObj2->SetKernelWidth(h);
	m_KernelObj3 = kFactory->CreateKernelObject();
	m_KernelObj3->SetKernelWidth(h);

	// Backward propagation
	TScalar dt = 1.0 / (numTimePoints - 1.0);

	for (long t = numTimePoints-1; t >= 1; t--)
	{
		VNLMatrixType dPos;
		VNLMatrixType dMom;
		this->ComputeUpdateAt(t, dPos, dMom);

		m_XiPosT[t-1] = m_XiPosT[t] + dPos*dt;
		m_XiMomT[t-1] = m_XiMomT[t] + dMom*dt;

		// Heun's method
		if (m_Def->ImprovedEuler())
		{
			VNLMatrixType dPos2;
			VNLMatrixType dMom2;
			this->ComputeUpdateAt(t-1, dPos2, dMom2);

			m_XiPosT[t-1] = m_XiPosT[t] + (dPos + dPos2) * (dt * 0.5);
			m_XiMomT[t-1] = m_XiMomT[t] + (dMom + dMom2) * (dt * 0.5);
		}

	} // for t

	delete m_KernelObj1;
	delete m_KernelObj2;
	delete m_KernelObj3;

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
BackwardGradientPropagatorUniversal<TScalar, Dimension>
::ComputeUpdateAt(unsigned int s, VNLMatrixType& dPos, VNLMatrixType &dMom)
 {

	VNLMatrixList PositionsT = m_Def->GetTrajectoryPositions();
	VNLMatrixList MomentasT = m_Def->GetTrajectoryMomentas();

	long numCP = PositionsT[0].rows();

	KernelType* momXiPosKernelObj = m_KernelObj1;
	KernelType* etaKernelObj = m_KernelObj2;
	KernelType* tmpKernelObj = m_KernelObj3;

	// Term 1
	etaKernelObj->SetSources(m_PointsT[s]);
	etaKernelObj->SetWeights(m_VectorsT[s]);

	VNLMatrixType dXi1 = etaKernelObj->ConvolveGradient(PositionsT[s], MomentasT[s]);

	// Term 2
	VNLMatrixType dXi2 = etaKernelObj->Convolve(PositionsT[s]);

	// Term 3
	VNLMatrixType AXiPos(numCP, Dimension*2, 0);
	AXiPos.set_columns(0, MomentasT[s]);
	AXiPos.set_columns(Dimension, m_XiPosT[s]);

	momXiPosKernelObj->SetSources(PositionsT[s]);
	momXiPosKernelObj->SetWeights(AXiPos);

	VNLMatrixType kAXiPos = momXiPosKernelObj->Convolve(PositionsT[s]);

	VNLMatrixList gradAXiPos =
			momXiPosKernelObj->ConvolveGradient(PositionsT[s]);

	VNLMatrixType dXi3(numCP, Dimension, 0);
	for (unsigned int i = 0; i < numCP; i++)
	{
		VNLMatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
		VNLMatrixType gradXiPos_i = gradAXiPos[i].get_n_rows(Dimension, Dimension);
		dXi3.set_row(i,
				gradMom_i.transpose() * m_XiPosT[s].get_row(i)
				+ gradXiPos_i.transpose() * MomentasT[s].get_row(i));
	}

	// Term 5
	VNLMatrixType dXi5 = kAXiPos.get_n_columns(Dimension, Dimension);

	// Term 6
	VNLMatrixType dXi6(numCP, Dimension, 0);

	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		VNLMatrixType W(numCP, Dimension, 0.0);
		for (unsigned int i = 0; i < numCP; i++)
			W.set_row(i, m_XiMomT[s](i,dim) * MomentasT[s].get_row(i));

		tmpKernelObj->SetSources(PositionsT[s]);
		tmpKernelObj->SetWeights(W);
		VNLMatrixType tmpgrad = tmpKernelObj->ConvolveGradient(PositionsT[s], dim);

		dXi6 += tmpgrad;
	}

	for (unsigned int i = 0; i < numCP; i++)
	{
		VNLMatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
		dXi6.set_row(i, dXi6.get_row(i) -
				gradMom_i * m_XiMomT[s].get_row(i));
	}
/*
	// Term 4
	VNLMatrixType dXi4(numCP, Dimension, 0);

	tmpKernelObj->SetSources(PositionsT[s]);
	tmpKernelObj->SetWeights(MomentasT[s]);
	std::vector< std::vector<VNLMatrixType> > HessA = tmpKernelObj->ConvolveHessian(PositionsT[s]);	
	for (unsigned int i = 0; i < numCP; i++)
	{
		VNLVectorType XiMomi = m_XiMomT[s].get_row(i);	
		VNLMatrixType Auxi(Dimension, Dimension, 0.0);
		for (unsigned int dim = 0; dim < Dimension; dim++)
			Auxi -= HessA[i][dim] * MomentasT[s](i,dim);

		dXi4.set_row(i, Auxi * XiMomi);
	}
	// dXi4 second term
	std::vector< VNLMatrixType > H(numCP);
	VNLMatrixType Zeros(Dimension, Dimension, 0.0);
	for (int i = 0; i < numCP; i++)
		H[i] = Zeros;

	for (unsigned r = 0; r < Dimension; r++)
	{
		VNLMatrixType Wd(numCP,Dimension,0.0);
		for (int i = 0; i < numCP; i++)
			Wd.set_row(i, MomentasT[s].get_row(i) * m_XiMomT[s](i,r) );

		tmpKernelObj->SetWeights(Wd);
		VNLMatrixType Hrr = tmpKernelObj->ConvolveHessian(PositionsT[s],r,r);
		for (int i = 0; i < numCP; i ++)
			for (unsigned int q = 0; q < Dimension; q++)
				H[i](r,q) += Hrr(i,q);

		for (unsigned q = r+1; q < Dimension; q++)
		{
			VNLMatrixType W(numCP, 2*Dimension, 0.0);
			for (int i = 0; i < numCP; i++)
				for (unsigned int d = 0; d < Dimension; d++)
				{
					W(i,d) = Wd(i,d);
					W(i, d + Dimension) = MomentasT[s](i,d) * m_XiMomT[s](i,q);
				}

			tmpKernelObj->SetWeights(W);
			VNLMatrixType Hrq = tmpKernelObj->ConvolveHessian(PositionsT[s],r,q);
			for (int i = 0; i < numCP; i ++)
				for (unsigned int d = 0; d < Dimension; d++)
				{
					H[i](q,d) += Hrq(i,d);
					H[i](r,d) += Hrq(i,d + Dimension);
				}
		}		
	}

	for (int i = 0; i < numCP; i++)
		dXi4.set_row(i, dXi4.get_row(i) + H[i] * MomentasT[s].get_row(i));
*/
	// Term 4
	tmpKernelObj->SetSources(PositionsT[s]);
	tmpKernelObj->SetWeights(MomentasT[s]);

	VNLMatrixType dXi4 = tmpKernelObj->ConvolveSpecialHessian(m_XiMomT[s]);	



	dPos = dXi1 + dXi3 + dXi4;
	dMom = dXi2 + dXi5 + dXi6;

 }



#endif /* _BackwardGradientPropagatorUniversal_txx */
