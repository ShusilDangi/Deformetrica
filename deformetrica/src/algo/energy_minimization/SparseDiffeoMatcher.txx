/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoMatcher_txx
#define _SparseDiffeoMatcher_txx

#include "KernelFactory.h"

#include "Diffeos.h"
#include "BackwardGradientPropagatorUniversal.h"

#include "vnl/algo/vnl_qr.h"

#include "vnl/vnl_math.h"

#include "writeMatrixDLM.txx"

#include <iostream>
#include <iomanip>
#include <sstream>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
SparseDiffeoMatcher<TScalar, Dimension>
::SparseDiffeoMatcher()
 {

	this->SetMaxIterations(10);
	m_MaxLineSearchIterations = 10;

	m_AdaptiveExpand = 1.2;
	m_AdaptiveShrink = 0.5;
	m_AdaptiveTolerance = 1e-4;

	m_InitialStepMultiplier = 10.0;

	m_SparsityPrior = 0.0;
	m_OptimizationMethod = null;

	m_freezeCP = false;
	m_NumberOfCPs = 0;

 }



template <class TScalar, unsigned int Dimension>
SparseDiffeoMatcher<TScalar, Dimension>
::~SparseDiffeoMatcher()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::SetMaxIterations(unsigned int n)
 {
	m_MaxIterations = n;
	m_ValuesHistory.resize(n+1);
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
typename SparseDiffeoMatcher<TScalar, Dimension>::DeformableMultiObjectType*
SparseDiffeoMatcher<TScalar, Dimension>
::GetDeformedObject()
 {

	m_Def->SetStartPositions(m_ControlPoints);
	m_Def->SetStartMomentas(m_InitialMomentas);
	m_Def->Update();

	if (m_Def->OutOfBox())
		std::cout << "out of box in GetDeformedImage (this should not be)" << std::endl;

	m_Source->SetDiffeos(m_Def);
	m_Source->Update();
	DeformableMultiObjectType* defSource = m_Source->GetDeformedObject();

	return defSource;

 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::WriteFlow(std::vector<std::string>& name, std::vector<std::string>& extension)
 {

	m_Def->SetStartPositions(m_ControlPoints);
	m_Def->SetStartMomentas(m_InitialMomentas);
	m_Def->Update();

	if (m_Def->OutOfBox())
		std::cout << "out of box in WriteFlow (this should not be)" << std::endl;

	m_Source->SetDiffeos(m_Def);
	m_Source->Update();

	m_Source->WriteFlow(name, extension);
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::Update()
 {

	this->InitializeControlPointsAndMomentas();

	VNLMatrixType& X = m_ControlPoints;
	VNLMatrixType& A = m_InitialMomentas;

	// run optimization method
	if ( (m_OptimizationMethod == GradientDescent) && (m_SparsityPrior > 1e-10) )
	{
		m_OptimizationMethod = ISTA;
		std::cout << "Warning: gradient descent does not work with L^1 penalty term. Optimization method switched to ISTA" << std::endl;
	}
	
	if (m_OptimizationMethod == F_ISTA)
		this->FISTA(X, A);
	else
		this->GradientDescentAndISTA(X, A);

	writeMatrixDLM<TScalar>("CP_final.txt", X);
	writeMatrixDLM<TScalar>("Mom_final.txt", A);
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::GradientDescentAndISTA(VNLMatrixType& X, VNLMatrixType& A)
 {

	VNLMatrixType Xnew;
	VNLMatrixType Anew;

	TScalar step;

	m_ValuesHistory[0] = this->ComputeFunctional(X, A);
	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();

	unsigned int iterRef = 0;

	m_ValuesHistory[0]->PrintIter(0);

	for (unsigned int iter = 0; iter < m_MaxIterations; iter++)
	{

		this->ComputeGradient(X, A);

		if (iter == 0)
		{
			TScalar maxGrad = 1e-20;
			for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			{
				TScalar g = 0;
				g = m_GradPos.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
				g = m_GradMom.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
			}
			TScalar initStep = m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
			step = initStep;
		}

		std::vector<int> activeCPMask(m_NumberOfCPs, 1);
		std::vector<int> activeCPMask_new(m_NumberOfCPs, 1);
		bool foundMin = false;

		for (unsigned int li = 0; li < m_MaxLineSearchIterations; li++)
		{
			std::cout << "stepsize = " << step << std::endl;

			this->GradientDescentStep(Xnew, X, Anew, A, activeCPMask_new, step);

			m_ValuesHistory[iter+1] = this->ComputeFunctional(Xnew, Anew);
			TScalar Q = lsqRef - m_ValuesHistory[iter+1]->GetTotalL2Value();
			if (m_OptimizationMethod == ISTA)
					Q += this->_QdiffTerm(Xnew, X, Anew, A, step);

			if (Q >= 0) //(J1 < Jcurr) if no sparsity
			{
				foundMin = true;
				break;
			}
			else
				step *= m_AdaptiveShrink;

		} // for li

		if (foundMin)
		{
			X = Xnew;
			A = Anew;
			activeCPMask = activeCPMask_new;
			step *= m_AdaptiveExpand;
			lsqRef = m_ValuesHistory[iter+1]->GetTotalL2Value();
		}

		if (!foundMin)
		{
			// Loop terminated without finding smaller functional
			std::cout << " number of loops exceeded " << std::endl;
			break;
		}

		if ( (m_ValuesHistory[iter]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue())
				<
				m_AdaptiveTolerance*(m_ValuesHistory[iterRef]->GetTotalValue()-m_ValuesHistory[iter+1]->GetTotalValue()) )
			break;

		m_ValuesHistory[iter+1]->PrintIter(iter+1);

		int NbActiveCP = 0;
		for (int i = 0; i < m_NumberOfCPs; i++)
			NbActiveCP += (activeCPMask[i] != 0);
		std::cout << NbActiveCP << " active control points" << std::endl;


	} // for iter

 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::FISTA(VNLMatrixType& X, VNLMatrixType& A)
 {
	VNLMatrixType Xtest;
	VNLMatrixType Atest;

	TScalar step;

	m_ValuesHistory[0] = this->ComputeFunctional(X, A);
	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();

	VNLMatrixType Xprev = X;
	VNLMatrixType Aprev = A;

	std::vector<int> activeCPMask(m_NumberOfCPs, 1);
	std::vector<int> activeCPMask_test(m_NumberOfCPs, 1);

	unsigned int freezeDirectionCounter = 0;
	unsigned int iterRef = 0;

	m_ValuesHistory[0]->PrintIter(0);

	TScalar tau = 1.0;
	for (unsigned int iter = 0; iter < m_MaxIterations; iter++)
	{

		// Compute gradient
		this->ComputeGradient(X, A);

		// Determine the initial stepsizes
		if (iter == 0)
		{
			TScalar maxGrad = 1e-20;
			for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			{
				TScalar g = 0;
				g = m_GradPos.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
				g = m_GradMom.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
			}
			TScalar initStep = m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
			step = initStep;
		}

		// Start line search
		bool minimtest = false;
		unsigned int lineIter = 0;
		for (; lineIter < m_MaxLineSearchIterations; lineIter++)
		{
			std::cout << "stepsize = " << step << std::endl;

			this->GradientDescentStep(Xtest, X, Atest, A, activeCPMask_test, step);

			m_ValuesHistory[iter+1] = this->ComputeFunctional(Xtest, Atest);

			TScalar Q = lsqRef - m_ValuesHistory[iter+1]->GetTotalL2Value()	+ this->_QdiffTerm(
					Xtest, X, Atest, A, step);

			if (Q >= 0)
			{
				minimtest = true;
				break;
			}
			else // update failed, continue line search from X - stepX*gradX, A - stepA*gradA, T - stepT*gradT
				step *= m_AdaptiveShrink;
		} 

		if (minimtest)
		{
			TScalar tau_next = 1.0;

			if ( (lineIter == 0)
					&& ((freezeDirectionCounter == 0) || (freezeDirectionCounter > NbIterFreeze)) )
			{
				// stepsize increased
				tau_next = (1 + sqrt(1.0 + 4.0*tau*tau / m_AdaptiveExpand)) / 2.0;
				if (freezeDirectionCounter > NbIterFreeze)
					freezeDirectionCounter = 0;
			}
			else // stepsize decreased
			{
				tau_next = (1 + sqrt(1.0 + 4.0*tau*tau)) / 2.0;
				freezeDirectionCounter++;
			}

			TScalar tau_scale = (tau-1.0) / tau_next;
			X = Xtest + (Xtest - Xprev) * tau_scale;
			A = Atest + (Atest - Aprev) * tau_scale;
			activeCPMask = activeCPMask_test;

			Xprev = Xtest;
			Aprev = Atest;

			tau = tau_next;

			// Check if we're in increasing trend
			if ((freezeDirectionCounter == 0) || (freezeDirectionCounter > NbIterFreeze))
				step *= m_AdaptiveExpand;

			// Update the reference L2 part of the functional
			FunctionalValuesType* values = this->ComputeFunctional(X, A);
			if ( values->IsOutOfBox() )
			{
				std::cerr << "out of box: needs to restart FISTA from this point." << std::endl;
				X = Xtest; A = Atest;
				lsqRef = m_ValuesHistory[iter+1]->GetTotalL2Value();
			}
			else
				lsqRef = values->GetTotalL2Value();

			// these values are not necessarily decreasing, whereas Functional(Xprev) is.
			m_ValuesHistory[iter+1]->PrintIter(iter+1);

			int NbActiveCP = 0;
			for (int i = 0; i < m_NumberOfCPs; i++)
				NbActiveCP += (activeCPMask[i] != 0);
			std::cout << NbActiveCP << " active control points" << std::endl;
		}

		if (!minimtest) // Inner loop terminated without finding smaller functional
			break;

		TScalar deltaF_curr = m_ValuesHistory[iter]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue();
		TScalar deltaF_ref = m_ValuesHistory[iterRef]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue();
		if (fabs(deltaF_curr) < m_AdaptiveTolerance*fabs(deltaF_ref))
		{
			std::cout << "Tolerance BREAK" << std::endl;
			std::cout << "FINAL VALUES: ";
			m_ValuesHistory[iter+1]->PrintIter(iter+1);
			break;
		}
	} // end iter

	X = this->_maskMatrix(X, activeCPMask);
	A = this->_maskMatrix(A, activeCPMask);

 }



template <class TScalar, unsigned int Dimension>
typename SparseDiffeoMatcher<TScalar, Dimension>::VNLMatrixType
SparseDiffeoMatcher<TScalar, Dimension>
::_maskMatrix(const VNLMatrixType& M, const std::vector<int>& mask)
 {
	unsigned int numElements = 0;
	for (unsigned int i = 0; i < mask.size(); i++)
		if (mask[i] != 0)
			numElements++;

	if (numElements == M.rows())
		return M;

	VNLMatrixType Z(numElements, M.columns());

	unsigned int r = 0;
	for (unsigned int i = 0; i < M.rows(); i++)
	{
		if (mask[i] != 0)
			Z.set_row(r++, M.get_row(i));
	}

	return Z;
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::GradientDescentStep(
		VNLMatrixType& Xtest, const VNLMatrixType& X,
		VNLMatrixType& Atest, const VNLMatrixType& A,
		std::vector<int>& activeCPMask,
		TScalar step)
{
	if (!m_freezeCP)
		Xtest = X - m_GradPos*step;
	else
		Xtest = X;


	if (m_SparsityPrior == 0)
		Atest = A - m_GradMom*step;
	else
	{
		activeCPMask.assign(m_NumberOfCPs, 0); // gives for each CP if it carries a non-zero momenta
		Atest = this->SoftThresholdUpdate(A, m_GradMom, step, activeCPMask);
	}
}



template <class TScalar, unsigned int Dimension>
typename SparseDiffeoMatcher<TScalar, Dimension>::VNLMatrixType
SparseDiffeoMatcher<TScalar, Dimension>
::SoftThresholdUpdate(const VNLMatrixType& X, VNLMatrixType& gradX, TScalar step, std::vector<int>& mask)
 {
	VNLMatrixType Xnew = X - step*gradX;
	TScalar thres = step * m_SparsityPrior;

	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
	{
		TScalar ai = Xnew.get_row(i).magnitude();
		TScalar ai_thres = ai;
		if (ai > thres)
			ai_thres -= thres;
		else if (ai < -thres)
			ai_thres += thres;
		else
			ai_thres = 0;

		mask[i] += (ai_thres != 0);

		if (ai != 0.0)
			Xnew.set_row(i, Xnew.get_row(i) * (ai_thres / ai));
	}

	return Xnew;	
 }



template <class TScalar, unsigned int Dimension>
TScalar
SparseDiffeoMatcher<TScalar, Dimension>
::_QdiffTerm(
		const VNLMatrixType& Xtest, const VNLMatrixType& X,
		const VNLMatrixType& Atest, const VNLMatrixType& A,
		TScalar step)
{
	VNLMatrixType diffX = Xtest - X;
	VNLMatrixType diffA = Atest - A;

	// first term of the form (xtest - x)^t * Grad
	TScalar termVal1 = 0;
	// second term of the form ||xtest - x||
	TScalar termVal2 = 0;

	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
	{
		VNLVectorType diffXi = diffX.get_row(i);
		VNLVectorType diffAi = diffA.get_row(i);
		termVal1 += dot_product(diffXi, m_GradPos.get_row(i));
		termVal1 += dot_product(diffAi, m_GradMom.get_row(i));
		termVal2 += dot_product(diffXi, diffXi);
		termVal2 += dot_product(diffAi, diffAi);
	}

	// sum all this
	TScalar val = termVal1 + termVal2 / (2.0 * step);

	return val;
}



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::InitializeControlPointsAndMomentas()
 {
	if (m_NumberOfCPs == 0)
	{
		// create initial control points
		this->GenerateInitialControlPoints();

		m_NumberOfCPs = m_ControlPoints.rows();
		std::cout << "Generated " << m_NumberOfCPs << " control points" << std::endl;
		std::cout << "Momentas (re)set to zero" << std::endl;

		m_InitialMomentas.set_size(m_NumberOfCPs, Dimension);
		m_InitialMomentas.fill(0);
	}
	else 
	{	
		std::cout << "Using predefined set of " << m_NumberOfCPs << " control points" << std::endl;

		if ( m_InitialMomentas.rows() == m_NumberOfCPs )
		{
			std::cout << "Using predefined set of momenta" << std::endl;
		}
		else
		{
			if ( m_InitialMomentas.rows() > 0 )
				std::cout << "Warning: number of initial momenta and control points mismatch: initial momenta set to zero" << std::endl;

			m_InitialMomentas.set_size(m_NumberOfCPs, Dimension);
			m_InitialMomentas.fill(0);
		}
	}
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::GenerateInitialControlPoints()
 {
	VNLMatrixType DD = m_Def->GetDataDomain();
	VNLVectorType Xmin = DD.get_column(0);
	VNLVectorType Xmax = DD.get_column(1);

	std::vector<VNLVectorType> pointList;
	VNLVectorType v(Dimension);
	switch (Dimension)
	{
	case 2:
	{
		for (TScalar x = Xmin[0]; x <= Xmax[0]; x += m_InitialCPSpacing)
			for (TScalar y = Xmin[1]; y <= Xmax[1]; y += m_InitialCPSpacing)
			{
				v[0] = x;
				v[1] = y;
				pointList.push_back(v);
			}
		break;
	}
	case 3:
	{
		for (TScalar x = Xmin[0]; x <= Xmax[0]; x += m_InitialCPSpacing)
			for (TScalar y = Xmin[1]; y <= Xmax[1]; y += m_InitialCPSpacing)
				for (TScalar z = Xmin[2]; z <= Xmax[2]; z += m_InitialCPSpacing)
				{
					v[0] = x;
					v[1] = y;
					v[2] = z;
					pointList.push_back(v);
				}
		break;
	}
	default:
		throw std::runtime_error("GenerateInitialCP not implemented in Dimensions other than 2 and 3");
		break;
	}

	unsigned int NumCPs = pointList.size();

	m_ControlPoints.set_size(NumCPs, Dimension);
	for (unsigned int i = 0; i < NumCPs; i++)
		m_ControlPoints.set_row(i, pointList[i]);

 }



template <class TScalar, unsigned int Dimension>
typename SparseDiffeoMatcher<TScalar, Dimension>::FunctionalValuesType*
SparseDiffeoMatcher<TScalar, Dimension>
::ComputeFunctional(VNLMatrixType& X0, VNLMatrixType& Mom0)
 {

	FunctionalValuesType* funcValues = new FunctionalValuesType(m_NumberOfObjects);

	m_Def->SetStartPositions(X0);
	m_Def->SetStartMomentas(Mom0);
	m_Def->Update();

	if (m_Def->OutOfBox())
	{
		funcValues->SetOutOfBox();
		return funcValues;
	}

	m_Source->SetDiffeos(m_Def);
	m_Source->Update();
	if (m_Source->OutOfBox())
	{
		funcValues->SetOutOfBox();
		return funcValues;
	}

	std::vector<TScalar> DataTerm = m_Source->ComputeMatch(m_Target);

	DeformableObjectList sourceList = m_Source->GetObjectList();
	for (int i = 0; i < m_NumberOfObjects; i++)
		DataTerm[i] /=  2.0 * sourceList[i]->GetDataSigmaSquared();

	// compute the regularity term: the squared norm of the initial velocity speed (halfed)
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject();
	momKernelObj->SetKernelWidth(m_Def->GetKernelWidth());
	momKernelObj->SetSources(X0);
	momKernelObj->SetWeights(Mom0);
	VNLMatrixType kMom = momKernelObj->Convolve(X0);

	TScalar Reg = 0.0;
	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
		Reg += 0.5 * dot_product(kMom.get_row(i), Mom0.get_row(i));

	delete momKernelObj;

	// compute the sparsity prior
	TScalar SparsityPrior = 0.0;
	for (int i = 0; i < m_NumberOfCPs; i++)
		SparsityPrior += Mom0.get_row(i).magnitude();

	SparsityPrior *= m_SparsityPrior;

	funcValues->SetDataTerm(DataTerm);
	funcValues->SetRegularity(Reg);
	funcValues->SetSparsity(SparsityPrior);
	funcValues->Update();

	return funcValues;
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoMatcher<TScalar, Dimension>
::ComputeGradient(VNLMatrixType& X0, VNLMatrixType& Mom0)
 {

	// Compute forward trajectory
	m_Def->SetStartPositions(X0);
	m_Def->SetStartMomentas(Mom0);
	m_Def->Update();

	if (m_Def->OutOfBox())
		std::cout << "out of box in Compute Gradient (this should not be...)" << std::endl;

	VNLMatrixList posT = m_Def->GetTrajectoryPositions();
	VNLMatrixList momT = m_Def->GetTrajectoryMomentas();

	// Compute eta(t)
	m_Source->SetDiffeos(m_Def);
	m_Source->Update();

	VNLMatrixList InitialConditions = m_Source->ComputeMatchGradient(m_Target);

	DeformableObjectList sourceList = m_Source->GetObjectList();	
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if (sourceList[i]->IsOfLandmarkKind())
			InitialConditions[i] /= 2.0 * sourceList[i]->GetDataSigmaSquared();
	}

	VNLMatrixList VectorsT;
	VNLMatrixList PointsT;

	m_Source->TransportAlongGeodesic(InitialConditions,VectorsT,PointsT);

	typedef BackwardGradientPropagatorUniversal<TScalar, Dimension> BackwardEstimator;
	BackwardEstimator* gradbwd = new BackwardEstimator();
	gradbwd->SetPointsTrajectory(PointsT);
	gradbwd->SetVectorsTrajectory(VectorsT);
	gradbwd->SetDiffeos(m_Def);
	gradbwd->Update();

	VNLMatrixType dPos = gradbwd->GetGradientPosAt(0);
	VNLMatrixType dMom = gradbwd->GetGradientMomAt(0);

	delete gradbwd;

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject();
	momKernelObj->SetKernelWidth(m_Def->GetKernelWidth());
	momKernelObj->SetSources(X0);
	momKernelObj->SetWeights(Mom0);

	VNLMatrixType kMom = momKernelObj->Convolve(X0);

	dMom += kMom;

	dPos = dPos + momKernelObj->ConvolveGradient(X0, Mom0);


	delete momKernelObj;

	m_GradPos = dPos;
	m_GradMom = dMom;

 }


#endif /* _SparseDiffeoMatcher_txx */
