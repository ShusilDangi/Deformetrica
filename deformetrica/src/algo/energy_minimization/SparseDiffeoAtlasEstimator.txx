/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoAtlasEstimator_txx
#define _SparseDiffeoAtlasEstimator_txx

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
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::SparseDiffeoAtlasEstimator()
 {
	this->SetMaxIterations(10);
	m_MaxLineSearchIterations = 10;

	m_AdaptiveExpand = 1.2;
	m_AdaptiveShrink = 0.5;
	m_AdaptiveTolerance = 1e-4;

	m_InitialStepMultiplier = 1.0;

	m_SparsityPrior = 0.0;
	SetWeights("uniform");
	m_OptimizationMethod = null;

	m_freezeCP = false;

	m_NumberOfCPs = 0;

	m_NumberOfThreads = 1;

	m_MT_FuncValues = 0;
 }



template <class TScalar, unsigned int Dimension>
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::~SparseDiffeoAtlasEstimator()
{

}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::SetMaxIterations(unsigned int n)
 {
	m_MaxIterations = n;
	m_ValuesHistory.resize(n+1);
	// for (int i = 0; i <= n; i++)
	// 	m_ValuesHistory[i].resize(n_Values, 0.0);
 }



template <class TScalar, unsigned int Dimension>
typename SparseDiffeoAtlasEstimator<TScalar, Dimension>::DeformableMultiObjectType*
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::GetTemplate()
 {

	m_Template->SetImageAndPointData(m_TemplateData);
	return m_Template;

 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


template <class TScalar, unsigned int Dimension>
typename SparseDiffeoAtlasEstimator<TScalar, Dimension>::DeformableMultiObjectList
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::GetDeformedObjects()
 {

	DeformableMultiObjectList defTemplate(m_NumberOfSubjects);

	for (int s = 0; s < m_NumberOfSubjects; s++)
	{
		m_Def->SetStartPositions(m_ControlPoints);
		m_Def->SetStartMomentas(m_InitialMomentas[s]);
		m_Def->Update();

		if (m_Def->OutOfBox())
			std::cout << "out of box in GetDeformedImage (this should not be)" << std::endl;

		m_Template->SetDiffeos(m_Def);
		m_Template->Update();

		defTemplate[s] = m_Template->GetDeformedObject();

	}

	return defTemplate;
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::WriteFlow(std::vector<std::string>& name, std::vector<std::string>& extension)
 {

	m_Template->SetImageAndPointData(m_TemplateData);

	for (int s = 0; s < m_NumberOfSubjects; s++)
	{
		m_Def->SetStartPositions(m_ControlPoints);
		m_Def->SetStartMomentas(m_InitialMomentas[s]);
		m_Def->Update();

		if (m_Def->OutOfBox())
			std::cout << "out of box in WriteFlow (this should not be)" << std::endl;

		m_Template->SetDiffeos(m_Def);
		m_Template->Update();

		std::vector<std::string> nameSubj(m_NumberOfObjects);
		for (int i = 0; i < m_NumberOfObjects; i++)
		{
			std::ostringstream oss;
			oss << name[i] << "_to_subject_" << s;
			nameSubj[i] = oss.str();
		}

		m_Template->WriteFlow(nameSubj, extension);
	}

 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::Update()
 {

	this->InitializeControlPointsAndMomentas();


	m_TemplateData = m_Template->GetImageAndPointData();

	VNLMatrixType& X = m_ControlPoints;
	VNLMatrixList& A = m_InitialMomentas;
	VNLMatrixList& T = m_TemplateData;

	// run optimization method
	if ( (m_OptimizationMethod == GradientDescent) && (m_SparsityPrior > 1e-10) )
	{
		m_OptimizationMethod = ISTA;
		std::cout << "Warning: gradient descent does not work with L^1 penalty term. Optimization method switched to ISTA" << std::endl;
	}

	if (m_OptimizationMethod == F_ISTA)
		this->FISTA(X, A, T);
	else
		this->GradientDescentAndISTA(X, A, T);

	writeMatrixDLM<TScalar>("CP_final.txt", X);
	writeMultipleMatrixDLM<TScalar>("MOM_final.txt",A);
 }



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::GradientDescentAndISTA(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T)
 {

	VNLMatrixType Xnew;
	VNLMatrixList Anew(m_NumberOfSubjects);
	VNLMatrixList Tnew(m_NumberOfObjects);

	m_GradPos.set_size(m_NumberOfCPs, Dimension);
	m_GradMom.resize(m_NumberOfSubjects);
	m_GradTemplate_L2.resize(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_GradTemplate_L2[i].set_size(m_TemplateData[i].rows(), m_TemplateData[i].columns());

	TScalar stepXA, stepT;

	m_ValuesHistory[0] = this->ThreadedComputeFunctional(X, A, T);
	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();

	unsigned int iterRef = 0;
	m_ValuesHistory[0]->PrintIter(0);

	for (unsigned int iter = 0; iter < m_MaxIterations; iter++)
	{

		// if (!(iter % 3))
		// {
		// 	std::ostringstream oss;
		// 	oss << "CP_iter" << std::setfill('0') << std::setw(4) << iter+1
		// 			<< ".txt" << std::ends;
		// 	writeMatrixDLM<TScalar>(oss.str().c_str(), X);
		// 
		// 	oss.str("");
		// 	oss << "Mom0_iter" << std::setfill('0') << std::setw(4) << iter+1
		// 			<< ".txt" << std::ends;
		// 	writeMultipleMatrixDLM<TScalar>(oss.str().c_str(), A);
		// 
		// 	DeformableMultiObjectType* templ = this->GetTemplate();
		// 
		// 	std::vector<std::string> templateName(m_NumberOfObjects);
		// 	for (int k = 0; k < m_NumberOfObjects; k++)
		// 	{
		// 		std::ostringstream oss;
		// 		oss << "Template_object_" << k << "_iter" << std::setfill('0') << std::setw(4) << iter+1 << ".vtk" << std::ends;
		// 		templateName[k] = oss.str();
		// 	}
		// 	templ->WriteDeformedMultiObjectAt(0, templateName);
		// }

		this->ThreadedComputeGradient(X, A, T);

		if (iter == 0)
		{
			TScalar maxGrad = 1e-20;
			for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			{
				TScalar g = 0;
				g = m_GradPos.get_row(i).magnitude();
				if (g > maxGrad)
					maxGrad = g;
				for (int s = 0; s < m_NumberOfSubjects; s++)
				{
					g = m_GradMom[s].get_row(i).magnitude();
					if (g > maxGrad)
						maxGrad = g;
				}
			}
			TScalar initStepXA = m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
			TScalar initStepT = initStepXA;

			stepXA = initStepXA;
			stepT = initStepT;
		}

		std::vector<int> activeCPMask(m_NumberOfCPs, 1);
		bool foundMin = false;

		for (unsigned int li = 0; li < m_MaxLineSearchIterations; li++)
		{
			std::cout << "stepsizeXA = " << stepXA << "\tstepsizeT = " << stepT << std::endl;

			this->ThreadedGradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);

			m_ValuesHistory[iter+1] = this->ThreadedComputeFunctional(Xnew, Anew, Tnew);
			std::cout << "Regularity = " << m_ValuesHistory[iter+1]->GetTotalRegularity() << " " << "Data Term = " << m_ValuesHistory[iter+1]->GetTotalDataTerm() << std::endl;

			TScalar Q = lsqRef - m_ValuesHistory[iter+1]->GetTotalL2Value();
			if (m_OptimizationMethod == ISTA)
					Q += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);


			if (Q >= 0) //(J1 < Jcurr) if no sparsity
			{
				foundMin = true;
				break;
			}
			else
			{
				// case test 1:
				stepXA *= m_AdaptiveShrink;
				this->ThreadedGradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);

				FunctionalValuesType* valuesTest1 = this->ThreadedComputeFunctional(Xnew, Anew, Tnew);

				TScalar Q1 = lsqRef - valuesTest1->GetTotalL2Value();
				if (m_OptimizationMethod == ISTA)
						Q1 += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);

				// case test 2:
				stepXA /= m_AdaptiveShrink;
				stepT *= m_AdaptiveShrink;
				this->ThreadedGradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);

				FunctionalValuesType* valuesTest2 = this->ThreadedComputeFunctional(Xnew, Anew, Tnew);

				TScalar Q2 = lsqRef - valuesTest2->GetTotalL2Value();
				if (m_OptimizationMethod == ISTA)
						Q2 += this->_QdiffTerm(Xnew, X, Anew, A, Tnew, T, stepXA, stepT);

				if ( (Q1 >= 0) || (Q2 >= 0) ) //( (Jtest1 < Jcurr) || (Jtest2 < Jcurr) )
				{
					if ( Q1 >= Q2)
					{
						stepXA *= m_AdaptiveShrink;
						stepT /= m_AdaptiveShrink;
						this->ThreadedGradientDescentStep(Xnew, X, Anew, A, Tnew, T, activeCPMask, stepXA, stepT);
						m_ValuesHistory[iter+1] = valuesTest1->Clone();
					}
					else
					{
						m_ValuesHistory[iter+1] = valuesTest2->Clone();
					}

					foundMin = true;
					break;
				}
				else
				{
					stepXA *= m_AdaptiveShrink;
				}

				delete valuesTest1;
				delete valuesTest2;

			}
		} // for li


		if (foundMin)
		{
			X = Xnew;
			A = Anew;
			T = Tnew;
			stepXA *= m_AdaptiveExpand;
			stepT *= m_AdaptiveExpand;
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
				m_AdaptiveTolerance*(m_ValuesHistory[iterRef]->GetTotalValue() - m_ValuesHistory[iter+1]->GetTotalValue()) )
			break;

		m_ValuesHistory[iter+1]->PrintIter(iter+1);

	} // for iter
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::FISTA(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T)
 {
	VNLMatrixType Xtest;
	VNLMatrixList Atest(m_NumberOfSubjects);
	VNLMatrixList Ttest(m_NumberOfObjects);

	m_GradPos.set_size(m_NumberOfCPs, Dimension);
	m_GradMom.resize(m_NumberOfSubjects);
	m_GradTemplate_L2.resize(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_GradTemplate_L2[i].set_size(m_TemplateData[i].rows(), m_TemplateData[i].columns());

	TScalar stepXA, stepT;

	m_ValuesHistory[0] = this->ThreadedComputeFunctional(X, A, T);
	TScalar lsqRef = m_ValuesHistory[0]->GetTotalL2Value();

	VNLMatrixType Xprev = X;
	VNLMatrixList Aprev = A;
	VNLMatrixList Tprev = T;

	std::vector<int> activeCPMask(m_NumberOfCPs, 1);
	std::vector<int> activeCPMask_test(m_NumberOfCPs, 1);

	TScalar Qarray[4];
	unsigned int freezeDirectionCounter = 0;
	unsigned int iterRef = 0;

	m_ValuesHistory[0]->PrintIter(0);

	TScalar tau = 1.0;
	for (unsigned int iter = 0; iter < m_MaxIterations; iter++)
	{

		this->ThreadedComputeGradient(X, A, T);

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
				for (int s = 0; s < m_NumberOfSubjects; s++)
				{
					g = m_GradMom[s].get_row(i).magnitude();
					if (g > maxGrad)
						maxGrad = g;
				}
			}
			TScalar initStepXA = m_InitialStepMultiplier * m_ValuesHistory[0]->GetTotalL2Value() / maxGrad / maxGrad;
			TScalar initStepT = initStepXA;

			stepXA = initStepXA;
			stepT = initStepT;
		}

		// Start line search
		bool minimtest = false;
		unsigned int lineIter = 0;
		for (; lineIter < m_MaxLineSearchIterations; lineIter++)
		{
			std::cout << "stepsizeXA = " << stepXA << "\tstepsizeT = " << stepT << std::endl;

			this->ThreadedGradientDescentStep(Xtest, X, Atest, A, Ttest, T, activeCPMask_test, stepXA, stepT);

			m_ValuesHistory[iter+1] = this->ThreadedComputeFunctional(Xtest, Atest, Ttest);

			TScalar Q_XAT = lsqRef -  m_ValuesHistory[iter+1]->GetTotalL2Value() + this->_ThreadedQdiffTerm(
					Xtest, X, Atest, A, Ttest, T, stepXA, stepT);

			if (Q_XAT >= 0)
				minimtest = true;
			else
			{
				// test case 0
				stepXA *= m_AdaptiveShrink;
				std::vector<int> activeCPMask_prime1(m_NumberOfCPs, 1);
				VNLMatrixType X_prime1;
				VNLMatrixList A_prime1(m_NumberOfSubjects);
				this->ThreadedGradientDescentStep(X_prime1, X, A_prime1, A, Ttest, T, activeCPMask_prime1, stepXA, stepT);
				FunctionalValuesType* values_XAprime1 = this->ThreadedComputeFunctional(X_prime1, A_prime1, Ttest);
				TScalar Q_X1 = lsqRef - values_XAprime1->GetTotalL2Value() + this->_ThreadedQdiffTerm(
						X_prime1, X, A_prime1, A, Ttest, T, stepXA, stepT);

				// test case 1
				stepXA *= m_AdaptiveShrink;
				std::vector<int> activeCPMask_prime2(m_NumberOfCPs, 1);
				VNLMatrixType X_prime2;
				VNLMatrixList A_prime2(m_NumberOfSubjects);
				this->ThreadedGradientDescentStep(X_prime2, X, A_prime2, A, Ttest, T, activeCPMask_prime2, stepXA, stepT);
				FunctionalValuesType* values_XAprime2 = this->ThreadedComputeFunctional(X_prime2, A_prime2, Ttest);
				TScalar Q_X2 = lsqRef - values_XAprime2->GetTotalL2Value() + this->_ThreadedQdiffTerm(
						X_prime2, X, A_prime2, A, Ttest, T, stepXA, stepT);

				// test case 2
				stepXA /= (m_AdaptiveShrink * m_AdaptiveShrink);
				stepT *= m_AdaptiveShrink;
				VNLMatrixList T_prime1(m_NumberOfObjects);
				this->ThreadedGradientDescentStep(Xtest, X, Atest, A, T_prime1, T, activeCPMask_test, stepXA, stepT);
				FunctionalValuesType* values_Tprime1 = this->ThreadedComputeFunctional(Xtest, Atest, T_prime1);
				TScalar Q_T1 = lsqRef - values_Tprime1->GetTotalL2Value() + this->_ThreadedQdiffTerm(
						Xtest, X, Atest, A, T_prime1, T, stepXA, stepT);

				// test case 3
				stepT *= m_AdaptiveShrink;
				VNLMatrixList T_prime2(m_NumberOfObjects);
				this->ThreadedGradientDescentStep(Xtest, X, Atest, A, T_prime2, T, activeCPMask_test, stepXA, stepT);
				FunctionalValuesType* values_Tprime2 = this->ThreadedComputeFunctional(Xtest, Atest, T_prime2);
				TScalar Q_T2 = lsqRef - values_Tprime2->GetTotalL2Value() + this->_ThreadedQdiffTerm(
						Xtest, X, Atest, A, T_prime2, T, stepXA, stepT);

				// Select best case
				Qarray[0] = Q_X1;
				Qarray[1] = Q_X2;
				Qarray[2] = Q_T1;
				Qarray[3] = Q_T2;

				int ind_max = 0;
				TScalar Q_max = Qarray[0];
				for (int j = 1; j < 4; j++)
					if (Qarray[j] > Q_max)
					{
						ind_max = j;
						Q_max = Qarray[j];
					}

				if (Q_max > 0)
				{
					minimtest = true;
					if (ind_max == 0)
					{
						Xtest = X_prime1;
						Atest = A_prime1;
						activeCPMask_test = activeCPMask_prime1;
						m_ValuesHistory[iter+1] = values_XAprime1->Clone();

						stepXA *= m_AdaptiveShrink;
						stepT /= (m_AdaptiveShrink * m_AdaptiveShrink);
					}
					else if (ind_max == 1)
					{
						Xtest = X_prime2;
						Atest = A_prime2;
						activeCPMask_test = activeCPMask_prime2;
						m_ValuesHistory[iter+1] = values_XAprime2->Clone();
						stepXA *= (m_AdaptiveShrink * m_AdaptiveShrink);
						stepT /= (m_AdaptiveShrink * m_AdaptiveShrink);
					}
					else if (ind_max == 2)
					{
						Ttest = T_prime1;
						m_ValuesHistory[iter+1] = values_Tprime1->Clone();
						stepT /= m_AdaptiveShrink;
					}
					else
					{
						Ttest = T_prime2;
						m_ValuesHistory[iter+1] = values_Tprime2->Clone();
					}
				} // if (Q_max > 0)
				else // update failed, continue line search from X - stepX*gradX, A - stepA*gradA, T - stepT*gradT
				{
					stepXA *= m_AdaptiveShrink;
					stepT /= m_AdaptiveShrink;
				}

				delete values_XAprime1;
				delete values_XAprime2;
				delete values_Tprime1;
				delete values_Tprime2;

			} // if (Q_XAT >= 0)

			if (minimtest)
				break;
		} // end line search


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
			for (unsigned int s = 0; s < m_NumberOfSubjects; s++)
				A[s] = Atest[s] + (Atest[s] - Aprev[s]) * tau_scale;
			for (unsigned int i = 0; i < m_NumberOfObjects; i++)
				T[i] = Ttest[i] + (Ttest[i] - Tprev[i]) * tau_scale;
			activeCPMask = activeCPMask_test;

			Xprev = Xtest;
			Aprev = Atest;
			Tprev = Ttest;

			tau = tau_next;

			// Check if we're in increasing trend
			if ((freezeDirectionCounter == 0) || (freezeDirectionCounter > NbIterFreeze))
			{
				stepXA *= m_AdaptiveExpand;
				stepT *= m_AdaptiveExpand;
			}

			// Update the reference L2 part of the functional
			FunctionalValuesType* values = this->ThreadedComputeFunctional(X, A, T);
			if (values->IsOutOfBox()) //m_Def->OutOfBox() || m_Template->OutOfBox())
			{
				std::cerr << "out of box: needs to restart FISTA from this point." << std::endl;
				X = Xtest; A = Atest; T = Ttest;
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
	for (unsigned int s = 0; s < m_NumberOfSubjects; s++)
		A[s] = this->_maskMatrix(A[s], activeCPMask);


 }



template <class TScalar, unsigned int Dimension>
typename SparseDiffeoAtlasEstimator<TScalar, Dimension>::VNLMatrixType
SparseDiffeoAtlasEstimator<TScalar, Dimension>
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
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::GradientDescentStep(
		VNLMatrixType& Xtest, const VNLMatrixType& X,
		VNLMatrixList& Atest, const VNLMatrixList& A,
		VNLMatrixList& Ttest, const VNLMatrixList& T,
		std::vector<int>& activeCPMask,
		TScalar stepXA, TScalar stepT)
{
	if (!m_freezeCP)
		Xtest = X - m_GradPos*stepXA;
	else
		Xtest = X;

	if (m_SparsityPrior == 0)
	{
		activeCPMask.assign(m_NumberOfCPs, 1);
		for (int s = 0; s < m_NumberOfSubjects; s++)
			Atest[s] = A[s] - m_GradMom[s] * stepXA;
	}
	else
	{
		activeCPMask.assign(m_NumberOfCPs, 0); // gives for each CP the number of non-zero momenta that it carries
		for (int s = 0; s < m_NumberOfSubjects; s++)
			Atest[s] = this->SoftThresholdUpdate(A[s], m_GradMom[s], stepXA, activeCPMask);
	}

	for (int i = 0; i < m_NumberOfObjects; i++)
		Ttest[i] = T[i] - m_GradTemplate_Sob[i]*stepT;
}



template <class TScalar, unsigned int Dimension>
typename SparseDiffeoAtlasEstimator<TScalar, Dimension>::VNLMatrixType
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::SoftThresholdUpdate(const VNLMatrixType& X, VNLMatrixType& gradX, TScalar step, std::vector<int>& mask)
 {
	VNLMatrixType Xnew = X - step*gradX;
	TScalar thres = step * m_SparsityPrior;

	for (unsigned int i = 0; i < m_NumberOfCPs; i++)
	{
		TScalar ai = Xnew.get_row(i).magnitude();
		TScalar ai_thres = ai;
		// soft-thresholding
		if (ai > thres)
			ai_thres -= thres;
		else if (ai < -thres)
			ai_thres += thres;
		else
			ai_thres = 0;

		// Hard-thresholding
		// if (fabs(ai) < thres)
		// 	ai_thres = 0;

		mask[i] += (ai_thres != 0);

		if (ai != 0.0)
			Xnew.set_row(i, Xnew.get_row(i) * (ai_thres / ai));
	}

	return Xnew;
 }



template <class TScalar, unsigned int Dimension>
TScalar
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::_QdiffTerm(
		const VNLMatrixType& Xtest, const VNLMatrixType& X,
		const VNLMatrixList& Atest, const VNLMatrixList& A,
		const VNLMatrixList& Ttest, const VNLMatrixList& T,
		TScalar stepXA, TScalar stepT)
{

	// for variables updated with a standard gradient step
	TScalar valT = 0.0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		for (unsigned int d = 0; d < Dimension; d++)
		{
			valT += dot_product(m_GradTemplate_L2[i].get_column(d),m_GradTemplate_Sob[i].get_column(d)); 
		}
	}
	valT *= ( -0.5 * stepT );

	TScalar valX = m_GradPos.fro_norm();
	valX *= ( -0.5 * stepXA * valX );

	// for variables updated with a soft-thresholded gradient step
	VNLMatrixList diffA(m_NumberOfSubjects);
	for (unsigned int s = 0; s < m_NumberOfSubjects; s++)
		diffA[s] = Atest[s] - A[s];

	// Compute first term of the form (xtest - x)^t * Grad
	TScalar termVal1 = 0.0;
	for (unsigned int s = 0; s < m_NumberOfSubjects; s++)
	{
		VNLMatrixType& diffA_s = diffA[s];
		const VNLMatrixType& gradA_s = m_GradMom[s];
		for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			termVal1 += dot_product(diffA_s.get_row(i), gradA_s.get_row(i));
	}

	// Compute second term of the form ||xtest - x||
	TScalar termVal2 = 0.0;
	for (unsigned int s = 0; s < m_NumberOfSubjects; s++)
	{
		VNLMatrixType& diffA_s = diffA[s];
		for (unsigned int i = 0; i < m_NumberOfCPs; i++)
			termVal2 += dot_product(diffA_s.get_row(i), diffA_s.get_row(i));
	}

	// sum all this
	TScalar val = termVal1 + termVal2 / (2.0 * stepXA) + valT + valX;

	return val;
}


template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ConvolveGradTemplate(VNLMatrixList& tempPoints)
 {

	if ( m_SmoothingKernelWidth < 1e-20 )
	{
		m_GradTemplate_Sob = m_GradTemplate_L2;
		return;
	}

	m_GradTemplate_Sob.resize(m_NumberOfObjects);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject();
	momKernelObj->SetKernelWidth(m_SmoothingKernelWidth);

	for (unsigned int i = 0; i < m_NumberOfObjects; i++)
	{
		if ( m_Template->GetObjectList()[i]->IsMesh() )
		{
			momKernelObj->SetSources(tempPoints[i]);
			momKernelObj->SetWeights(m_GradTemplate_L2[i]);

			m_GradTemplate_Sob[i] = momKernelObj->Convolve(tempPoints[i]);
		}
		else
			m_GradTemplate_Sob[i] = m_GradTemplate_L2[i];

	}

	delete momKernelObj;
	return;		
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::InitializeControlPointsAndMomentas()
 {
	if (m_NumberOfCPs == 0)
	{
		// create initial control points
		this->GenerateInitialControlPoints();

		m_NumberOfCPs = m_ControlPoints.rows();
		std::cout << "Generated " << m_NumberOfCPs << " control points" << std::endl;
		std::cout << "Momentas (re)set to zero" << std::endl;

		m_InitialMomentas.resize(m_NumberOfSubjects);
		for (int s = 0; s < m_NumberOfSubjects; s++)
		{
			m_InitialMomentas[s].set_size(m_NumberOfCPs, Dimension);
			m_InitialMomentas[s].fill(0);
		}
	}
	else 
	{	
		std::cout << "Using predefined set of " << m_NumberOfCPs << " control points" << std::endl;

		if ( (m_InitialMomentas.size() == m_NumberOfSubjects) && (m_InitialMomentas[0].rows() == m_NumberOfCPs) )
		{
			std::cout << "Using predefined set of momenta" << std::endl;
		}
		else
		{
			if (m_InitialMomentas.size() > 0)
				std::cout << "Warning: initial momenta file has incompatible number of subjects and/or vectors. Initial momenta reset to zero" << std::endl;

			m_InitialMomentas.resize(m_NumberOfSubjects);
			for (int s = 0; s < m_NumberOfSubjects; s++)
			{
				m_InitialMomentas[s].set_size(m_NumberOfCPs, Dimension);
				m_InitialMomentas[s].fill(0);
			}
		}
	}
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
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
typename SparseDiffeoAtlasEstimator<TScalar, Dimension>::FunctionalValuesType*
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ComputeFunctional(VNLMatrixType& X0, VNLMatrixList& Mom0, VNLMatrixList& TempL)
 {
	m_Template->SetImageAndPointData(TempL);

	FunctionalValuesType* funcValues = new FunctionalValuesType(m_NumberOfObjects, m_NumberOfSubjects, m_Weights);

	for (int s = 0; s < m_NumberOfSubjects; s++)
	{
		std::vector<TScalar> funcValuesSubject;
		bool oob = this->ComputeFunctionalSubject(X0, Mom0[s], m_TargetList[s], funcValuesSubject);
		// if (m_Def->OutOfBox() || m_Template->OutOfBox()): with cloning m_Def has never been updated: m_Def->OutOfBox() has never been set up
		if (oob)
		{
			funcValues->SetOutOfBox();
			return funcValues;
		}

		for (int i = 0; i < m_NumberOfObjects; i++)
			funcValues->SetDataTerm(funcValuesSubject[i],s,i);

		funcValues->SetRegularity(funcValuesSubject[m_NumberOfObjects],s);
		funcValues->SetSparsity(funcValuesSubject[m_NumberOfObjects+1],s);
	}

	DeformableObjectList templateList = m_Template->GetObjectList();

	funcValues->Update();

	return funcValues;
 }



template <class TScalar, unsigned int Dimension>
bool
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ComputeFunctionalSubject(VNLMatrixType& X0, VNLMatrixType& Mom0,
		DeformableMultiObjectType* targetObj, std::vector<TScalar>& funcValuesSubject)
 {

	DiffeosType* subjectDef = m_Def->Clone();
	subjectDef->SetStartPositions(X0);
	subjectDef->SetStartMomentas(Mom0);
	subjectDef->Update();

	if (subjectDef->OutOfBox())
	{
		delete subjectDef;
		return true;
	}

	DeformableMultiObjectType* subjectTemplate = m_Template->Clone();
	subjectTemplate->SetDiffeos(subjectDef);
	subjectTemplate->Update();

	if (subjectTemplate->OutOfBox())
	{
		delete subjectTemplate;
		delete subjectDef;
		return true;
	}

	funcValuesSubject.resize(m_NumberOfObjects+2);

	std::vector<TScalar> DataTerm = subjectTemplate->ComputeMatch(targetObj);
	DeformableObjectList templateList = subjectTemplate->GetObjectList();
	for (int i = 0; i < m_NumberOfObjects; i++)
		funcValuesSubject[i] = DataTerm[i] / ( 2.0 * templateList[i]->GetDataSigmaSquared() );


	// compute the regularity term: the squared norm of the initial velocity speed (halfed)
	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();

	KernelType* momKernelObj  = kfac->CreateKernelObject();
	momKernelObj->SetKernelWidth(subjectDef->GetKernelWidth());
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

	funcValuesSubject[m_NumberOfObjects] = Reg;
	funcValuesSubject[m_NumberOfObjects+1] = SparsityPrior;

	delete subjectTemplate;
	delete subjectDef;

	// return funcValues;
	return false;
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ComputeGradient(VNLMatrixType& X0, VNLMatrixList& Mom0, VNLMatrixList& TempL)
 {
	m_Template->SetImageAndPointData(TempL);

	m_GradPos.fill(0.0);
	for (int i = 0; i < m_NumberOfObjects; i++)
		m_GradTemplate_L2[i].fill(0.0);

	for (int s = 0; s < m_NumberOfSubjects; s++)
	{
		VNLMatrixType dPos;
		VNLMatrixType dMom;
		VNLMatrixList dTempL(m_NumberOfObjects);
		this->ComputeGradientSubject(X0, Mom0[s], dPos, dMom, dTempL, m_TargetList[s]);

		m_GradPos += m_Weights[s]*dPos;
		dMom *= m_Weights[s];
		m_GradMom[s] = dMom;

		for (int i = 0; i < m_NumberOfObjects; i++)
			m_GradTemplate_L2[i] += (m_Weights[s]*dTempL[i]);
			// m_GradTemplate_L2[i] += (dTempL[i] / m_NumberOfSubjects);
	}

	// m_GradPos /= m_NumberOfSubjects;
	this->ConvolveGradTemplate(TempL);

 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ComputeGradientSubject(VNLMatrixType& X0, VNLMatrixType& Mom0,
		VNLMatrixType& dPos, VNLMatrixType& dMom, VNLMatrixList& dTempL,
		DeformableMultiObjectType* targetObj)
 {
	// Compute forward trajectory
	DiffeosType* subjectDef = m_Def->Clone();
	subjectDef->SetStartPositions(X0);
	subjectDef->SetStartMomentas(Mom0);
	subjectDef->Update();

	if (subjectDef->OutOfBox())
		std::cerr << "out of box in Compute Gradient (this should not be...)" << std::endl;

	VNLMatrixList posT = subjectDef->GetTrajectoryPositions();
	VNLMatrixList momT = subjectDef->GetTrajectoryMomentas();

	// Compute eta(t)
	DeformableMultiObjectType* subjectTemplate = m_Template->Clone();
	subjectTemplate->SetDiffeos(subjectDef);
	subjectTemplate->Update();

	VNLMatrixList InitialConditions = subjectTemplate->ComputeMatchGradient(targetObj);

	DeformableObjectList templateList = subjectTemplate->GetObjectList();	
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		if (templateList[i]->IsOfLandmarkKind())
			InitialConditions[i] /= 2.0 * templateList[i]->GetDataSigmaSquared();

		else
			throw std::runtime_error("Unknown object type");
	}

	VNLMatrixList VectorsT;
	VNLMatrixList PointsT;

	subjectTemplate->TransportAlongGeodesic(InitialConditions,VectorsT,PointsT);

	typedef KernelFactory<TScalar, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kfac = KernelFactoryType::Instantiate();
	KernelType* momKernelObj  = kfac->CreateKernelObject();

	// std::cout << "compute gradient w.r.t template data" << std::endl;
	int counter = 0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		int NumPoints_i = templateList[i]->GetNumberOfPoints();

		if (templateList[i]->IsOfLandmarkKind())
					dTempL[i] = VectorsT[0].get_n_rows(counter, NumPoints_i);
		else
			throw std::runtime_error("Unknown object type");

		counter += NumPoints_i;
	}

	typedef BackwardGradientPropagatorUniversal<TScalar, Dimension> BackwardEstimator;
	BackwardEstimator* gradbwd = new BackwardEstimator();
	gradbwd->SetPointsTrajectory(PointsT);
	gradbwd->SetVectorsTrajectory(VectorsT);
	gradbwd->SetDiffeos(subjectDef);
	gradbwd->Update();

	dPos = gradbwd->GetGradientPosAt(0);
	dMom = gradbwd->GetGradientMomAt(0);

	delete gradbwd;

	// gradient of the regularity term
	momKernelObj->SetKernelWidth(subjectDef->GetKernelWidth());
	momKernelObj->SetSources(X0);
	momKernelObj->SetWeights(Mom0);

	VNLMatrixType kMom = momKernelObj->Convolve(X0);

	dMom += kMom;

	dPos = dPos + momKernelObj->ConvolveGradient(X0, Mom0);

	delete momKernelObj;

	delete subjectTemplate;
	delete subjectDef;

 }



//
// Multi-threaded components
//

template <class TScalar, unsigned int Dimension>
typename SparseDiffeoAtlasEstimator<TScalar, Dimension>::FunctionalValuesType*
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ThreadedComputeFunctional(VNLMatrixType& X0, VNLMatrixList& Mom0, VNLMatrixList& TempL)
 {
	if (m_NumberOfThreads < 2)
	{
		return this->ComputeFunctional(X0, Mom0, TempL);
	}

	m_MT_SubjectCounter = 0;

	m_MT_X0 = X0;
	m_MT_Mom0 = Mom0;
	m_Template->SetImageAndPointData(TempL);

	delete m_MT_FuncValues;
	m_MT_FuncValues = new FunctionalValuesType(m_NumberOfObjects, m_NumberOfSubjects, m_Weights);

	itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

	//unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
	unsigned int numThreads = m_NumberOfThreads;
	if (numThreads > m_NumberOfSubjects)
		numThreads = m_NumberOfSubjects;
	if (numThreads < 2)
		numThreads = 2;

	threader->SetNumberOfThreads(numThreads);
	threader->SetSingleMethod(
			&SparseDiffeoAtlasEstimator::_functionalThread, (void*)this);
	threader->SingleMethodExecute();

	if (m_MT_FuncValues->IsOutOfBox())
	{
		return new FunctionalValuesType(*m_MT_FuncValues);
	}

	DeformableObjectList templateList = m_Template->GetObjectList();

	FunctionalValuesType* fVal = new FunctionalValuesType(*m_MT_FuncValues);
	fVal->Update();

	return fVal;
 }



template <class TScalar, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::_functionalThread(void* arg)
 {
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast<ThreadInfoType*>( arg );
	SparseDiffeoAtlasEstimator* obj = static_cast<SparseDiffeoAtlasEstimator*>(
			infoStruct->UserData);

	unsigned int s = obj->m_NumberOfSubjects;

	while (true)
	{
		obj->m_Mutex.Lock();
		s = obj->m_MT_SubjectCounter++;
		obj->m_Mutex.Unlock();

		if (s >= obj->m_NumberOfSubjects)
			break;

		std::vector<TScalar> values;
		bool oob = obj->ComputeFunctionalSubject(
				obj->m_MT_X0, obj->m_MT_Mom0[s], obj->m_TargetList[s], values);

		if (oob)
		{
			obj->m_Mutex.Lock();
			obj->m_MT_FuncValues->SetOutOfBox();
			obj->m_Mutex.Unlock();
			break;
		}

		for (int i = 0; i < obj->m_NumberOfObjects; i++)
		{
			obj->m_Mutex.Lock();
			obj->m_MT_FuncValues->SetDataTerm(values[i],s,i);
			obj->m_Mutex.Unlock();
		}
		obj->m_Mutex.Lock();
		obj->m_MT_FuncValues->SetRegularity(values[obj->m_NumberOfObjects],s);
		obj->m_MT_FuncValues->SetSparsity(values[obj->m_NumberOfObjects+1],s);
		obj->m_Mutex.Unlock();

	}

	return ITK_THREAD_RETURN_VALUE;
 }



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ThreadedComputeGradient(VNLMatrixType& X0, VNLMatrixList& Mom0, VNLMatrixList& TempL)
 {
	if (m_NumberOfThreads < 2)
	{
		this->ComputeGradient(X0, Mom0, TempL);
		return;
	}

	m_MT_SubjectCounter = 0;

	m_MT_X0 = X0;
	m_MT_Mom0 = Mom0;
	m_Template->SetImageAndPointData(TempL);

	m_GradPos.fill(0.0);
	m_GradMom.resize(m_NumberOfSubjects);
	m_GradTemplate_L2.resize(m_NumberOfObjects);
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		m_GradTemplate_L2[i].set_size(m_TemplateData[i].rows(), m_TemplateData[i].columns());
		m_GradTemplate_L2[i].fill(0.0);
	}

	itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

	//unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
	unsigned int numThreads = m_NumberOfThreads;
	if (numThreads > m_NumberOfSubjects)
		numThreads = m_NumberOfSubjects;
	if (numThreads < 2)
		numThreads = 2;

	threader->SetNumberOfThreads(numThreads);
	threader->SetSingleMethod(
			&SparseDiffeoAtlasEstimator::_gradientThread, (void*)this);
	threader->SingleMethodExecute();

	// m_GradPos /= m_NumberOfSubjects;
	this->ConvolveGradTemplate(TempL);	

 }



template <class TScalar, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::_gradientThread(void* arg)
 {
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast<ThreadInfoType*>( arg );
	SparseDiffeoAtlasEstimator* obj = static_cast<SparseDiffeoAtlasEstimator*>(
			infoStruct->UserData);

	unsigned int s = obj->m_NumberOfSubjects;

	while (true)
	{
		obj->m_Mutex.Lock();
		s = obj->m_MT_SubjectCounter++;
		obj->m_Mutex.Unlock();

		if (s >= obj->m_NumberOfSubjects)
			break;

		VNLMatrixType dPos;
		VNLMatrixType dMom;
		VNLMatrixList dTempL(obj->m_NumberOfObjects);
		obj->ComputeGradientSubject(
				obj->m_MT_X0, obj->m_MT_Mom0[s], dPos, dMom, dTempL,
				obj->m_TargetList[s]);

		obj->m_Mutex.Lock();
		obj->m_GradPos += obj->m_Weights[s]*dPos;
		obj->m_Mutex.Unlock();

		// dMom /= obj->m_NumberOfSubjects;
		dMom *= obj->m_Weights[s];

		obj->m_Mutex.Lock();
		obj->m_GradMom[s] = dMom;
		obj->m_Mutex.Unlock();

		for (int i = 0; i < obj->m_NumberOfObjects; i++)
		{
			obj->m_Mutex.Lock();
			obj->m_GradTemplate_L2[i] += (obj->m_Weights[s]*dTempL[i]);
			// obj->m_GradTemplate_L2[i] += (dTempL[i] / obj->m_NumberOfSubjects);
			obj->m_Mutex.Unlock();
		}

	}

	return ITK_THREAD_RETURN_VALUE;
 }



template <class TScalar, unsigned int Dimension>
TScalar
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::_ThreadedQdiffTerm(
		const VNLMatrixType& Xtest, const VNLMatrixType& X,
		const VNLMatrixList& Atest, const VNLMatrixList& A,
		const VNLMatrixList& Ttest, const VNLMatrixList& T,
		TScalar stepXA, TScalar stepT)
{
	if (m_NumberOfThreads < 2)
	{
		return this->_QdiffTerm(Xtest, X, Atest, A, Ttest, T, stepXA, stepT);
	}

	//
	// Handle the global parameters (X and T)
	//

	TScalar valT = 0.0;
	for (int i = 0; i < m_NumberOfObjects; i++)
	{
		for (unsigned int d = 0; d < Dimension; d++)
		{
			valT += dot_product(m_GradTemplate_L2[i].get_column(d),m_GradTemplate_Sob[i].get_column(d)); 
		}
	}
	valT *= ( -0.5 * stepT );

	TScalar valX = m_GradPos.fro_norm();
	valX *= ( -0.5 * stepXA * valX );

	TScalar val = valT + valX;

	//
	// Process the per-subject A[s] in separate threads
	//

	m_MT_SubjectCounter = 0;

	m_MT_Qdiff_Atest = Atest;
	m_MT_Qdiff_A = A;
	m_MT_StepXA = stepXA;

	m_MT_Qdiff_Value = 0;

	itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

	//unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
	unsigned int numThreads = m_NumberOfThreads;
	if (numThreads > m_NumberOfSubjects)
		numThreads = m_NumberOfSubjects;
	if (numThreads < 2)
		numThreads = 2;

	threader->SetNumberOfThreads(numThreads);
	threader->SetSingleMethod(
			&SparseDiffeoAtlasEstimator::_QdiffThread, (void*)this);
	threader->SingleMethodExecute();

	val += m_MT_Qdiff_Value;

	return val;
}



template <class TScalar, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::_QdiffThread(void* arg)
{
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast<ThreadInfoType*>( arg );
	SparseDiffeoAtlasEstimator* obj = static_cast<SparseDiffeoAtlasEstimator*>(
			infoStruct->UserData);

	unsigned int s = obj->m_NumberOfSubjects;

	VNLMatrixList& Atest = obj->m_MT_Qdiff_Atest;
	VNLMatrixList& A = obj->m_MT_Qdiff_A;
	VNLMatrixList& gradA = obj->m_GradMom;

	TScalar stepXA = obj->m_MT_StepXA;

	while (true)
	{
		obj->m_Mutex.Lock();
		s = obj->m_MT_SubjectCounter++;
		obj->m_Mutex.Unlock();

		if (s >= obj->m_NumberOfSubjects)
			break;

		VNLMatrixType diffA_s =  Atest[s] - A[s];

		const VNLMatrixType& gradA_s = gradA[s];

		TScalar val1 = 0;
		for (unsigned int i = 0; i < obj->m_NumberOfCPs; i++)
			val1 += dot_product(diffA_s.get_row(i), gradA_s.get_row(i));

		TScalar val2 = 0;
		for (unsigned int i = 0; i < obj->m_NumberOfCPs; i++)
			val2 += dot_product(diffA_s.get_row(i), diffA_s.get_row(i));
		val2 /= (2.0 * stepXA);

		obj->m_Mutex.Lock();
		obj->m_MT_Qdiff_Value += val1 + val2;
		obj->m_Mutex.Unlock();
	}

	return ITK_THREAD_RETURN_VALUE;
}



template <class TScalar, unsigned int Dimension>
void
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::ThreadedGradientDescentStep(
		VNLMatrixType& Xtest, const VNLMatrixType& X,
		VNLMatrixList& Atest, const VNLMatrixList& A,
		VNLMatrixList& Ttest, const VNLMatrixList& T,
		std::vector<int>& activeCPMask,
		TScalar stepXA, TScalar stepT)
{
	if (m_NumberOfThreads < 2)
	{
		this->GradientDescentStep(Xtest, X, Atest, A, Ttest, T, activeCPMask, stepXA, stepT);
		return;
	}

	Xtest = X - m_GradPos*stepXA;

	for (int i = 0; i < m_NumberOfObjects; i++)
		Ttest[i] = T[i] - m_GradTemplate_Sob[i]*stepT;

	m_MT_SubjectCounter = 0;

	m_MT_Mom0 = A;
	m_MT_StepXA = stepXA;

	m_MT_ActiveCPMask.assign(m_NumberOfCPs, 0);

	itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

	//unsigned int numThreads = threader->GetGlobalMaximumNumberOfThreads() / 4 * 3;
	unsigned int numThreads = m_NumberOfThreads;
	if (numThreads > m_NumberOfSubjects)
		numThreads = m_NumberOfSubjects;
	if (numThreads < 2)
		numThreads = 2;

	threader->SetNumberOfThreads(numThreads);
	threader->SetSingleMethod(
			&SparseDiffeoAtlasEstimator::_descentStepThread, (void*)this);
	threader->SingleMethodExecute();

	Atest = m_MT_Mom0;

	if (m_SparsityPrior == 0)
		activeCPMask.assign(m_NumberOfCPs, 1);
	else
		activeCPMask = m_MT_ActiveCPMask;
}



template <class TScalar, unsigned int Dimension>
ITK_THREAD_RETURN_TYPE
SparseDiffeoAtlasEstimator<TScalar, Dimension>
::_descentStepThread(void* arg)
 {
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast<ThreadInfoType*>( arg );
	SparseDiffeoAtlasEstimator* obj = static_cast<SparseDiffeoAtlasEstimator*>(
			infoStruct->UserData);

	unsigned int s = obj->m_NumberOfSubjects;

	VNLMatrixList& A = obj->m_MT_Mom0;
	VNLMatrixList& gradA = obj->m_GradMom;

	TScalar stepXA = obj->m_MT_StepXA;

	while (true)
	{
		obj->m_Mutex.Lock();
		s = obj->m_MT_SubjectCounter++;
		obj->m_Mutex.Unlock();

		if (s >= obj->m_NumberOfSubjects)
			break;

		std::vector<int> cpMask_s;

		VNLMatrixType Anew_s;
		if (obj->m_SparsityPrior == 0)
		{
			Anew_s = A[s] - gradA[s] * stepXA;
		}
		else
		{
			cpMask_s.assign(obj->m_NumberOfCPs, 0);
			Anew_s = obj->SoftThresholdUpdate(A[s], gradA[s], stepXA, cpMask_s);
		}

		obj->m_Mutex.Lock();
		A[s] = Anew_s;
		obj->m_Mutex.Unlock();

		if (obj->m_SparsityPrior != 0)
		{
			for (unsigned int i = 0; i < obj->m_NumberOfCPs; i++)
			{
				obj->m_Mutex.Lock();
				obj->m_MT_ActiveCPMask[i] += cpMask_s[i];
				obj->m_Mutex.Unlock();
			}
		}
	}

	return ITK_THREAD_RETURN_VALUE;

 }



#endif /* _SparseDiffeoAtlasEstimator_txx */
