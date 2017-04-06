/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AtlasFunctionalValues_h
#define _AtlasFunctionalValues_h

#include <vector>
#include <vnl/vnl_vector.h>

#include <string>
#include <sstream>
#include <iostream>

/**
 *	\brief 		Values of the functional during atlas estimation.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The AtlasFunctionalValues class stores the values of the different
 *				parts of the functional at a given iteration. The functional is :
 *				\f[
 *				\functional =
 *				\sum_{i=1}^{N} \left\lbrace \sum_{k=1}^{\nbobj}\frac{1}{2\sigma_k^2} D(X_{0,k}^i(1),X_i)^2 +
 *				\frac{1}{2}\sum_{p,q} (\alpha_{0,p}^i)^T K(c_{0,p},c_{0,q})\alpha_{0,q}^i +
 *				\sparsityprior \sum_{p} \norm{\alpha_{0,p}^i} \right\rbrace
 *				\f]
 *				where :
 *				\li \f$\frac{1}{2\sigma_k^2} D(X_{0,k}^i(1),X_i)^2\f$ denotes the fidelity-to-data term of the k-th object of the i-th subject (AtlasFunctionalValues::m_DataTerm) ;
 *				\li \f$\frac{1}{2}\sum_{p,q} (\alpha_{0,p}^i)^T K(c_{0,p},c_{0,q})\alpha_{0,q}^i \f$ denotes the regularity term of the i-th deformation, which deforms the template objects to the objects of the i-th subject (AtlasFunctionalValues::m_Regularity) ;
 *				\li \f$\sparsityprior \sum_{p} \norm{\alpha_{0,p}^i}\f$ denotes the \f$\Lone\f$ prior (AtlasFunctionalValues::m_Sparsity) . \n
 *				This class goes with the SparseDiffeoAtlasEstimator class.
 */
template <class TScalar>
class AtlasFunctionalValues
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type (std).
	typedef std::vector<TScalar> VectorType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Constructor with initialization of the number of deformable objects
	/// (resp. the number of subjects) to \e nobj (resp. \e to nsubj).
	AtlasFunctionalValues(int nobj, int nsubj, std::string weights)
	{
		m_NumberOfObjects = nobj;
		m_NumberOfSubjects = nsubj;
		SetWeights(weights);

		m_DataTerm.resize(m_NumberOfSubjects);
		for (int s = 0; s < m_NumberOfSubjects; s++)
			m_DataTerm[s].resize(m_NumberOfObjects,0.0);

		m_Regularity.resize(m_NumberOfSubjects, 0.0);
		m_Sparsity.resize(m_NumberOfSubjects, 0.0);

		m_DataTermPerObject.resize(m_NumberOfObjects, vnl_huge_val(1.0));
		m_TotalDataTerm = vnl_huge_val(1.0);
		m_TotalRegularity = vnl_huge_val(1.0);
		m_TotalSparsity = vnl_huge_val(1.0);
		m_TotalValue = vnl_huge_val(1.0);
		m_TotalL2Value = vnl_huge_val(1.0);

		m_OutOfBox = false;
	}

	/// Copy constructor.
	AtlasFunctionalValues(const AtlasFunctionalValues& other)
	{
		m_NumberOfObjects = other.m_NumberOfObjects;
		m_NumberOfSubjects = other.m_NumberOfSubjects;

		m_DataTerm = other.m_DataTerm;
		m_Regularity = other.m_Regularity;
		m_Weights = other.m_Weights;
		m_Sparsity = other.m_Sparsity;

		m_OutOfBox = other.m_OutOfBox;

		m_DataTermPerObject = other.m_DataTermPerObject;
		m_TotalDataTerm = other.m_TotalDataTerm;
		m_TotalRegularity = other.m_TotalRegularity;
		m_TotalSparsity = other.m_TotalSparsity;
		m_TotalValue = other.m_TotalValue;
		m_TotalL2Value = other.m_TotalL2Value;
	}

	/// Make a copy of the object.
	AtlasFunctionalValues* Clone() { return new AtlasFunctionalValues(*this); }

	~AtlasFunctionalValues() {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the whole fidelity-to-data terms.
	inline std::vector< VectorType > GetDataTerm() { return m_DataTerm; }
	/// Sets the fidelity-to-data terms of the \e s-th subject to \e dt.
	inline void SetDataTerm(VectorType dt, int s) { m_DataTerm[s] = dt; }
	/// Sets the fidelity-to-data term of the \e s-th subject to \e r.
	inline void SetDataTerm(TScalar r, int s, int i) { m_DataTerm[s][i] = r; }

	/// Sets the regularity term of the \e s-th subject to \e r.
	inline void SetRegularity(TScalar r, int s) { m_Regularity[s] = r; }

	/// Sets the \f$\Lone\f$ prior of the \e s-th subject to \e r.
	inline void SetSparsity(TScalar r, int s) { m_Sparsity[s] = r; }

	/// Returns true if any point of the trajectory is out of box, false otherwise.
	inline bool IsOutOfBox() { return m_OutOfBox; }
	/// Sets the boolean OutOfBox to true.
	inline void SetOutOfBox() { m_OutOfBox = true; }

	/// Returns a vector of \f$\nbobj\f$fidelity-to-data terms.
	inline VectorType GetDataTermPerObject() { return m_DataTermPerObject; }
	/// Returns the total fidelity-to-data term.
	inline TScalar GetTotalDataTerm() { return m_TotalDataTerm; }

	/// Returns the total regularity term.
	inline TScalar GetTotalRegularity() { return m_TotalRegularity; }

	/// Returns the total regularity term.
	inline TScalar GetTotalSparsity() { return m_TotalSparsity; }

	/// Returns the value of the functional with the \f$\Lone\f$ prior.
	inline TScalar GetTotalValue() { return m_TotalValue; }

	/// Returns the value of the functional without the \f$\Lone\f$ prior.
	inline TScalar GetTotalL2Value() { return m_TotalL2Value; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the regularity weights for the \e s-th subject according to the paramsDiffeo file
	void SetWeights(std::string weights)
	{
		// Equal weights to each subject by default
		m_Weights.resize(m_NumberOfSubjects, 1.0/m_NumberOfSubjects);
		if(weights!="uniform")
		{
			std::stringstream ss(weights);
			TScalar w;
			int count = 0;
			while(ss >> w)
			{
				m_Weights[count] = w;
				if(ss.peek() == ',' | ss.peek() == ' ')
					ss.ignore();
				++count;
			}

			if(count != m_NumberOfSubjects)
			{
				// If the number of provided weights does not match the number of subjects, use uniform weights
				std::cout << "Number of weights not equal to the number of patients!" << std::endl;
				std::cout << "Using uniform Regularity weights" << std::endl;
				for(int j=0; j<m_Weights.size(); j++)
					m_Weights[j] = 1.0/m_NumberOfSubjects;
        	}
		}
	}

	/// Updates the different total values and the fidelity-to-data term per object.
	void Update()
	{
		if (m_OutOfBox)
			return;


		m_DataTermPerObject.resize(m_NumberOfObjects,0.0);
		m_TotalDataTerm = 0.0;
		for (int i = 0; i < m_NumberOfObjects; i++)
		{
			TScalar val = 0.0;
			for (int s = 0; s < m_NumberOfSubjects; s++)
				val += m_DataTerm[s][i];

			val = m_Weights[i]*val;					// Weighting Data Term according to provided weights
			m_DataTermPerObject[i] = val;
			m_TotalDataTerm += val;
		}

		m_TotalRegularity = 0.0;
		m_TotalSparsity = 0.0;
		for (int s = 0; s < m_NumberOfSubjects; s++)
		{
			m_TotalRegularity += m_Weights[s]*m_Regularity[s];		// Weighting Regularity Term according to provided weights
			m_TotalSparsity += m_Sparsity[s];
		}
		// m_TotalRegularity /= m_NumberOfSubjects;
		m_TotalSparsity /= m_NumberOfSubjects;

		m_TotalValue = m_TotalDataTerm + m_TotalRegularity + m_TotalSparsity;
		m_TotalL2Value = m_TotalDataTerm + m_TotalRegularity;
	}

	/// Displays the different parts of the functional with an indication of iteration \e iter.
	void PrintIter(const int iter)
	{
		std::cout << "Iter " << iter << "  >> Objective = " << m_TotalValue <<
				"\t(Data Term = " << m_TotalDataTerm <<
				"\tRegularity = " << m_TotalRegularity <<
				"\tSparsityPrior = " << m_TotalSparsity <<
				")" << std::endl;
	}



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Number of deformable objects.
	int m_NumberOfObjects;

	/// Number of subjects.
	int m_NumberOfSubjects;

	/// See AbstractDeformations::m_OutOfBox for details.
	bool m_OutOfBox;

	/// Vector of \f$\nbsubj\f$ vectors of \f$\nbobj\f$ fidelity-to-data terms.
	std::vector< VectorType > m_DataTerm;

	/// Vector of \f$\nbsubj\f$ regularity terms.
	VectorType m_Regularity;

	/// Vector of \f$\nbsubj\f$ regularity weights.
	VectorType m_Weights;

	/// Vector of \f$\nbsubj\f$ \f$\Lone\f$ prior terms.
	VectorType m_Sparsity;

	/// Vector of \f$\nbobj\f$ fidelity-to-data terms.
	VectorType m_DataTermPerObject;

	/// Double sum on the subjects and the deformable objects of the fidelity-to-data terms.
	TScalar m_TotalDataTerm;

	/// Sum on the subjects of the regularity terms.
	TScalar m_TotalRegularity;

	/// Sum on the subjects of the \f$\Lone\f$ priors.
	TScalar m_TotalSparsity;

	/// Total value of the functional (\f$\Lone\f$ prior included).
	TScalar m_TotalValue;

	/// Total value of the functional (\f$\Lone\f$ prior not included).
	TScalar m_TotalL2Value;



}; /* class AtlasFunctionalValues */


#endif /* _AtlasFunctionalValues_h */
