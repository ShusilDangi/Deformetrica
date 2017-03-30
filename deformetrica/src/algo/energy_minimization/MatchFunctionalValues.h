/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _MatchFunctionalValues_h
#define _MatchFunctionalValues_h

#include <vector>
#include <vnl/vnl_vector.h>


/**
 *	\brief 		Values of the functional during registration
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The MatchFunctionalValues class stores the values of the different
 *				parts of the functional at a given iteration. We recall that the functional is :
 *				\f[
 *				\functional =  \sum_{k=1}^{\nbobj}\left\lbrace
 *				\frac{1}{2\sigma_k^2} D(S_{0,k}(1),T_k)^2 + \sum_{p,q} (\alpha_{0,p})^T K(c_{0,p},c_{0,q})\alpha_{0,q} + \sparsityprior \sum_{p} \norm{\alpha_{0,p}^i}
 *				\right\rbrace
 *				\f]
 *				where :
 *				\li \f$\frac{1}{2\sigma_k^2} D(S_{0,k}(1),T_k)^2\f$ denotes the fidelity-to-data term between the k-th deformed source object and the k-th target object (MatchFunctionalValues::m_DataTerm) ;
 *				\li \f$\frac{1}{2}\sum_{p,q} (\alpha_{0,p})^T K(c_{0,p},c_{0,q})\alpha_{0,q} \f$ denotes the regularity term of the deformation, which deforms the source objects to the target objects (MatchFunctionalValues::m_Regularity) ;
 *				\li \f$\sparsityprior \sum_{p} \norm{\alpha_{0,p}}\f$ denotes the \f$\Lone\f$ prior (MatchFunctionalValues::m_Sparsity) . \n
 *				This class goes with the SparseDiffeoMatcher class.
 */
template <class TScalar>
class MatchFunctionalValues
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

	/// Constructor with initialization of the number of deformable objects to \e nobj.
	MatchFunctionalValues(int nobj)
	{
		m_NumberOfObjects = nobj;

		m_DataTerm.resize(m_NumberOfObjects, 0);			
		m_Regularity = vnl_huge_val(1.0);
		m_Sparsity = vnl_huge_val(1.0);

		m_TotalDataTerm = vnl_huge_val(1.0);
		m_TotalValue = vnl_huge_val(1.0);
		m_TotalL2Value = vnl_huge_val(1.0);

		m_OutOfBox = false;
	}

	/// Copy constructor.
	MatchFunctionalValues(const MatchFunctionalValues& other)
	{
		m_NumberOfObjects = other.m_NumberOfObjects;

		m_DataTerm = other.m_DataTerm;
		m_Regularity = other.m_Regularity;
		m_Sparsity = other.m_Sparsity;

		m_OutOfBox = other.m_OutOfBox;

		m_TotalDataTerm = other.m_TotalDataTerm;
		m_TotalValue = other.m_TotalValue;
		m_TotalL2Value = other.m_TotalL2Value;
	}

	/// Make a copy of the object.
	MatchFunctionalValues* Clone() { return new MatchFunctionalValues(*this); }

	~MatchFunctionalValues() {}



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the vector of the fidelity-to-data term.
	inline VectorType GetDataTerm() { return m_DataTerm; }
	/// Sets the \e i-th fidelity-to-data term to \e r.
	inline void SetDataTerm(TScalar r, int i) { m_DataTerm[i] = r; }
	/// Sets the whole fidelity-to-data terms to \e t.
	inline void SetDataTerm(VectorType t) { m_DataTerm = t; }

	/// Returns the regularity term.
	inline TScalar GetRegularity() { return m_Regularity; }
	/// Sets the regularity term to \e r.
	inline void SetRegularity(TScalar r) { m_Regularity = r; }

	/// Returns the \f$\Lone\f$ prior.
	inline TScalar GetSparsity() { return m_Sparsity; }
	/// Set the \f$\Lone\f$ prior to \e r.
	inline void SetSparsity(TScalar r) { m_Sparsity = r; }

	/// Returns true if any point of the trajectory is out of box, false otherwise.
	inline bool IsOutOfBox() { return m_OutOfBox; }
	/// Sets the boolean OutOfBox to true.
	inline void SetOutOfBox() { m_OutOfBox = true; }

	/// Returns the total fidelity-to-data term.
	inline TScalar GetTotalDataTerm() { return m_TotalDataTerm; }

	/// Returns the value of the functional with the \f$\Lone\f$ prior.
	inline TScalar GetTotalValue() { return m_TotalValue; }

	/// Returns the value of the functional without the \f$\Lone\f$ prior.
	inline TScalar GetTotalL2Value() { return m_TotalL2Value; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates the different total values.
	void Update()
	{
		if (m_OutOfBox)
			return;

		m_TotalDataTerm = 0.0;
		for (int i = 0; i < m_NumberOfObjects; i++)
			m_TotalDataTerm += m_DataTerm[i];

		m_TotalValue = m_TotalDataTerm + m_Regularity + m_Sparsity;
		m_TotalL2Value = m_TotalDataTerm + m_Regularity;
	}


	/// Displays the different parts of the functional with an indication of iteration \e iter.
	void PrintIter(const int iter)
	{
		std::cout << "Iter " << iter << "  >> Objective = " << m_TotalValue <<
				"\t(Data Term = " << m_TotalDataTerm <<
				"\tRegularity = " << m_Regularity <<
				"\tSparsityPrior = " << m_Sparsity << ")" << std::endl;
	}



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Number of deformable objects.
	int m_NumberOfObjects;

	/// See AbstractDeformations::m_OutOfBox for details.
	bool m_OutOfBox;

	/// Vector containing the fidelity-to-data term of each deformable object.
	VectorType m_DataTerm;

	/// Value of the regularity term appearing in the functional.
	TScalar m_Regularity;

	/// Value of the \f$\Lone\f$ penalty.
	TScalar m_Sparsity;

	/// Sum of the fidelity-to-data term of each deformable object.
	TScalar m_TotalDataTerm;

	/// Total value of the functional (\f$\Lone\f$ prior included).
	TScalar m_TotalValue;

	/// Total value of the functional (\f$\Lone\f$ prior not included).
	TScalar m_TotalL2Value;



}; /* class MatchFunctionalValues */


#endif /* _MatchFunctionalValues_h */
