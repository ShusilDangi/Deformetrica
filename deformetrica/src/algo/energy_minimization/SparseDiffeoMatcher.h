/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoMatcher_h
#define _SparseDiffeoMatcher_h

#include "itkImage.h"
#include "itkVector.h"
#include "Diffeos.h"
#include "DeformableMultiObject.h"

#include "readMatrixDLM.txx"

#include "MatchFunctionalValues.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include <vector>

/**
 *	\brief 		Registration class
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The SparseDiffeoMatcher class enables to register two configurations of objects embedded in the ambient space.\n \n
 */
template <class TScalar, unsigned int Dimension>
class SparseDiffeoMatcher
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type
	typedef vnl_vector<TScalar> VNLVectorType;
	/// Matrix type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// List of matrices type.
	typedef std::vector<VNLMatrixType> VNLMatrixList;

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	/// Vector of deformable objects type.
	typedef std::vector<DeformableObjectType*> DeformableObjectList;
	/// Deformable multi-object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;

	/// Functional values type.
	typedef MatchFunctionalValues<TScalar> FunctionalValuesType;

	/// Deformation type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;

	/// Type of optimization method.
	typedef enum
	{
		null,			/*!< Null value. */
		GradientDescent,/*!< Usual line search gradient descent (see GradientDescentAndISTA()). */
		ISTA,			/*!< Line search gradient descent compatible with \f$L^1\f$ penalty terms (see GradientDescentAndISTA()). */
		F_ISTA,			/*!< Faster gradient scheme built on ISTA and Nesterov's scheme (see FISTA()). */
	} OptimizationMethodType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	SparseDiffeoMatcher();

	~SparseDiffeoMatcher();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the maximum of descent iterations to \e n.
	inline void SetMaxIterations(unsigned int n);

	/// Sets the maximum number of iterations to \e n.
	inline void SetMaxLineSearchIterations(unsigned int n) { m_MaxLineSearchIterations = n; }

	/// Sets the expand parameter to \e d.
	inline void SetAdaptiveExpand(TScalar d) { m_AdaptiveExpand = d; }

	/// Sets the shrink parameter to \e d.
	inline void SetAdaptiveShrink(TScalar d) { m_AdaptiveShrink = d; }

	/// Sets the tolerance parameter to \e d.
	inline void SetAdaptiveTolerance(TScalar d) { m_AdaptiveTolerance = d; }

	/// Sets the initial step to \e d.
	inline void SetInitialStepMultiplier(TScalar d) { m_InitialStepMultiplier = d; }

	/// Set the initial spacing of the control points to \e d.
	inline void SetInitialCPSpacing(TScalar d) { m_InitialCPSpacing = d; }

	/// Initializes the control points by reading the file \e fn.
	inline void SetInitialCPPosition(std::string fn)
	{
		if (strlen(fn.c_str())) // null file name means no file to be read
		{
			m_ControlPoints = readMatrixDLM<TScalar>(fn.c_str());
			m_NumberOfCPs = m_ControlPoints.rows();
		}
	}

	/// Initializes the momentas by reading the file \e fn.
	inline void SetInitialMomentas(std::string fn)
	{
		if (strlen(fn.c_str())) // null file name means no file to be read
		{
			m_InitialMomentas = readMatrixDLM<TScalar>(fn.c_str());
		}
	}

	/// Sets the deformable objects of the target to \e obj
	/// (initialization of the number of objects at the same time).
	inline void SetTarget(DeformableMultiObjectType* obj) {
		m_Target = obj; m_NumberOfObjects = obj->GetNumberOfObjects();
	}

	/// Sets the deformable objects of the source to \e obj.
	inline void SetSource(DeformableMultiObjectType* obj) { m_Source = obj; }

	/// Returns the deformation.
	inline DiffeosType* GetDiffeos() const { return m_Def; }
	/// Sets the deformation to \e def.
	inline void SetDiffeos(DiffeosType* def) {m_Def = def;}

	/// Selects the Gradient Descent optimization method.
	inline void SetGradientDescent() { m_OptimizationMethod = GradientDescent; }
	/// Selects the ISTA algorithm.
	inline void SetISTA() { m_OptimizationMethod = ISTA; }
	/// Selects the FISTA algorithm.
	inline void SetFISTA() { m_OptimizationMethod = F_ISTA; }

	/// Freezes the control points.
	inline void SetFreezeCP() { m_freezeCP = true; }
	/// Unfreezes the control point.
	inline void UnsetFreezeCP() { m_freezeCP = false; }

	/// Sets the sparsity prior to \e d.
	inline void SetSparsityPrior(TScalar d) { m_SparsityPrior = d; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the list of deformed objects at final time.
	DeformableMultiObjectType* GetDeformedObject();

	/// Calls the method DeformableMultiObject::WriteFlow after updating the flow
	/// w.r.t. the new control points and momentas.
	void WriteFlow(std::vector<std::string>& name, std::vector<std::string>& extension);

	/// Performs the optimization method (according to SparseDiffeoMatcher::m_OptimizationMethod).
	void Update();



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Initializes the control points and momentas.
	/// \warning	If there are no control points, the GenerateInitialControlPoints() method will be
	///				automatically called.
	void InitializeControlPointsAndMomentas();

	/// Initializes the control points at the node of a regular lattice covering the data domain.
	/// Lattice step is given by SparseDiffeoMatcher::m_InitialCPSpacing.
	void GenerateInitialControlPoints();

	/**
	 *	\brief		Performs line search gradient descent.
	 *
	 *	\details	This method implements gradient descent with a line search strategy.
	 *  		\li	If m_ISTA = 0, gradient descent step with current stepsize is accepted if energy
	 *  			is decreased (\f$f(x^{n+1}) - f(x^n) < 0\f$) ;
	 *  		\li If m_ISTA = 1 and m_sparsityPrior = 0.0 (no \f$L^1\f$ penalty), gradient descent step
	 *  			with current stepsize is accepted if energy difference is such that :
	 *  			\f[
	 *  			f(x^{n+1}) - f(x^n) < -0.5 * \verb#stepsize# * ||\nabla_{x^n} f||^2 ~;
	 *  			\f]
	 *  		\li If m_ISTA = 1 and m_sparsityPrior > 0, gradient descent step with current stepsize
	 *				is accepted provided that :
	 *  			\f[
	 *  			f_2(x^{n+1}) - f_2(x^n) < \verb#_QdiffTerm#(x^n,x^{n+1}) ,
	 *  			\f]
	 *  			where \f$f_2\f$ denotes the quadratic part of the functional and momenta updates are soft-thresholded.
	 *
	 *  See details in : \n
	 *  A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems,
	 *  A. Beck & M. Teboulle, SIAM J. Imaging Sciences 2(1), 183--202 (2009)
	 *
	 *	\param[in, out]	X	Control points to be optimized or not (depending on
	 *						SparseDiffeoMatcher::m_freezeCP).
	 *	\param[in, out]	A	Momentas to be optimized.
	 */
	void GradientDescentAndISTA(VNLMatrixType&X, VNLMatrixType& A);

	/**
	 *	\brief		Performs FISTA.
	 *
	 *	\details	This optimization method combines the ISTA method (condition for accepting
	 *				gradient step) with an improved gradient descent scheme introduced by Nesterov :\n
	 *  		\li	If m_sparsityPrior = 0.0, the Nesterov scheme is applied as explained in :\n
	 *				A method of solving a convex programming problem with convergence rate O (1/k2),
	 *				Y. Nesterov, Soviet Mathematics Doklady 27 (2), 372-376 (1983) ;
	 *  		\li	If m_SparsityPrior > 0, the method implements the FISTA method as explained in :\n
	 *  			A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems,
	 *  			A. Beck & M. Teboulle, SIAM J. Imaging Sciences 2(1), 183--202 (2009) .
	 *
	 *	\param[in, out]	X	Control points positions (optimized or not, depending on
	 *						SparseDiffeoMatcher::m_freezeCP).
	 *	\param[in, out]	A	Momentas to be optimized.
	 */
	void FISTA(VNLMatrixType& X, VNLMatrixType& A);

	/// Compute value of the functional at a given step of the optimization.
	FunctionalValuesType* ComputeFunctional(VNLMatrixType& X, VNLMatrixType& A);

	/// Computes the gradient of the functional with respect to the positions of the control points \e X and the momentas \e A
	void ComputeGradient(VNLMatrixType& X, VNLMatrixType& A);

	/// Masks matrix rows. Rows whose index are masked (i.e. mask[i] = 1) are removed from matrix \e M
	VNLMatrixType _maskMatrix(const VNLMatrixType& M, const std::vector<int>& mask);

	/**
	 * \brief		Computes a gradient descent step for test values of the variables and step sizes
	 *
	 * \details		If there is no sparsity prior and if control points can move, this method computes :
	 *				\f[ x^{n+1} \leftarrow x^n - \tau \left(\nabla_{x}\functional\right) ,\f]
	 *				\f[ \alpha^{n+1} \leftarrow \alpha^n - \tau \left(\nabla_{\alpha}\functional\right), \f]
	 *				where \f$\tau\f$ is the current step-size of the gradient descent. \n
	 *				If the sparsity prior is positive, then the method SoftThresholdUpdate() is called
	 *				to update the momentas using a soft-thresholding function. If control points can't move, then \f$x^{n+1} \leftarrow x^n\f$.
	 *
	 * \param[out]	Xtest			Updated value of the control points \f$\left(x_i^{n+1}\right)_i\f$.
	 * \param[in]	X				Current value of the control points \f$\left(x_i^n\right)_i\f$.
	 * \param[out]	Atest			Updated value of the momentas \f$\left(\alpha_i^{n+1}\right)_i\f$.
	 * \param[in]	A				Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 * \param[in]	activeCPMask	Mask of active control points (those carrying non-zero momenta).
	 * \param[in]	step			Current stepsize.
	 */
	void GradientDescentStep(
			VNLMatrixType& Xtest, const VNLMatrixType& X,
			VNLMatrixType& Atest, const VNLMatrixType& A,
			std::vector<int>& activeCPMask,
			TScalar step);

	/**
	 * \brief		Soft-Threshold momentas updates.
	 *
	 * \details		This method soft-threshold the gradient descent step of the momentas. Has no effect if m_sparsityPrior = 0.
	 * 				The momentas are updated according to :
	 * 				\f[
	 * 				\alpha^{n+1}_{0,p} \leftarrow S_{\tau\sparsityprior}\left(\norm{\alpha^n_{0,k}-\tau\nabla_{\alpha}\functional}\right)
	 * 				\frac{\alpha^n_{0,k}-\tau\nabla_{\alpha}\functional}{\norm{\alpha^n_{0,k}-\tau\nabla_{\alpha}\functional}},
	 * 				\f]
	 * 				where \f$\tau\f$ is the current step-size of the gradient descent and \f$S\f$ is
	 * 				the thresholding function \f$S_{\lambda}(x)=\max(0,x-\lambda)+\min(0,x+\lambda)\f$.
	 *
	 * \param[in]	X		Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 * \param[in]	gradX	Gradient of the functional at current iteration.
	 * \param[in]	step	Descent step size for the momentas.
	 * \param[out]	mask	Mask of active control points.
	 * \return		The soft-thresholded momentas.
	 */
	VNLMatrixType SoftThresholdUpdate(const VNLMatrixType& X, VNLMatrixType& gradX, TScalar step, std::vector<int>& mask);

	/**
	 * \brief	Computes the condition for accepting gradient descent step for the F/ISTA methods.
	 *
	 * \details	F/ISTA methods require the current stepsize to satisfy
	 * 			\f$f(x^{n+1}) - f(x^n) < \verb#_QdiffTerm#(x^{n+1},x^{n})\f$ :
	 * 		\li	If m_sparsityPrior = 0.0, \f$\verb#_QdiffTerm# = -0.5 * stepsize * ||\nabla_{x^n} f||^2\f$ ;
	 * 		\li	If m_sparsityPrior > 0.0, _QdiffTerm() takes into account the soft-thresholded values of \f$x^n\f$.
	 *
	 * See details in \n
	 * A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems,
	 * A. Beck & M. Teboulle, SIAM J. Imaging Sciences 2(1), 183--202 (2009)
	 *
	 * \param[out]	Xtest	Updated value of the control points \f$\left(x_i^{n+1}\right)_i\f$.
	 * \param[in]	X		Current value of the control points \f$\left(x_i^n\right)_i\f$.
	 * \param[out]	Atest	Updated value of the momentas \f$\left(\alpha_i^{n+1}\right)_i\f$.
	 * \param[in]	A		Current value of the momentas \f$\left(\alpha_i^n\right)_i\f$.
	 * \param[in]	step	Current stepsize.
	 * \return		The value of QdiffTerm.
	 */
	TScalar _QdiffTerm(
			const VNLMatrixType& Xtest, const VNLMatrixType& X,
			const VNLMatrixType& Atest, const VNLMatrixType& A,
			TScalar step);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Initial spacing of the control points.
	TScalar m_InitialCPSpacing;

	/// Number of control points (and momentas).
	int m_NumberOfCPs;

	/// Maximum number of descent iterations.
	unsigned int m_MaxIterations;

	/// Maximum number of iterations in the line search procedure.
	unsigned int m_MaxLineSearchIterations;

	/// Shrink parameter for the line search.
	TScalar m_AdaptiveShrink;

	/// Expand parameter for the line search.
	TScalar m_AdaptiveExpand;

	/// The algorithm stops when \f$F(end-1)-F(end) < \verb#m_AdaptiveTolerance# * \left( F(0)-F(end) \right)\f$
	TScalar m_AdaptiveTolerance;

	/// FISTA parameter. Can't increase stepsize during NbIterFreeze iterations.
	static const unsigned int NbIterFreeze = 3;

	/// Initial step multiplier.
	TScalar m_InitialStepMultiplier;

	/// Coefficient added to the \f$L^1\f$ penalty.
	TScalar m_SparsityPrior;

	/// Type of optimization method used.
	OptimizationMethodType m_OptimizationMethod;

	/// Boolean which indicates if the control points can move or not.
	bool m_freezeCP;

	/// Gradient of the functional w.r.t. the positions
	/// (Size of \f$\nabla_{x}\functional\f$: NumberOfCPs x Dimension).
	VNLMatrixType m_GradPos;

	/// Gradient of the functional w.r.t. the momentas
	/// (Size of \f$\nabla_{\alpha}\functional\f$: NumberOfCPs x Dimension).
	VNLMatrixType m_GradMom;

	/// Pointer on the collection of deformable objects of the target.
	DeformableMultiObjectType* m_Target;

	/// Pointer on the collection of deformable objects of the source.
	DeformableMultiObjectType* m_Source;

	/// Number of deformable objects in the target and the source.
	int m_NumberOfObjects;

	/// Pointer to a deformation.
	DiffeosType* m_Def;

	/// Matrix containing the coordinates of the control points (Size : NumberOfCPs x Dimension).
	VNLMatrixType m_ControlPoints;

	/// Matrix containing the coordinates of the momentas (Size : NumberOfCPs x Dimension).
	VNLMatrixType m_InitialMomentas;

	/// Vector containing a history of the the values of the functional during optimization method.
	std::vector< FunctionalValuesType* > m_ValuesHistory;

}; /* class SparseDiffeoMatcher */


#ifndef MU_MANUAL_INSTANTIATION
#include "SparseDiffeoMatcher.txx"
#endif


#endif /* _SparseDiffeoMatcher_h */
