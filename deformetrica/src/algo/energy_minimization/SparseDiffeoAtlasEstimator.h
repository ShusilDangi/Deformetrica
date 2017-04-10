/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoAtlasEstimator_h
#define _SparseDiffeoAtlasEstimator_h

#include "itkVector.h"

#include "itkSimpleFastMutexLock.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "readMatrixDLM.txx"

#include "Diffeos.h"
#include "DeformableMultiObject.h"

#include "AtlasFunctionalValues.h"

#include <vector>

/**
 *	\brief 		Atlas construction class
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The SparseDiffeoAtlasEstimator class enables to build atlases from collections of objects configurations.\n \n
 *				Given a population of M objects and N subjects, it estimates a common template complex T made of M template objects
 *				and its deformations to each shape complex in the population. \n \n
 */
template <class TScalar, unsigned int Dimension>
class SparseDiffeoAtlasEstimator
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Vector type (std).
	typedef std::vector<TScalar> VectorType;
	/// Vector type
	typedef vnl_vector<TScalar> VNLVectorType;
	/// Matrix type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// List of matrices type.
	typedef std::vector<VNLMatrixType> VNLMatrixList;

	/// Deformable object type.
	typedef DeformableObject<TScalar, Dimension> DeformableObjectType;
	/// List of deformable objects type.
	typedef std::vector<DeformableObjectType*> DeformableObjectList;
	/// Deformable multi-object type.
	typedef DeformableMultiObject<TScalar, Dimension> DeformableMultiObjectType;
	/// List of deformable multi-objects type.
	typedef std::vector<DeformableMultiObjectType*> DeformableMultiObjectList;

	/// Deformation type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;

	/// Functional values type.
	typedef AtlasFunctionalValues<TScalar> FunctionalValuesType;

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

	SparseDiffeoAtlasEstimator();

	~SparseDiffeoAtlasEstimator();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the size of the smoothing kernel to \e d. See ConvolveGradTemplate().
	inline void SetSmoothingKernelWidth(TScalar d) { m_SmoothingKernelWidth = d; }

	/// Sets the maximum of descent iterations to \e n.
	void SetMaxIterations(unsigned int n);

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
			m_InitialMomentas = readMultipleMatrixDLM<TScalar>(fn.c_str());
	}

	/// Set the sparsity prior to \e d.
	inline void SetSparsityPrior(TScalar d) { m_SparsityPrior = d; }

	/// Sets the list of targets to \e obj.
	inline void SetTargetList(DeformableMultiObjectList& obj)
	{
		m_TargetList = obj;
		m_NumberOfSubjects = obj.size();
		m_NumberOfObjects = obj[0]->GetNumberOfObjects();
	}

	/// Returns the template.
	DeformableMultiObjectType* GetTemplate();
	/// Sets the template to \e obj.
	inline void SetTemplate(DeformableMultiObjectType* obj)	{ m_Template = obj; }

	/// Returns the deformation.
	inline DiffeosType* GetDiffeos() const { return m_Def; }
	/// Sets the deformation to \e def.
	inline void SetDiffeos(DiffeosType* def) { m_Def = def; }

	/// Selects the Gradient Descent optimization method.
	inline void SetGradientDescent() { m_OptimizationMethod = GradientDescent; }
	/// Selects the ISTA algorithm.
	inline void SetISTA() { m_OptimizationMethod = ISTA; }
	/// Selects the FISTA algorithm.
	inline void SetFISTA() { m_OptimizationMethod = F_ISTA; }

	/// Disables optimization in control point positions.
	inline void SetFreezeCP() { m_freezeCP = true; }
	/// Enables optimization in control point positions.
	inline void UnsetFreezeCP() { m_freezeCP = false; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the weights for the \e s-th subject according to the paramsDiffeo file
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

	/// Returns the list of deformed objects at final time.
	DeformableMultiObjectList GetDeformedObjects();

	/// Calls the method DeformableMultiObject::WriteFlow for each subject
	/// after updating the flow w.r.t. the new control points and momentas.
	void WriteFlow(std::vector<std::string>& name, std::vector<std::string>& extension);

	/// Set the number of threads to \e n.
	void SetNumberOfThreads(unsigned int n) { m_NumberOfThreads = n; }

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
	/// Lattice step is given by SparseDiffeoAtlasEstimator::m_InitialCPSpacing.
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
	 *						SparseDiffeoAtlasEstimator::m_freezeCP).
	 *	\param[in, out]	A	Momentas to be optimized.
	 *	\param[in, out]	T	Template to be optimized.
	 */
	void GradientDescentAndISTA(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T);

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
	 *						SparseDiffeoAtlasEstimator::m_freezeCP).
	 *	\param[in, out]	A	Momentas to be optimized.
	 *	\param[in, out]	T	Template to be optimized.
	 */
	void FISTA(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T);

	/// Compute value of the functional at a given step of the optimization.
	FunctionalValuesType* ComputeFunctional(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T);

	/**
	 * \brief		Compute terms of the functional corresponding to a given subject.
	 *
	 * \details		This method performs the functional according to the given control points and momentas
	 * 				for a specific subject.
	 *
	 * \param[in]	X0			Current positions of the control points.
	 * \param[in]	Mom0		Current values of the momentas for the given subject.
	 * \param[in]	targetObj	Collection of target objects of the given subject.
	 * \param[out]	funcValues	Computes and stores each term of the functional.
	 */
	bool ComputeFunctionalSubject(VNLMatrixType& X0, VNLMatrixType& Mom0,
			DeformableMultiObjectType* targetObj, std::vector<TScalar>& funcValues);

	/// Computes the gradient of the functional with respect to the positions of the control points \e X,
	/// the momentas \e A and the template objects \e T
	void ComputeGradient(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T);

	/**
	 * \brief		Computes the gradient of the terms in the functional corresponding to a given subject.
	 *
	 * \details		This method the gradient of the functional w.r.t. control point positions, momenta
	 * 				and template shapes for a given subject
	 *
	 * \param[in]	X			Current positions of the control points.
	 * \param[in]	A			Current values of the momentas for the given subject.
	 * \param[out]	gradX		Gradient of the functional w.r.t. control points.
	 * \param[out]	gradA		Gradient of the functional w.r.t. momentas.
	 * \param[out]	gradT		Gradient of the functional w.r.t. template.
	 * \param[in]	targetObj	Collection of target objects of the given subject
	 */
	void ComputeGradientSubject(VNLMatrixType& X, VNLMatrixType& A,
			VNLMatrixType& gradX, VNLMatrixType& gradA, VNLMatrixList& gradT,
			DeformableMultiObjectType* targetObj);

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
	 * \param[out]	Ttest			Updated value of the template \f$\left(x_0\right)^{n+1}\f$.
	 * \param[in]	T				Current value of the template \f$\left(x_0\right)^n\f$.
	 * \param[in]	activeCPMask	Mask of active control points (those carrying non-zero momenta).
	 * \param[in]	stepXA			Step-size for momentas and control point positions update.
	 * \param[in]	stepT			Step-size for template update.
	 */
	void GradientDescentStep(
			VNLMatrixType& Xtest, const VNLMatrixType& X,
			VNLMatrixList& Atest, const VNLMatrixList& A,
			VNLMatrixList& Ttest, const VNLMatrixList& T,
			std::vector<int>& activeCPMask,
			TScalar stepXA, TScalar stepT);

	/// Masks matrix rows. Rows whose index are masked (i.e. mask[i] = 1) are removed from matrix \e M.
	VNLMatrixType _maskMatrix(const VNLMatrixType& M, const std::vector<int>& mask);

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
	 * \param[out]	Ttest	Updated value of the template \f$\left(x_0\right)^{n+1}\f$.
	 * \param[in]	T		Current value of the template \f$\left(x_0\right)^n\f$.
	 * \param[in]	stepXA	Current stepsize for the control points and the momentas.
	 * \param[in]	stepT	Current stepsize for the template.
	 * \return		The value of QdiffTerm.
	 */
	TScalar _QdiffTerm(
			const VNLMatrixType& Xtest, const VNLMatrixType& X,
			const VNLMatrixList& Atest, const VNLMatrixList& A,
			const VNLMatrixList& Ttest, const VNLMatrixList& T,
			TScalar stepXA, TScalar stepT);

	/**
	 * \brief	Convolves the gradient w.r.t. template meshes with a Gaussian kernel.
	 *
	 * \details	This method convolves the gradient w.r.t. template meshes with a Gaussian kernel
	 * 			with standard deviation SparseDiffeoAtlasEstimator::m_SmoothingKernelWidth.
	 * 			It guarantees that the topology of the mesh is preserved along the iterations
	 * 			of the gradient descent (no self-intersection).
	 * 		\li	It has no effect if m_SmoothingKernelWidth = 0.0 ;
	 * 		\li	It has no effect on non-mesh data (landmarks, point cloud, images) .
	 *
	 * \param[in]	temp gradient of the functional w.r.t. positions of the vertices of template meshes.
	 */
	void ConvolveGradTemplate(VNLMatrixList& temp);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Initial spacing of the control points (needed if they are automatically generated before optimization).
	TScalar m_InitialCPSpacing;

	/// Number of control points (and momentas).
	int m_NumberOfCPs;

	/// Type of optimization method used.
	OptimizationMethodType m_OptimizationMethod;

	/// Size of the smoothing kernel used to smooth the gradient w.r.t. vertex position of template meshes.
	/// See ConvolveGradTemplate().
	TScalar m_SmoothingKernelWidth;

	/// Maximum number of descent iterations.
	unsigned int m_MaxIterations;

	/// Maximum number of iterations in the line search procedure.
	unsigned int m_MaxLineSearchIterations;

	/// Shrink parameter for the line search.
	TScalar m_AdaptiveShrink;
	/// Expand parameter for the line search.
	TScalar m_AdaptiveExpand;
	/// The algorithm stops when \f$F(end-1)-F(end) < \verb#m_AdaptiveTolerance# * \left( F(0)-F(end) \right)\f$.
	TScalar m_AdaptiveTolerance;

	/// FISTA parameter. Can't increase stepsize during NbIterFreeze iterations.
	static const unsigned int NbIterFreeze = 3;

	/// Initial step multiplier.
	TScalar m_InitialStepMultiplier;

	/// Weight of the \f$L^1\f$ penalty in the functional.
	TScalar m_SparsityPrior;

	/// Weights for each subject.
	VectorType m_Weights;

	/// Boolean which indicates if use the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) or not.
	bool m_UseFISTA;

	/// Indicates if control points are optimized or not.
	bool m_freezeCP;

	/// Gradient of the functional w.r.t. the CP position.
	VNLMatrixType m_GradPos;

	/// Gradient of the functional w.r.t. momentas
	VNLMatrixList m_GradMom;
	
	/// Gradient w.r.t. template objects (using a \f$L^2\f$ metric)
	VNLMatrixList m_GradTemplate_L2;
	/// Gradient w.r.t. template objects (using a Sobolev metric), i.e. the convolution of the \f$L^2\f$
	/// gradient with a Gaussian kernel. See ConvolveGradTemplate() for details.
	VNLMatrixList m_GradTemplate_Sob;

	/// \f$\nbsubj\f$ collections of deformable target objects.
	DeformableMultiObjectList m_TargetList;
	/// Pointer on the collection of deformable target objects.
	DeformableMultiObjectType* m_Template;

	/// Number of deformable objects of the target and the source.
	int m_NumberOfObjects;
	/// Number of subjects.
	int m_NumberOfSubjects;

	/// Pointer to a deformation.
	DiffeosType* m_Def;

	/// Matrix of the coordinates of the control points (Size : NumberOfCPs x Dimension).
	VNLMatrixType m_ControlPoints;
	/// List of matrices of the coordinates of the momenta (List size: NumberOfSubjects, Matrices size: NumberOfCPs x Dimension).
	VNLMatrixList m_InitialMomentas;
	/// List of matrices of the template point coordinates (List size: NumberOfObjects, Matrices size: NumberOfPoints x Dimension).
	VNLMatrixList m_TemplateData;

	/// Vector containing a history of the the values of the functional during optimization method.
	std::vector< FunctionalValuesType* > m_ValuesHistory;


protected:

	/// \cond HIDE_FOR_DOXYGEN

	FunctionalValuesType* ThreadedComputeFunctional(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T);

	void ThreadedComputeGradient(VNLMatrixType& X, VNLMatrixList& A, VNLMatrixList& T);

	TScalar _ThreadedQdiffTerm(
			const VNLMatrixType& Xtest, const VNLMatrixType& X,
			const VNLMatrixList& Atest, const VNLMatrixList& A,
			const VNLMatrixList& Ttest, const VNLMatrixList& T,
			TScalar stepXA, TScalar stepT);

	void ThreadedGradientDescentStep(
			VNLMatrixType& Xtest, const VNLMatrixType& X,
			VNLMatrixList& Atest, const VNLMatrixList& A,
			VNLMatrixList& Ttest, const VNLMatrixList& T,
			std::vector<int>& activeCPMask,
			TScalar stepXA, TScalar stepT);

	static ITK_THREAD_RETURN_TYPE _functionalThread(void* arg);

	static ITK_THREAD_RETURN_TYPE _gradientThread(void* arg);

	static ITK_THREAD_RETURN_TYPE _QdiffThread(void* arg);

	static ITK_THREAD_RETURN_TYPE _descentStepThread(void* arg);

	unsigned int m_NumberOfThreads;

	itk::SimpleFastMutexLock m_Mutex;

	unsigned int m_MT_SubjectCounter;

	FunctionalValuesType* m_MT_FuncValues;

	VNLMatrixType m_MT_X0;
	VNLMatrixList m_MT_Mom0;

	VNLMatrixList m_MT_Qdiff_Atest;
	VNLMatrixList m_MT_Qdiff_A;
	TScalar m_MT_Qdiff_Value;

	TScalar m_MT_StepXA;

	std::vector<int> m_MT_ActiveCPMask;

	/// \endcond


}; /* class SparseDiffeoAtlasEstimator */


#ifndef MU_MANUAL_INSTANTIATION
#include "SparseDiffeoAtlasEstimator.txx"
#endif


#endif /* _SparseDiffeoAtlasEstimator_h */
