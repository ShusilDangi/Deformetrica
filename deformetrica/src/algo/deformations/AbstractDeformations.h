/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _AbstractDeformations_h
#define _AbstractDeformations_h

#include "vnl/vnl_matrix.h"

#include <vector>


/**
 *	\brief 		Abstract deformations.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The AbstractDeformations class encodes diffeomorphic deformations of the ambient 2D or 3D space, which hen deforms objects embedded into it.\n\n
 */
template <class TScalar, unsigned int Dimension>
class AbstractDeformations
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Possible type of deformations.
	typedef enum
	{
		null,			/*!< Default value. */
		Diffeos			/*!< Standard deformation (see Diffeos). */
	} DeformationsType;

	/// Vector type.
	typedef vnl_vector<TScalar> VNLVectorType;
	/// Matrix type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// List of matrices type.
	typedef std::vector<VNLMatrixType> VNLMatrixList;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	AbstractDeformations();
	///	Copy constructor.
	AbstractDeformations(const AbstractDeformations& other);

	///	Makes a copy of the object.
	virtual AbstractDeformations* Clone() = 0;

	virtual ~AbstractDeformations();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns true if it is of Diffeos type, false otherwise.
	inline bool isADiffeo() const { return (m_Type == Diffeos); }
	/// Sets type of the deformation to Diffeos.
	inline void SetDiffeosType() { m_Type = Diffeos; }

	/// Returns the number of time points between \f$t_0\f$ and \f$t_n\f$.
	inline unsigned int GetNumberOfTimePoints() const { return m_NumberOfTimePoints; }
	/// Sets the number of time points between \f$t_0\f$ and \f$t_n\f$.
	inline void SetNumberOfTimePoints(unsigned int n) { m_NumberOfTimePoints = n; m_Modified = true;}

	///	Returns the size of the kernel.
	inline TScalar GetKernelWidth() const { return m_KernelWidth; }
	/// Sets the size of the kernel to \e kernelWidth.
	inline void SetKernelWidth(TScalar kernelWidth) { m_KernelWidth = kernelWidth; m_Modified = true; }

	/// Sets the initial positions to \e X.
	inline void SetStartPositions(const VNLMatrixType& X) { m_StartPositions = X; m_Modified = true; }

	/// Sets the initial momentas to \e A.
	inline void SetStartMomentas(const VNLMatrixType& A) { m_StartMomentas = A; m_Modified = true; }

	/// Returns the list of the different positions between \f$t_0\f$ and \f$t_n\f$.
	inline VNLMatrixList GetTrajectoryPositions() const { return m_PositionsT; }

	/// Returns the list of the different momentas between \f$t_0\f$ and \f$t_n\f$.
	inline VNLMatrixList GetTrajectoryMomentas() const { return m_MomentasT; }

	/// Returns true if use improved Euler's method, false otherwise.
	inline bool ImprovedEuler() const { return m_UseImprovedEuler; }
	/// Sets standard Euler's method.
	void UseStandardEuler() { m_UseImprovedEuler = false; m_Modified = true; }
	/// Sets improved Euler's method.
	void UseImprovedEuler() { m_UseImprovedEuler = true; m_Modified = true; }

	/// Returns the data domain.
	inline VNLMatrixType GetDataDomain() const { return m_DataDomain; }
	/// Sets the data domain to \e domain.
	inline void SetDataDomain(VNLMatrixType& domain) { m_DataDomain = domain; }

	/// Returns true if any point is out of the bounding box, false otherwise.
	inline bool OutOfBox() const { return m_OutOfBox; }

	/// Returns the padding factor.
	inline TScalar GetPaddingFactor() const { return m_PaddingFactor; }
	/// Sets the padding factor to \e paddingFactor.
	inline void SetPaddingFactor(TScalar paddingFactor) { m_PaddingFactor = paddingFactor; }

	/// Returns true if one of the parameters of the deformation has changed, false otherwise.
	inline bool IsModified() const { return m_Modified; }
	/// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
	inline void SetModified() { m_Modified = true; }
	/// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
	inline void UnsetModified() { m_Modified = false; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 *	\brief		Checks if a set of points is out of box or not.
	 *
	 *	\details	This function enables to see if the coordinates of X at time t is out of box or not.
	 *
	 *	\param[in]	X	List of matrices containing coordinates at different time steps.
	 *	\param[in]	t	Time index for X[t].
	 *	\return		True if X[t] is out of box, false otherwise.
	 */
	bool CheckBoundingBox(VNLMatrixList& X, int t);

	/**
	 *	\brief		Updates the trajectory.
	 *
	 *	\details	If the class has been modified, this method computes trajectory thanks to the
	 *				Shoot() method.
	 */
	virtual void Update();

	/// Reverse the direction of the flow.
	/// \warning   Be careful, this is not the flow of the inverse map.
	virtual void ReverseFlow() = 0;



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 *	\brief		Computes the trajectory.
	 *
	 *	\details	According to the type of deformation, this method computes the trajectory between
	 *				\f$ t_0 \f$ and \f$ t_n \f$ of the initial momentas and positions.
	 */
	virtual void Shoot() = 0;

	/**
	 *	\brief		Initializes the bounding box.
	 *
	 *	\details	This method initialize the AbstractDeformations::m_BoundingBox attribute according to
	 *				the data domain, the padding factor and the size of the kernel.
	 */
	void InitBoundingBox();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Type of the deformation.
	DeformationsType m_Type;

	///	Number of time points between \f$ t_0 \f$ and \f$ t_n \f$.
	unsigned int m_NumberOfTimePoints;

	///	Matrix containing initial positions (Size : N x Dimension).
	VNLMatrixType m_StartPositions;

	///	Matrix containing initial momentas (Size : N x Dimension).
	VNLMatrixType m_StartMomentas;

	///	Matrix containing the min (resp. the max) of the initial positions at first (resp. second) line.
	VNLMatrixType m_DataDomain;

	///	Size of the kernel.
	TScalar m_KernelWidth;

	///	List of m_NumberOfTimePoints matrices containing positions.
	VNLMatrixList m_PositionsT;

	///	List of m_NumberOfTimePoints matrices containing momentas.
	VNLMatrixList m_MomentasT;

	/// Boolean which indicates if we use improved Euler's method or not.
	bool m_UseImprovedEuler;

	///	Boolean which avoids computing the trajectory (via Update()) if no parameter has changed.
	bool m_Modified;

	///	Boolean which prevents computations if any point (e.g. the coordinates of trajectory) is outside the bounding box.
	bool m_OutOfBox;

	///	Multiplier coefficient for the creation of the bounding box ("boundingBox = dataDomain +/- paddingFactor*kernelWidth/2").
	TScalar m_PaddingFactor;

	///	Box where any trajectory must not exit (Size : 2 x Dimension with the min (resp. the max) at first (resp. second) line).
	VNLMatrixType m_BoundingBox;



}; /* class AbstractDeformations */


#ifndef MU_MANUAL_INSTANTIATION
#include "AbstractDeformations.txx"
#endif


#endif /* _AbstractDeformations_h */
