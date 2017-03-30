/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _Diffeos_h
#define _Diffeos_h

#include "AbstractDeformations.h"

/**
 *	\brief 		Standard diffeomorphisms.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The Diffeos class inherited from AbstractDeformations represents the standard
 *				deformation which is usually employed in minimization problems such as registration
 *				or atlas construction. Deformations are encoded by control points and momentum vectors attached to them, which move in time accroding a to an Hamiltonian set of equations.\n \n
 */
template <class TScalar, unsigned int Dimension>
class Diffeos : public AbstractDeformations<TScalar, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Abstract deformation type.
	typedef AbstractDeformations<TScalar, Dimension> Superclass;

	/// Vector type.
	typedef typename Superclass::VNLVectorType VNLVectorType;
	/// Matrix type.
	typedef typename Superclass::VNLMatrixType VNLMatrixType;
	/// List of matrices type.
	typedef typename Superclass::VNLMatrixList VNLMatrixList;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	Diffeos();

	/// Copy constructor.
	Diffeos(const Diffeos& other);

	virtual Diffeos* Clone() { return new Diffeos(*this); };

	virtual ~Diffeos();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void ReverseFlow();



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual void Shoot();



}; /* class Diffeos */


#ifndef MU_MANUAL_INSTANTIATION
#include "Diffeos.txx"
#endif


#endif /* _Diffeos_h */

