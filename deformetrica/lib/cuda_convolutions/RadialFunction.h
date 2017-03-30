/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _RadialFunction_h
#define _RadialFunction_h

#include <cuda.h>

/**
 *  \brief      Radial function.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.1
 *
 *  \details    The RadialFunction class defines the different evaluation methods when working with kernels of the form :
 *              \f[
 *              K(x,y) = f(|x-y|)) = f(u).
 *              \f]
 *              In other words, the radial functions are translation-invariant isotropic kernels.\n\n
 *
 *              See child classes for examples of radial functions.
 */
template < typename TYPE >
class RadialFunction
{

public:

	/// Computes \f$f(u)\f$.
    virtual __device__ __forceinline__ TYPE Eval(TYPE u) = 0;

    /// Computes the derivative of \f$f(u)\f$.
    virtual __device__ __forceinline__ TYPE Diff(TYPE u) = 0;

    /// Computes the second derivative of \f$f(u)\f$.
    virtual __device__ __forceinline__ TYPE Diff2(TYPE u) = 0;

    /// Stores in \e d1 (resp. \e d2) the first (resp. second) derivative of \f$f(u)\f$.
    virtual __device__ __forceinline__ void DiffDiff2(TYPE u, TYPE* d1, TYPE* d2)
    {
        *d1 = Diff(u);
        *d2 = Diff2(u);
    }

}; /* class RadialFunction */


#endif /* _RadialFunction_h */
