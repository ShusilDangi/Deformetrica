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

#ifndef _GaussFunction_h
#define _GaussFunction_h

#include "RadialFunction.h"

/**
 *  \brief      Gaussian function.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.1
 *
 *  \details    The GaussianFunction class inherited from RadialFunction implements the operations for
 *              a radial with a Gaussian kernel.
 */
template < typename TYPE >
class GaussFunction : public RadialFunction<TYPE>
{

    TYPE Sigma, ooSigma2, ooSigma4;

public:

    GaussFunction() { }

    /// Contructor with parameter \e Sigma.
    GaussFunction(TYPE sigma)
    {
        Sigma = sigma;
        ooSigma2 = 1.0/(Sigma*Sigma);
        ooSigma4 = 1.0/(Sigma*Sigma*Sigma*Sigma);
    }

    __device__ __forceinline__ TYPE Eval(TYPE r2)
    {
        return exp(-r2*ooSigma2);
    }

    __device__ __forceinline__ TYPE Diff(TYPE r2)
    {
        return - ooSigma2 * exp(-r2*ooSigma2);
    }

    __device__ __forceinline__ TYPE Diff2(TYPE r2)
    {
        return ooSigma4 * exp(-r2*ooSigma2);
    }

    __device__ __forceinline__ void DiffDiff2(TYPE r2, TYPE* d1, TYPE* d2)
    {
        *d1 = - ooSigma2 * exp(-r2*ooSigma2);
        *d2 = - ooSigma2 * *d1;
    }

}; /* class GaussFunction */

#endif /* _GaussFunction_h */
