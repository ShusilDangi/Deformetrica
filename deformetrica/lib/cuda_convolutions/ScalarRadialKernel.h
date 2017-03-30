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

#ifndef _ScalarRadialKernel_h
#define _ScalarRadialKernel_h

#include "RadialFunction.h"



/**
 *  \brief      Scalar radial kernel.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.1
 *
 *  \details    The ScalarRadialKernel class implements standard operations with kernels
 *              when dealing with convolutions. More particulary, the kernels are
 *              translation-invariant isotropic scalar kernels. \n\n
 *
 *              In the following, \f$r\f$ will denote \f$||x_i-y_j||^2\f$.
 */
template < typename TYPE, int DIMPOINT, int DIMVECT, class RADIAL_FUN >
class ScalarRadialKernel
{

    /// Our radial function \f$f(r)\f$.
    RADIAL_FUN Rfun;

public:

    /// Constructor with initialization of the radial function to \e rfun.
    ScalarRadialKernel(RADIAL_FUN rfun)
    {
        Rfun = rfun;
    }

    /// Computes \f$\gamma_i += f(r)\beta_j\f$.
    __device__ __forceinline__ void Eval(TYPE* gammai, TYPE* xi, TYPE* yj, TYPE* betaj)
    {
        TYPE r2 = 0.0f;
        TYPE temp;
        for(int k=0; k<DIMPOINT; k++)
        {
            temp =  yj[k]-xi[k];
            r2 += temp*temp;
        }
        TYPE s = Rfun.Eval(r2);
        for(int k=0; k<DIMVECT; k++)
            gammai[k] += s * betaj[k];
    }

    /// Computes \f$\gamma_i += 2 (\alpha_i^T \beta_j) f'(r)\beta_j\f$.
    __device__ __forceinline__ void Grad1(TYPE* gammai, TYPE* alphai, TYPE* xi, TYPE* yj, TYPE* betaj)
    {
        TYPE r2 = 0.0f, sga = 0.0f;
        TYPE xmy[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            xmy[k] =  xi[k]-yj[k];
            r2 += xmy[k]*xmy[k];
        }
        for(int k=0; k<DIMVECT; k++)
            sga += betaj[k]*alphai[k];
        TYPE s = 2.0 * sga * Rfun.Diff(r2);
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] += s * xmy[k];
    }

    /// Computes \f$\gamma_i += 2 (\beta_j^T \alpha_i + \beta_i^T \alpha_j) f'(r)\beta_j\f$.
    __device__ __forceinline__ void Grad(TYPE* gammai, TYPE* xi, TYPE* xj, TYPE* alphai, TYPE* alphaj, TYPE* betai, TYPE* betaj)
    {
        TYPE r2 = 0.0f, sga = 0.0f;
        TYPE ximxj[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            ximxj[k] =  xi[k]-xj[k];
            r2 += ximxj[k]*ximxj[k];
        }
        for(int k=0; k<DIMVECT; k++)
            sga += betaj[k]*alphai[k] + betai[k]*alphaj[k];
        TYPE s = 2.0 * sga * Rfun.Diff(r2);
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] += s * ximxj[k];
    }

    /// Computes \f$\gamma_i += 4 (\beta_j^T \beta_i) f'(r) (\eta_i-\eta_j) + 8 (\beta_i^T \beta_j) f"(r) \left[(x_i-x_j)^T(\eta_i-\eta_j)\right](x_i-x_j)\f$.
    __device__ __forceinline__ void GradDiff(TYPE* gammai, TYPE* xi, TYPE* xj, TYPE* betai, TYPE* betaj, TYPE* etai, TYPE* etaj)
    {
        TYPE r2 = 0.0f, bidbj = 0.0f, dotex = 0.0f;
        TYPE ximxj[DIMPOINT], eimej[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            ximxj[k] =  xi[k]-xj[k];
            r2 += ximxj[k]*ximxj[k];
            eimej[k] =  etai[k]-etaj[k];
            dotex += ximxj[k]*eimej[k];
        }
        for(int k=0; k<DIMVECT; k++)
            bidbj += betai[k]*betaj[k];
        TYPE d1, d2;
        Rfun.DiffDiff2(r2,&d1,&d2);
        d1 *= 4 * bidbj;
        d2 *= 8 * bidbj * dotex;
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] += d1 * eimej[k] + d2 * ximxj[k];
    }

    /// Computes \f$\gamma_i += 2 f'(r) \left[(x_i-x_j)^T(\eta_i-\eta_j)\right] \beta_j\f$.
    __device__ __forceinline__ void Diff(TYPE* gammai, TYPE* xi, TYPE* xj, TYPE* betaj, TYPE* etai, TYPE* etaj)
    {
        TYPE r2 = 0.0f, dotex = 0.0f;
        TYPE ximxj[DIMPOINT], eimej[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            ximxj[k] =  xi[k]-xj[k];
            r2 += ximxj[k]*ximxj[k];
            eimej[k] =  etai[k]-etaj[k];
            dotex += ximxj[k]*eimej[k];
        }
        TYPE s = Rfun.Diff(r2) * 2.0 * dotex;
        for(int k=0; k<DIMVECT; k++)
            gammai[k] += s * betaj[k];
    }



}; /* class ScalarRadialKernel */



#endif /* _ScalarRadialKernel_h */
