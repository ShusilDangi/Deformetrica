#ifndef _FastGauss_txx
#define _FastGauss_txx

// Implementation for Fast Gauss Transform

#include "FastGauss.h"

#include <climits>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define FGT_MANAGE_MEMORY 1

////////////////////////////////////////////////////////////////////
// Constructor 1
//
// PURPOSE                                                    
// -------   
// Initialize the class without the parameter tol. 
// Allocate memory for the fast Gauss algorithm
// Determine the memory space needed. 
// Compute the coefficients of the Taylor expansions.
//
// PARAMETERS                                                      
// ----------
// dim		--> Dimension of problem
// u		--> Weight matrix with du rows and ns columns
// du		--> Row dimension of matrix u (default = 1)
// x		--> Source matrix with dim rows and ns columns
// ns		--> Column dimension of matrix x (no. of sources)
// h		--> Bandwidth
// p		--> Order of truncation
// K		--> Number of centers
// e		--> Ratio of far field
//
// OUTPUTS
// -------
// xc		<--	center of the sources
// A_k		<-- Coefficients of the Taylor expansion with wts u
////////////////////////////////////////////////////////////////////

template <class TScalar>
FastGauss<TScalar>
::FastGauss(int dim, int dw, const TScalar *u, const TScalar *x, int ns, 
		TScalar h, int p, int K, TScalar e, TScalar* ptr)
		{
	// interface:
	Dim = dim;
	NSources = ns;	
	DWeights = dw;
	ptRadius = ptr;
	bandwid = h;
	pterms = p;
	K_N = K;
	ratio_far = e;

	// Make this class manage its own memory?
#if FGT_MANAGE_MEMORY
	pWeights = new TScalar[ns*dw];
	for (int i = 0; i < (ns*dw); i++)
		pWeights[i] = u[i];
	pSources = new TScalar[ns*dim];
	for (int i = 0; i < (ns*dim); i++)
		pSources[i] = x[i];
#else
	pWeights = (TScalar*)u;
	pSources = (TScalar*)x;
#endif

	pd = nchoosek(pterms+Dim-1, Dim); //pd = C_dim^(dim+pterms-1)

	A_k = new TScalar[pd*K_N*DWeights];

	xc = new TScalar[K_N*Dim];
	xcind = new int[K_N];
	indx = new int[NSources];
	xheads = new int[K_N];
	xboxsz = new int[K_N];

	// Divide the source space into K_N parts using K-center algorithm
	rx2 = dkcenter(Dim,NSources,K_N,pSources,xcind,indx,xboxsz);
	*ptRadius = rx2;  // return radius of clustering

	// computer the center of the sources
	Compute_Centers(xc,pSources,Dim,NSources,indx,K_N,xboxsz);

	// Build hashing tables for the nearest neighbor searching

	TaylorExpansion(Dim, pterms);
}


template <class TScalar>
FastGauss<TScalar>
::FastGauss(const FastGauss& o)
 {
	// interface:
	Dim = o.Dim;
	NSources = o.NSources;
	DWeights = o.DWeights;
	bandwid = o.bandwid;
	pterms = o.pterms;
	K_N = o.K_N;
	ratio_far = o.ratio_far;

	// Make this class manage its own memory?
#if FGT_MANAGE_MEMORY
	pWeights = new TScalar[NSources*DWeights];
	for (int i = 0; i < (NSources*DWeights); i++)
		pWeights[i] = o.pWeights[i];
	pSources = new TScalar[NSources*Dim];
	for (int i = 0; i < (NSources*Dim); i++)
		pSources[i] = o.pSources[i];
#else
	pWeights = o.pWeights;
	pSources = o.pSources;
#endif

	pd = o.pd;

	A_k = new TScalar[pd*K_N*DWeights];
	for (int i = 0; i < pd*K_N*DWeights; i++)
		A_k[i] = o.A_k[i];

	xc = new TScalar[K_N*Dim];
	for (int i = 0; i < K_N*Dim; i++)
		xc[i] = o.xc[i];

	xcind = new int[K_N];
	for (int i = 0; i < K_N; i++)
		xcind[i] = o.xcind[i];

	indx = new int[NSources];
	for (int i = 0; i < NSources; i++)
		indx[i] = o.indx[i];

	xheads = new int[K_N];
	for (int i = 0; i < K_N; i++)
		xheads[i] = o.xheads[i];

	xboxsz = new int[K_N];
	for (int i = 0; i < K_N; i++)
		xboxsz[i] = o.xboxsz[i];

	rx2 = o.rx2;
	*ptRadius = rx2;  // return radius of clustering

 }



template <class TScalar>
FastGauss<TScalar>
::~FastGauss()
{
	delete []A_k;
	delete []xc;
	delete []xcind;
	delete []indx;
	delete []xheads;
	delete []xboxsz;

#if FGT_MANAGE_MEMORY
	delete [] pSources;
	delete [] pWeights;
#endif
}


//////////////////////////////////////////////////////////////
// n choose k
//
// PURPOSE                                                    
// -------   
// Compute the combinatorial number $C_k^n$.
//////////////////////////////////////////////////////////////
template <class TScalar>
int FastGauss<TScalar>::nchoosek(int n, int k)
{
	int n_k = n - k;

	if (k < n_k)
	{
		k = n_k;
		n_k = n - k;
	}

	int  nchsk = 1; 
	for ( int i = 1; i <= n_k; i++)
	{
		nchsk *= (++k);
		nchsk /= i;
	}

	return nchsk;
}



////////////////////////////////////////////////////////////////////
// Taylor Expansion
//
// PURPOSE                                                    
// -------   
// Compute the Taylor expansions from the sources.
//
// PARAMETERS                                                      
// ----------
// d		--> Dimension of problem
// p		--> Order of truncation
//
// OUTPUTS
// -------
// xc		<--	center of the sources
// A_k		<-- Coefficients of the Taylor expansion with wts u
// B_k		<-- Coefficients of the Taylor expansion without wts
////////////////////////////////////////////////////////////////////
template <class TScalar>
void FastGauss<TScalar>::TaylorExpansion(int d, int p)
{
	// compute C_k(p, p)
	TScalar *C_k = new TScalar[pd];
	Compute_C_k(d, p, C_k);

	// compute A_k(p, p) and B_k(p, p)	
	Compute_A_k(xc, C_k);

	delete[] C_k;

	return;
}



template <class TScalar>
void FastGauss<TScalar>
::Compute_Centers(
		TScalar *xc, TScalar *x,
		int d, int n, int *ind, int K, int *bxsz)
{

	for (int i1 = 0; i1 < d*K; i1++)
		xc[i1] = 0.0;


	for (int i2 = 0, nd = 0; i2 < n; i2++,nd+=d)
	{
		int ibase = ind[i2]*d;
		for (int j = 0; j < d; j++)
			xc[ibase+j] += x[nd+j];
	}

	for (int i3 = 0, ibase = 0; i3 < K; i3++,ibase+=d)
		for (int j = 0; j < d; j++)
			xc[ibase+j] /= bxsz[i3];

}

////////////////////////////////////////////////////////////////////
// Compute Taylor Expansion Coefficients
//
// PURPOSE                                                    
// -------   
// Compute the coefficients of the Taylor expansions.
//
// PARAMETERS                                                      
// ----------
// d		--> Dimension of problem
// p		--> Order of truncation
//
// OUTPUTS
// -------
// C_k		<-- Coefficients of the Taylor expansion
////////////////////////////////////////////////////////////////////
template <class TScalar>
void FastGauss<TScalar>::Compute_C_k(int d, int p, TScalar *C_k)
{

	int *heads = new int[d+1];
	int *cinds = new int[pd];

	for (int i = 0; i < d; i++)
		heads[i] = 0;
	heads[d] = INT_MAX;

	cinds[0] = 0;
	C_k[0] = 1.0;
	for (int k=1, t=1, tail=1; k < p; k++, tail=t)
	{
		for (int i = 0; i < d; i++)
		{
			int head = heads[i];
			heads[i] = t;
			for ( int j = head; j < tail; j++, t++)
			{
				cinds[t] = (j < heads[i+1])? cinds[j] + 1 : 1;
				C_k[t] = 2.0 * C_k[j];
				C_k[t] /= (TScalar) cinds[t];
			}
		}
	}

	delete []cinds;
	delete []heads;
	return;
}


////////////////////////////////////////////////////////////////////
// Compute Taylor Expansion
//
// PURPOSE                                                    
// -------   
// Compute the Taylor expansions from the sources.
//
// PARAMETERS                                                      
// ----------
// xc		--> Center of the sources
// C_k		--> Coefficents of the Taylor expansion
//
// OUTPUTS
// -------
// A_k		<-- Coefficients of the Taylor expansion with wts u
////////////////////////////////////////////////////////////////////
template <class TScalar>
void FastGauss<TScalar>::Compute_A_k(TScalar *xc, TScalar *C_k)
{
	TScalar *dx = new TScalar[Dim];
	TScalar *prods = new TScalar[pd];

	int *heads = new int[Dim];

	for (int i = 0; i < pd*K_N*DWeights; i++)
		A_k[i] = 0.0;


	for (int n=0; n < NSources; n++)
	{
		int nbase = n*Dim;
		int ix2c = indx[n];
		int ix2cbase = ix2c*Dim;
		TScalar sum = 0.0;
		for (int i = 0; i < Dim; i++)
		{
			dx[i] = (pSources[nbase+i] - xc[ix2cbase+i])/bandwid;
			sum -= dx[i] * dx[i];
			heads[i] = 0;
		}

		prods[0] = exp(sum);
		for (int k=1, t=1, tail=1; k < pterms; k++, tail=t)
		{
			for (int i = 0; i < Dim; i++)
			{
				int head = heads[i];
				heads[i] = t;
				for ( int j = head; j < tail; j++, t++)
					prods[t] = dx[i] * prods[j];
			} // for i
		} // for k
		for (int dw = 0; dw < DWeights; dw++)
		{
			for (int i = 0; i < pd; i++)
				A_k[(ix2c*DWeights+dw)*pd+i] += pWeights[n*DWeights+dw]*prods[i];
		}
	} // for n
	for (int dw = 0; dw < DWeights; dw++)
		for (int k = 0; k < K_N; k++)
		{
			for (int i=0; i < pd; i++)
				A_k[(k*DWeights+dw)*pd+i] *= C_k[i];
		}

	delete []dx;
	delete []heads;
	delete []prods;

	return;
}




////////////////////////////////////////////////////////////////////
// Compute the summation
//
// PURPOSE                                                    
// -------   
// Compute the summatioin of the Taylor expansions from the sources.
//
// PARAMETERS                                                      
// ----------
// y		--> Target points
// mt		--> Column dimension of the target points
//
// OUTPUTS
// -------
// vi		<-- Results (DWeights-by-mt vector)
////////////////////////////////////////////////////////////////////
template <class TScalar>
void FastGauss<TScalar>
::Compute_v_i(const TScalar *y, const int mt, TScalar *v_i)
 {
	TScalar *dy = new TScalar[Dim];
	TScalar *prods = new TScalar[pd];
	int *heads = new int[Dim];

	for (int m=0; m < mt; m++)
	{	
		for (int dw = 0; dw < DWeights; dw++)
			v_i[m*DWeights+dw] = 0.0;

		int mbase = m*Dim;

		for (int kn=0; kn < K_N; kn++)
		{
			int xbase = kn*Dim;
			TScalar sum2 = 0.0;
			for (int i = 0; i < Dim; i++)
			{
				dy[i] = (y[mbase+i] - xc[xbase+i])/bandwid;
				sum2 += dy[i] * dy[i];
			}

			if (sum2 > ratio_far) continue; //skip to next kn

			for (int i1 = 0; i1 < Dim; i1++)
				heads[i1] = 0;

			prods[0] = exp(-sum2);		
			for (int k=1, t=1, tail=1; k < pterms; k++, tail=t)
			{
				for (int i = 0; i < Dim; i++)
				{
					int head = heads[i];
					heads[i] = t;
					for ( int j = head; j < tail; j++, t++)
						prods[t] = dy[i] * prods[j];
				} // for i
			}// for k
			for (int dw = 0; dw < DWeights; dw++)
			{
				for (int i = 0; i < pd; i++)
					v_i[m*DWeights+dw] += A_k[(kn*DWeights+dw)*pd+i]*prods[i];
			}// for dw
		}// for kn
	}//for m

	delete []dy;
	delete []prods;
	delete []heads;

	return;
 }

// kcenter.cpp -  K-center algorithm, also minmax algorithm
//
// Using Gonzalez's algorithm to get a 2-approximation of the 
// optimal solution.
//
// This is a cpp file for fast Gauss transform algorithm.
// Author: Changjiang Yang.
//

// $Revision: 1.60 $
// $Date: Mon Jul  8 18:11:31 EDT 2002 $

// Add another overload function which generates the backward mapping
// from centers to points. Feb. 10, 2003

// sdist is the square of the distance of two vectors(TScalar)
template <class TScalar>
inline TScalar
FastGauss<TScalar>
::sdist(const int d, const TScalar *x, const TScalar *y)
 {
	TScalar t, s = 0.0;
	for (int i = d; i != 0; i--)
	{
		t = *x++ - *y++;
		s += t * t;
	}

	return s;
 }


// ddist is the square of the distance of two vectors(TScalar)
template <class TScalar>
inline TScalar
FastGauss<TScalar>::ddist(const int d, const TScalar *x, const TScalar *y)
{
	TScalar t, s = 0.0;
	for (int i = d; i != 0; i--)
	{
		t = *x++ - *y++;
		s += t * t;
	}
	return s;
}



//copy vector to another vector
template <class TScalar>
void FastGauss<TScalar>::dcopy(const int d, const TScalar *x, TScalar *y)
{
	for (int i = d; i != 0; i--, x++, y++)
		*y = *x;
	return;
}

//find the largest element from a vector, in case of no intel mkl library.
template <class TScalar>
int FastGauss<TScalar>::idmax(int n, TScalar *x)
{
	int k = 0;
	TScalar t = *x++;
	for (int i = 1; i < n; i++, x++)
		if( t < *x )
		{
			t = *x;
			k = i;
		}
	return k;

}


//K is the number of centers.
template <class TScalar>
TScalar
FastGauss<TScalar>
::dkcenter(const int dim, const int n, const int K,
		const TScalar *x, TScalar *cn, int *cind)
		//dimension,number of nodes, number of centers,
		//nodes(n), centers(K), index to center(n, range:K)
		{

	//const int incu = 1; //unit increment
	TScalar *dist_C = new TScalar[n]; //distances to the center

	//PP
	// randomly pick one node as the first center.
	//int ind = rand() % n;
	int ind = 0;

	// copy the ind-th node to the first center.
	dcopy(dim, x + ind*dim, cn);

	// compute the distances from each node to the first center.
	for (int i = 0; i < n; i++)
	{
		dist_C[i] = (i==ind)? 0.0:ddist(dim, cn, x+i*dim);
		cind[i] = ind;
	}

	for(int k = 1; k < K; k++)
	{	 
		//find the maximum of vector dist_C, i.e., find the node
		//that is farthest away from C
		ind = idmax(n,dist_C);

		cn += dim; //current center.
		//copy the ind-th node to the current center.
		dcopy(dim, x + ind*dim, cn);
		//update the distances from each point to the current center.
		for (int j = 0; j < n; j++)
		{
			TScalar d = (j==ind)? 0.0:ddist(dim, cn, x+j*dim);
			if (d < dist_C[j])
			{
				dist_C[j] = d;
				cind[j] = ind;
			}
		}
	}

	//find the radius of the k-center algorithm.
	ind = idmax(n,dist_C);

	TScalar radius = dist_C[ind];

	delete []dist_C;

	return sqrt(radius);
		}



// Overloaded function without copying the centers, instead 
// using indices to the centers. Include the numbers of nodes
// belonged in each center.
template <class TScalar>
TScalar
FastGauss<TScalar>
::dkcenter(const int dim, const int n, const int K,
		const TScalar *x, int *cn, int *cind, int *cnn)
		//dim--dimension, n--number of nodes, K--number of centers,
		//x--nodes(n), cnt--centers(K), cind--index to center(n)
		//cnn--numbers of nodes(K)
		{
	//const int incu = 1; //unit increment
	TScalar *dist_C = new TScalar[n]; //distances to the center

	//PP
	// randomly pick one node as the first center.
	//int ind = rand() % n;
	int ind = 0;

	// add the ind-th node to the first center.
	*cn++ = ind;

	// compute the distances from each node to the first center.
	const TScalar *x_ind, *x_j;
	x_ind = x + ind*dim;
	x_j = x;
	for (int j = 0; j < n; x_j += dim, j++)
	{
		dist_C[j] = (j==ind)? 0.0:ddist(dim, x_j, x_ind);
		cind[j] = 0;
	}

	for(int i = 1; i < K; i++)
	{	 
		//find the maximum of vector dist_C, i.e., find the node
		//that is farthest away from C
		ind = idmax(n,dist_C);
		*cn++ = ind; //add the ind-th node to the current center.

		//update the distances from each point to the current center.
		x_ind = x + ind*dim;
		x_j = x;
		for (int j = 0; j < n; x_j += dim, j++)
		{
			TScalar d = (j==ind)? 0.0:ddist(dim, x_j, x_ind);
			if (d < dist_C[j])
			{
				dist_C[j] = d;
				cind[j] = i;
			}
		}
	}

	//find the radius of the k-center algorithm.
	ind = idmax(n,dist_C);

	TScalar radius = dist_C[ind];

	delete []dist_C;

	for (int i1 = 0; i1 < K; i1++)
		cnn[i1] = 0;

	for (int i2 = 0; i2 < n; i2++)
		cnn[cind[i2]]++;

	return sqrt(radius);

		}


// Overloaded function without copying the centers, instead 
// using indices to the centers.
template <class TScalar>
TScalar
FastGauss<TScalar>
::dkcenter(const int dim, const int n, const int K,
		const TScalar *x, int *cn, int *cind)
		//dim--dimension, n--number of nodes, K--number of centers,
		//x--nodes(n), cnt--centers(K), cind--index to center(n)
		{

	//const int incu = 1; //unit increment
	TScalar *dist_C = new TScalar[n]; //distances to the center

	//PP
	// randomly pick one node as the first center.
	//int ind = rand() % n;
	int ind = 0;

	// add the ind-th node to the first center.
	*cn++ = ind;

	// compute the distances from each node to the first center.
	const TScalar *x_ind, *x_j;
	x_ind = x + ind*dim;
	x_j = x;
	for (int j = 0; j < n; x_j += dim, j++)
	{
		dist_C[j] = (j==ind)? 0.0:ddist(dim, x_j, x_ind);
		cind[j] = 0;
	}

	for(int i = 1; i < K; i++)
	{	 
		//find the maximum of vector dist_C, i.e., find the node
		//that is farthest away from C
		ind = idmax(n,dist_C);
		*cn++ = ind; //add the ind-th node to the current center.

		//update the distances from each point to the current center.
		x_ind = x + ind*dim;
		x_j = x;
		for (int j = 0; j < n; x_j += dim, j++)
		{
			TScalar d = (j==ind)? 0.0:ddist(dim, x_j, x_ind);
			if (d < dist_C[j])
			{
				dist_C[j] = d;
				cind[j] = i;
			}
		}
	}

	//find the radius of the k-center algorithm.
	ind = idmax(n,dist_C);

	TScalar radius = dist_C[ind];

	delete []dist_C;

	return sqrt(radius);
		}

/*
// Overloaded function without copying the centers, instead 
// using indices to the centers, and get the hist information
// and the indices from centers to nodes.
template <class TScalar>
TScalar
FastGauss<TScalar>
::dkcenter(const int dim, const int n, const int K,
				const TScalar *x, int *cn, int *cind, 
				int *heads, int *boxsz, int *nind)
//dim--dimension, n--number of nodes, K--number of centers, 
//x--nodes(n), cnt--centers(K), cind--index to center(n)
//heads--(K), boxsz--(K), nind--index to point(n)
{
	TScalar r = dkcenter(dim,n,K,x,cn,cind);

	for (int i = 0; i < K; i++)
		boxsz[i] = 0;

	for (int i1 = 0; i1 < n; i1++)
		boxsz[cind[i1]]++;


	int ct = 0;
	for (int i2 = 0; i2 < K; i2++)
	{
		ct += boxsz[i2];
		heads[i2] = ct;
	}

	for (int i3 = n-1; i3 >= 0; i3--)
		nind[--heads[cind[i3]]] = i;

	return r;
}
 */



#endif /* _FastGauss_txx */
