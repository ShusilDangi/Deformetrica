#ifndef _FastGauss_h
#define _FastGauss_h

/**
 *	\brief 		Fast Gauss Transform interface.
 *
 *	\details	The FastGauss class is an adaptation of the Fast Gauss Transform initially
 *				developed by Changjiang Yang. \n
 *				(License : http://www.umiacs.umd.edu/~yangcj/license.html)
 */
template <class TScalar>
class FastGauss
{

public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Copy constructor.
	FastGauss(const FastGauss& o);

	/**
	 *	\brief		Constructor.
	 *
	 *	\details	This constructor initialize the class without the parameter \e tol, allocates memory
	 *				for the fast Gauss algorithm, determines the memory space needed, and computes the
	 *				coefficients of the Taylor expansions.
	 *
	 *	\param[in]	dim			Dimension of problem.
	 *	\param[in]	dw			Row dimension of matrix \e u (default = 1).
	 *	\param[in]	u			Weight matrix with \e du rows and \e ns columns.
	 *	\param[in]	x			Source matrix with \e dim rows and \e ns columns.
	 *	\param[in]	ns			Column dimension of matrix \e x (no. of sources).
	 *	\param[in]	h			Bandwidth.
	 *	\param[in]	p			Order of truncation.
	 *	\param[in]	K			Number of centers.
	 *	\param[in]	e			Ratio of far field.
	 *	\param[in]	ptRadius	Pointer to radius of clustering.
	 */
	FastGauss(int dim, int dw, const TScalar *u, const TScalar *x, int ns,
			TScalar h, int p, int K, TScalar e, TScalar* ptRadius);

	/// Makes a copy of the object.
	FastGauss* Clone() const { return new FastGauss(*this); }

	~FastGauss();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns the sources.
	TScalar* GetSources() { return pSources; }

	/// Returns the weights.
	TScalar* GetWeights() { return pWeights; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * \brief		Computes the summation.
	 *
	 * \details		Computes the summation of the Taylor expansions from the sources.
	 * \param[in]	y		Target points.
	 * \param[in]	mt		Column dimension of the target points.
	 * \param[out]	v_i		Results (DWeights-by-mt vector).
	 */
	void Compute_v_i(const TScalar *y, const int mt, TScalar *v_i);



private:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	int nchoosek(int n, int k);

	void TaylorExpansion(int d, int p);
	void Compute_C_k(int d, int p, TScalar *C_k);
	void Compute_A_k(TScalar *xc, TScalar *C_k);
	void Compute_Centers(TScalar *xc, TScalar *x, int d, int n, int *ind, int K, int *bxsz);


	// K-centers algorithm
	TScalar dkcenter(const int dim, const int n, const int K,
			const TScalar *x, TScalar *cn, int *cind);

	TScalar dkcenter(const int dim, const int n, const int K,
			const TScalar *x, int *cn, int *cind);

	TScalar dkcenter(const int dim, const int n, const int K,
			const TScalar *x, int *cn, int *cind, int *cnn);

	TScalar dkcenter(const int dim, const int n, const int K,
			const TScalar *x, int *cn, int *cind,
			int *heads, int *boxsz, int *nind);

	inline TScalar sdist(const int d, const TScalar *x, const TScalar *y);
	inline TScalar ddist(const int d, const TScalar *x, const TScalar *y);
	void dcopy(const int d, const TScalar *x, TScalar *y);
	int idmax(int n, TScalar *x);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// interface:

	/// Dimension of the points.
	int Dim;

	/// Number of sources.
	int NSources;

	/// Dimension of weights.
	int DWeights;

	/// Pointer to sources, (Dim*Nsources).
	TScalar *pSources;

	/// Pointer to sigmas (Dim*Nsources).
	TScalar *pSigmas;
	/// Pointer to weights (DWeights*NSources).
	TScalar *pWeights;
	/// Pointer to radius of clustering (output).
	TScalar *ptRadius;


	/// Bandwidth.
	TScalar bandwid;
	/// Desired precision.
	TScalar epsilon;
	/// Number of truncation terms.
	int pterms;

	// internal data:
	/// Number of coefficients.
	int pd;

	/// Pointer to coefficients.
	TScalar *A_k;

	int K_N;
	TScalar ratio_far;

	TScalar *xc;
	/// Center of the sources.
	int *xcind;
	int *indx;
	int *xheads;
	int *xboxsz;
	TScalar rx2;



}; /* class FastGauss */


#ifndef MU_MANUAL_INSTANTIATION
#include "FastGauss.txx"
#endif


#endif /* _FastGauss_h */
