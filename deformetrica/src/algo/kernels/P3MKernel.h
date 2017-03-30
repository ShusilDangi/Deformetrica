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

#ifndef _P3MKernel_h
#define _P3MKernel_h

#include "ExactKernel.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkSimpleFastMutexLock.h"

/**
 *	\brief 		P3M kernels.
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The P3MKernel class inherited from AbstractKernel implements operations with kernels using approximations based on FFTs and projections/interpolation on regular lattices.
 */
template <class TScalar, unsigned int PointDim>
class P3MKernel: public ExactKernel<TScalar, PointDim>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Exact kernel type.
	typedef ExactKernel<TScalar, PointDim> Superclass;

	/// Vector type
	typedef typename Superclass::VNLVectorType VNLVectorType;
	/// Matrix type.
	typedef typename Superclass::VNLMatrixType VNLMatrixType;

	/// Image type (itk).
	typedef itk::Image<TScalar, PointDim> ImageType;
	/// Image pointer type (itk).
	typedef typename ImageType::Pointer ImagePointer;
	/// Image index type (itk).
	typedef typename ImageType::IndexType ImageIndexType;
	/// Image point type (itk).
	typedef typename ImageType::PointType ImagePointType;
	/// Image region type (itk).
	typedef typename ImageType::RegionType ImageRegionType;
	/// Image size type (itk).
	typedef typename ImageType::SizeType ImageSizeType;
	/// Image spacing type (itk).
	typedef typename ImageType::SpacingType ImageSpacingType;

	/// Complex image type (itk).
	typedef itk::Image< std::complex<TScalar>, PointDim> ComplexImageType;
	/// Complex image pointer type (itk).
	typedef typename ComplexImageType::Pointer ComplexImagePointer;

	/// Image iterator type (itk).
	typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;

	/// List of images type (itk).
	typedef std::vector<ImagePointer> ImageListType;
	/// Image matrix type (itk).
	typedef std::vector<ImageListType> ImageMatrixType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	P3MKernel();
	/// Copy constructor.
	P3MKernel(const P3MKernel& o);
	/// See AbstractKernel::AbstractKernel(const VNLMatrixType& X, double h) for details.
	P3MKernel(const VNLMatrixType& X, double h);
	/// See AbstractKernel::AbstractKernel(const VNLMatrixType& X, const VNLMatrixType& W, double h) for details.
	P3MKernel(const VNLMatrixType& X, const VNLMatrixType& W, double h);

	virtual P3MKernel* Clone() const;

	virtual ~P3MKernel();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Set the origin of the grid to \e d.
	void SetGridOrigin(const ImagePointType& p) { m_GridOrigin = p; }

	/// Set the grid spacing to \e d.
	void SetGridSpacing(const ImageSpacingType& s) { m_GridSpacing = s; }

	/// Set the size of the grid to \e size.
	void SetGridSize(const ImageSizeType& size) { m_GridSize = size; }

	/// Set the size of the grid padding.
	void SetGridPadding(long n) { m_GridPadding = n; }

	/// Set the min and max values of the data domain in each dimension.
	void SetDataDomain(VNLMatrixType DD) { m_DataDomain = DD; }

	/// Set the ratio between the grid step and the kernel width. 0.2 gives a relative approximation error of about 5%, use 0.3 for increased speed
	void SetWorkingSpacingRatio(const TScalar d) { m_WorkingSpacingRatio = d; }

	/// Sets the padding factor to \e d. It will enlarge the grid by m_PaddingFactor x m_KernelWidth to avoid side effects
	/// (FFTs have circular boundary conditions). It is also used to define a bounding box.
	void SetPaddingFactor(const TScalar d) { m_PaddingFactor = d; }



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// Automatically called before convolution calls
	/// Splat vectors at grid points, compute its FFT, and call BuildFFTKernel() and BuildFFTGradientKernels()
	virtual void UpdateGrids();
	/// Essentially call BuildFFTHessianKernels()
	virtual void UpdateHessianGrids();

	void SetNearThresholdScale(double s) { m_NearThresholdScale = s; }

	virtual VNLMatrixType Convolve(const VNLMatrixType& X);

	virtual std::vector<VNLMatrixType> ConvolveGradient(const VNLMatrixType& X);

	virtual VNLVectorType ConvolveGradient(const VNLMatrixType& X, unsigned int k, unsigned int dp);
	virtual VNLMatrixType ConvolveGradient(const VNLMatrixType& X, unsigned int dim);

	virtual std::vector< std::vector<VNLMatrixType> > ConvolveHessian(const VNLMatrixType & X);

	virtual VNLVectorType ConvolveHessian(const VNLMatrixType & X, unsigned int k, unsigned int dp, unsigned int dq);
	virtual VNLMatrixType ConvolveHessian(const VNLMatrixType& X, unsigned int row, unsigned int col);

	/// Compute the FFT of the kernel sampled at grid points.
	ComplexImagePointer BuildFFTKernel();
	/// Compute the FFT of each dimension of the gradient of the kernel.
	std::vector<ComplexImagePointer> BuildFFTGradientKernels();
	/// Compute the FFT of each entry of the Hessian matrix of the kernel.
	std::vector<ComplexImagePointer> BuildFFTHessianKernels();

	// Methods to override the FFT kernel used
	void UseFFTKernel(ComplexImagePointer kernImg) { m_FFTKernel = kernImg; }
	void UseFFTGradientKernels(std::vector<ComplexImagePointer>& klist) { m_FFTGradientKernels = klist; }
	void UseFFTHessianKernels(std::vector<ComplexImagePointer>& klist) { m_FFTHessianKernels = klist; }



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	void Init();

	void ClearGrids();
	/// Set grids for splatting
	void DetermineGrids();

	VNLVectorType Interpolate(const VNLVectorType& x, const std::vector<ImagePointer>& imgs);

	void SplatToGrid(ImageType* mesh, const VNLMatrixType& X, const VNLVectorType& values);

	ImagePointer ApplyKernelFFT(ComplexImageType* kernelImg, ImageType* img);
	ImagePointer ApplyKernelFFT(ComplexImageType* kernelImg, ComplexImageType* img);


	/// Get weights for linear interpolation of data at grid points
	void inline _getInterpolationWeightsAndGridPoints(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VNLVectorType& x);
	/// Specialized linear interpolation weights for 2D.
	void inline _getInterpolationWeightsAndGridPoints_2(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VNLVectorType& x);
	/// Specialized linear interpolation weights for 3D.
	void inline _getInterpolationWeightsAndGridPoints_3(
			std::vector<TScalar>& weights,
			std::vector<ImageIndexType>& gridIndices,
			const VNLVectorType& x);

	long inline _getPointID(const VNLVectorType& x);
	long inline _getPointID(const ImageIndexType& x);



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	VNLMatrixType m_DataDomain;

	ImagePointType m_GridOrigin;

	ImageSpacingType m_GridSpacing;

	ImageSizeType m_GridSize;

	long m_GridPadding;

	TScalar m_WorkingSpacingRatio;

	TScalar m_PaddingFactor;


	/// \cond HIDE_FOR_DOXYGEN

	ComplexImagePointer m_FFTKernel;
	std::vector<ComplexImagePointer> m_FFTGradientKernels;
	std::vector<ComplexImagePointer> m_FFTHessianKernels;

	// Splatted values after convolution
	std::vector<ImagePointer> m_MeshList;
	std::vector<ComplexImagePointer> m_MeshListFFT;

	ImageMatrixType m_MeshGradientList;
	std::vector< ImageMatrixType > m_MeshHessianList;

	bool m_HessianUpdated;

	std::vector<unsigned int> m_SourceIDs;

	std::vector< std::vector<unsigned int> > m_IDLookupTable;

	TScalar m_NearThresholdScale;

	static itk::SimpleFastMutexLock m_FFTMutex;

	static itk::SimpleFastMutexLock m_CacheMutex;

	// Cache for BuildFFT...Kernel()
	static std::vector<TScalar> m_CacheFFTKernelWidths;
	static std::vector<ImageSizeType> m_CacheFFTKernelSizes;
	static std::vector<ComplexImagePointer> m_CacheFFTKernelImages;

	static std::vector<TScalar> m_CacheFFTGradientKernelWidths;
	static std::vector<ImageSizeType> m_CacheFFTGradientKernelSizes;
	static std::vector< std::vector<ComplexImagePointer> >  m_CacheFFTGradientKernelImages;

	static std::vector<TScalar> m_CacheFFTHessianKernelWidths;
	static std::vector<ImageSizeType> m_CacheFFTHessianKernelSizes;
	static std::vector< std::vector<ComplexImagePointer> > m_CacheFFTHessianKernelImages;

	/// \endcond

}; /* class P3MKernel */


#ifndef MU_MANUAL_INSTANTIATION
#include "P3MKernel.txx"
#endif


#endif /* _P3MKernel_h */
