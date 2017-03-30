/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _DeformableObject_h
#define _DeformableObject_h

#include "Diffeos.h"

#include "vnl/vnl_matrix.h"

#include <vector>



/**
 *	\brief 		Deformable objects
 *
 *	\copyright		Inria and the University of Utah
 *	\version		Deformetrica 2.0
 *
 *	\details	The DeformableObject class encodes an object embedded in the current 2D or 3D space. This object is designed to deform diffeomorphic deformations. \n \n
 *			See DeformableObject::DeformableObjectType for the list of available object types.
 */
template <class TScalar, unsigned int Dimension>
class DeformableObject
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	///	Possible type of deformable object.
	typedef enum
	{
		null,			/*!< Null value. */
		Landmark,		/*!< Landmark (see Landmark). */
		OrientedPolyLine,	/*!< Current representation of a curve (see OrientedPolyLine). */
		OrientedSurfaceMesh,	/*!< Current representation of a surface (see OrientedSurfaceMesh). */
		NonOrientedPolyLine,	/*!< Varifold representation of a curve (see NonOrientedPolyLine). */
		NonOrientedSurfaceMesh,	/*!< Varifold representation of a surface (see NonOrientedSurfaceMesh). */
		PointCloud		/*!< Point cloud (see PointCloud). */
	} DeformableObjectType;

	/// Vector type.
	typedef vnl_vector<TScalar> VNLVectorType;
	/// Matrix type.
	typedef vnl_matrix<TScalar> VNLMatrixType;
	/// List of matrices type.
	typedef std::vector<VNLMatrixType> VNLMatrixList;

	/// Deformation type.
	typedef Diffeos<TScalar, Dimension> DiffeosType;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	DeformableObject();
	/// Copy constructor.
	DeformableObject(const DeformableObject& other);

	/// Makes a copy of the object.
	virtual DeformableObject* Clone() = 0;

	virtual ~DeformableObject();



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Returns \f$\sigma^2\f$.
	inline TScalar GetDataSigmaSquared() { return m_DataSigmaSquared; }
	/// Sets \f$\sigma^2\f$ to \e d*d.
	inline void SetDataSigma(TScalar d) { m_DataSigmaSquared = d*d; }
	/// Sets \f$\sigma^2\f$ to \e d2.
	inline void SetDataSigmaSquared(TScalar d2) { m_DataSigmaSquared = d2; }

	/// Sets the deformation to \e def.
	inline void SetDiffeos(DiffeosType* def) { m_Def = def; m_Modified = true; }

	/// Returns true if the parameters of the deformation have changed, false otherwise.
	inline bool isModified() { return m_Modified; }
	/// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
	inline void SetModified() { m_Modified = true; }
	/// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
	inline void SetUnModified() { m_Modified = false; }

	/// Returns the type of the deformable object.
	inline DeformableObjectType GetType() { return m_Type; }
	/// Returns true if the deformable object is of landmark kind, false otherwise.
	inline bool IsOfLandmarkKind() { return ( (m_Type == Landmark) || (m_Type == PointCloud) ||
			(m_Type == OrientedPolyLine) || (m_Type == OrientedSurfaceMesh) ||
			(m_Type == NonOrientedPolyLine) || (m_Type == NonOrientedSurfaceMesh)  ); }
	/// Returns true if the deformable object is a mesh, false otherwise.
	inline bool IsMesh() {return ( (m_Type == OrientedPolyLine) || (m_Type == NonOrientedPolyLine) ||
			(m_Type == OrientedSurfaceMesh) || (m_Type == NonOrientedSurfaceMesh) ); }
	/// Sets the type of the deformable object to Landmark.
	inline void SetLandmarkType() { m_Type = Landmark; }
	/// Sets the type of the deformable object to PointCloud.
	inline void SetPointCloudType() { m_Type = PointCloud; }
	/// Sets the type of the deformable object to OrientedPolyLine.
	inline void SetOrientedPolyLineType() { m_Type = OrientedPolyLine; }
	/// Sets the type of the deformable object to OrientedSurfaceMesh.
	inline void SetOrientedSurfaceMeshType() { m_Type = OrientedSurfaceMesh; }
	/// Sets the type of the deformable object to NonOrientedPolyLine.
	inline void SetNonOrientedPolyLineType() { m_Type = NonOrientedPolyLine; }
	/// Sets the type of the deformable object to NonOrientedSurfaceMesh.
	inline void SetNonOrientedSurfaceMeshType() { m_Type = NonOrientedSurfaceMesh;}

	/// Returns true if any point is out of the bounding box, false otherwise.
	inline bool OutOfBox() { return m_OutOfBox; }

	/// Returns the number of points of the object.
	virtual int GetNumberOfPoints() const = 0;



	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Sets the deformation (resp. the trajectory) to \e def (resp. to \e traj).
	/// \warning Use with care. Use this->SetDiffeos(def); this->Update(); as much as possible.
	virtual void SetDiffeoAndPointTrajectory(DiffeosType* def, VNLMatrixList& traj) = 0;

	/// Computes the trajectory of the deformable object.
	virtual void Update() = 0;

	/// Returns the flow i.e. the trajectory of the deformable object under the deformation \e m_Def.
	virtual VNLMatrixList GetFlow() = 0;

	/// Returns the norm of the difference at final time.
	virtual TScalar ComputeMatch(DeformableObject* target) = 0;
	/**
	 *	\brief		Returns the norm between itself and the target.
	 *
	 *	\details	Denoting by \e S the source (i.e. the class itself) and by \e T the target object, this
	 *				method computes the following norm : \f[\|S-T\| ,\f] where the norm depends on the child
	 *				class.
	 *
	 *	\param[in]	target	Target of same type as the deformable object.
	 *	\param[in]	t	Time index.
	 *	\return		The norm of the difference.
	 */
	virtual TScalar ComputeMatch(DeformableObject* target, unsigned int t) = 0;

	/// Returns the gradient of the norm of the difference at final time.
	virtual VNLMatrixType ComputeMatchGradient(DeformableObject* target) = 0;
	/**
	 *	\brief		Returns the gradient of the norm between itself and the target.
	 *
	 *	\details	Denoting by \e S the source (i.e. the class itself) and by \e T the target object, this
	 *				method computes the following gradient : \f[ \nabla\|S-T\| , \f] where the norm
	 *				depends on the child class. The gradient is taken with respect to the position of the points (or vertices) of the source object.
	 *
	 *	\param[in]	target	Target of same type as the deformable object.
	 *	\param[in]	t	Time index.
	 *	\return		The gradient of the RKHS-norm of the difference.
	 */
	virtual VNLMatrixType ComputeMatchGradient(DeformableObject* target, unsigned int t) = 0;

	/// Compute the adjoint equation of the flow equation. It is used to transport the gradient of the dissimilarity metric from endpoint back to t=0
	virtual void TransportAlongGeodesic(VNLMatrixType& Eta0, VNLMatrixList& EtaT, VNLMatrixList& YT) = 0;

	/// Saves in \e filename the deformable object (for Landmark type, it is saved in *.vtk format).
	virtual void WriteDeformedObjectAt(unsigned int t, std::string filename) = 0;

	/// Saves the flow of the deformable object in m_NumberOfTimePoints files
	/// (Name of the file : ${name}__t_${indexTimePoints}${extension}).
	virtual void WriteFlow(std::string name, std::string extension)
	{
		//std::cout << "Writing flow..." << std::endl;
		this->Update();

		for (int t = 0; t < m_Def->GetNumberOfTimePoints(); t++)
		{
			std::ostringstream oss;
			oss << name << "__t_" << t << extension << std::ends;
			this->WriteDeformedObjectAt(t, oss.str());
		}
		//std::cout << "Done writing flow..." << std::endl;
	}



protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Type of the deformable object.
	DeformableObjectType m_Type;

	/// \f$\sigma^2\f$ appearing in the fidelity-to-data term.
	TScalar m_DataSigmaSquared;

	/// Entity representing the deformation.
	DiffeosType* m_Def;

	/// Boolean which avoids re-updating trajectory if no parameter has changed.
	bool m_Modified;

	///	Boolean which prevents computations if any point is outside the bounding box. \n
	/// See AbstractDeformations::m_BoundingBox.
	bool m_OutOfBox;


}; /* class DeformableObject */


#ifndef MU_MANUAL_INSTANTIATION
#include "DeformableObject.txx"
#endif


#endif /* _DeformableObject_h */
