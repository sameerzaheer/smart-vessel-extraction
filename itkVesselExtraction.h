/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    itkVesselExtraction.h
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVesselExtraction_h
#define __itkVesselExtraction_h

//SameerW2012
//===============================================================================
#include "itkAffineTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkFixedArray.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkLineConstIterator.h"
#include <math.h>
#include "itkMatrix.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkNormalVariateGenerator.h"
#include "itkResampleImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include <itkSymmetricEigenAnalysis.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//#include "FastNorm.h" //can go faster with this
#include "vnl/vnl_matrix.h"

#include <vtkCellArray.h>
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkDataArray.h"
#include <vtkFloatArray.h>
#include "vtkIdTypeArray.h"
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <algorithm>
#include <string>
#include<iostream>
#include<fstream>

#include "itkVesselExtractionTools.h"
//===============================================================================

#include "itkPointSetToPointSetMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkPointSet.h"

#include "itkBSplineResampleImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkImage.h"
#include "itkSize.h"
#include "itkRandomImageSource.h"
//#include "itkFilterWatcher.h"
#include "vnl/vnl_sample.h"


//For PointSet
//==============================================================================
#include "itkPointSet.h"

#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/vnl_sample.h"
#include "vnl/vnl_diag_matrix.h"

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include <testlib/testlib_test.h>
#include <vnl/algo/vnl_real_eigensystem.h>
//==============================================================================



//VTK
//==============================================================================
#include "vtkMath.h"

#include "vtkCellLocator.h"
//#include "vtkSTLReader.h"
//#include "vtkSTLWriter.h"
#include "vtkPolyDataNormals.h"
//#include "vtkThinPlateSplineTransform.h"
//#include "vtkTransformPolyDataFilter.h"
#include "vtkPointDataToCellData.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"

#include "vtkCellData.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"
//#include "vtkIdType.h"

#include "vtkMetaImageReader.h"
#include "vtkMetaImageWriter.h"
#include "vtkOutputWindow.h"
#include "vtkObjectFactory.h"
#include "vtkDataObject.h"
#include "vtkImageData.h"
#include "vtkImageMathematics.h"

#include <vtkMINCImageReader.h>
#include <vtkMINCImageWriter.h> 
//==============================================================================




//==============================================================================
#include "itkBSplineResampleImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkImage.h"
#include "itkSize.h"
#include "itkRandomImageSource.h"
//#include "itkFeatureDetect3DImageFilterCIGITI.h"
//#include "itkFilterWatcher.h"

#include "vnl/vnl_sample.h"
//==============================================================================


/** Note:  This is the same test used for the itkBSplineResampleImageFunctionTest
*        It is duplicated here because it excercises the itkBSplineDecompositionFilter
*        and demonstrates its use.
*/
#define SPLINECURVENUMPOINTS 300
#define PI 3.141592
#define detCZeroThreshold 1e-10
#define ShortestDistThreshold 10

#define SHORTESTDISTTHRESHOLD_ED 1000

//only for test data
#define SPLINECURVENUMPOINTS_TEST SPLINECURVENUMPOINTS


namespace itk
{
	/** \class VesselExtraction
	* \brief Computes the minimum distance between a moving point-set
	*  and a fixed point-set. A vector of minimum closest point distance is
	*  created for each point in the moving point-set.
	*  No correspondance is needed.
	*  For speed consideration, the point-set with the minimum number of points
	*  should be used as the moving point-set.
	*  If the number of points is high, the possibility of setting a distance map
	*  should improve the speed of the closest point computation.
	*
	*  Reference: "A Method for Registration of 3-D Shapes",
	*             IEEE PAMI, Vol 14, No. 2, February 1992
	*
	* \ingroup RegistrationMetrics
	*/
	template < class TFixedPointSet, class TMovingPointSet, 
	class TDistanceMap = ::itk::Image<unsigned short,::itk::GetPointSetDimension<TMovingPointSet>::PointDimension> >

	//class ITK_EXPORT VesselExtraction : public PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet>
	class ITK_EXPORT VesselExtraction : public VesselExtractionTools< TFixedPointSet, TMovingPointSet, TDistanceMap >
	{
	public:

		/** Standard class typedefs. */
		typedef VesselExtraction                              Self;
		typedef PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet >  Superclass;

		typedef SmartPointer<Self>         Pointer;
		typedef SmartPointer<const Self>   ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(VesselExtraction, Object);

		/** Types transferred from the base class */
		//==========================================================================
		typedef typename Superclass::TransformType              TransformType;
		typedef typename Superclass::TransformPointer           TransformPointer;
		typedef typename Superclass::TransformParametersType    TransformParametersType;
		typedef typename Superclass::TransformJacobianType      TransformJacobianType;

		typedef typename Superclass::MeasureType                MeasureType;
		typedef typename Superclass::DerivativeType             DerivativeType;
		typedef typename Superclass::FixedPointSetType          FixedPointSetType;
		typedef typename Superclass::MovingPointSetType         MovingPointSetType;
		typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
		typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;

		typedef typename Superclass::PointIterator              PointIterator;
		typedef typename Superclass::PointDataIterator          PointDataIterator;
		//==========================================================================

		typedef TDistanceMap                                    DistanceMapType;
		typedef typename DistanceMapType::ConstPointer          DistanceMapPointer;



		//Types from itk::VesselExtractionTools
		//==========================================================================
		//typedef VesselExtractionTools VesselExtractionToolsType;
		//typedef typename VesselExtractionToolsType::MType              MType;
		//==========================================================================

		/** Get the number of values */
		unsigned int GetNumberOfValues() const;

		/** Get the derivatives of the match measure. */
		void GetDerivative( const TransformParametersType & parameters,
			DerivativeType & Derivative ) const;

		/**  Get the value for single valued optimizers. */
		MeasureType GetValue( const TransformParametersType & parameters ) const;

		/**  Get value and derivatives for multiple valued optimizers. */
		void GetValueAndDerivative( const TransformParametersType & parameters,
			MeasureType& Value, DerivativeType& Derivative ) const;

		/** Set/Get the distance map */
		itkSetConstObjectMacro(DistanceMap,DistanceMapType);
		itkGetConstObjectMacro(DistanceMap,DistanceMapType);

		/** Set/Get if the distance should be squared. Default is true for computation speed */
		itkSetMacro(ComputeSquaredDistance,bool);
		itkGetConstMacro(ComputeSquaredDistance,bool);
		itkBooleanMacro(ComputeSquaredDistance);
		//---------------------------------------


		//EdW2012
		//void TransformPointsSomePointsVTK( const std::string filename ) const;
		int Test( int argc, char *argv[] );



	protected:
		VesselExtraction();
		virtual ~VesselExtraction() {};

		/** PrintSelf funtion */
		void PrintSelf(std::ostream& os, Indent indent) const;

	private:
		VesselExtraction(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented

		DistanceMapPointer m_DistanceMap;
		bool               m_ComputeSquaredDistance;


		/**** VesselTree ****/
	public:
		int maxJ;
	
		/*Tools*/
	public:

		int read_file(char * argv[]);
		int read_data(char * argv[]);
		int generate_info();
		void INITIALIZE();
		int tracking_LOOP(int n, vtkPoints *centerPoints, vtkFloatArray* radiusArray);
		void Vessel_tracking(int n, vtkPoints *smoothedPoints, vtkFloatArray* radiusArray);
		void Vessel_tracking_orientation(int n, vtkPoints *trimmedPoints, vtkFloatArray* radiusArray);
		int Extraction_MAIN(int argc, char * argv[]);
	};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVesselExtraction.hxx"
#endif

#endif
