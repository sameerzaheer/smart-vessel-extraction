/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    itkVesselExtractionTools.h
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVesselTools_h
#define __itkVesselTools_h

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
#include <time.h>
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

#include "VesselInfoNew.h"
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
	/** \class VesselExtractionTools
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

	class ITK_EXPORT VesselExtractionTools : public PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet>
	{
	public:

		/** Standard class typedefs. */
		typedef VesselExtractionTools                              Self;
		typedef PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet >  Superclass;

		typedef SmartPointer<Self>         Pointer;
		typedef SmartPointer<const Self>   ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(VesselExtractionTools, Object);

		/** Types transferred from the base class */
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

		typedef TDistanceMap                                    DistanceMapType;
		typedef typename DistanceMapType::ConstPointer          DistanceMapPointer;


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
		//static int Test( int argc, char *argv[] );
		int Test( int argc, char *argv[] );



	protected:
		VesselExtractionTools();
		virtual ~VesselExtractionTools() {};

		/** PrintSelf funtion */
		void PrintSelf(std::ostream& os, Indent indent) const;

	private:
		VesselExtractionTools(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented

		DistanceMapPointer m_DistanceMap;
		bool               m_ComputeSquaredDistance;


		/**** VesselTree ****/

		// // Global variables // //
	public:
		//Image declarations
		typedef float				PixelType;
		const static unsigned int			Dimension = 3;
		//typedef Image< PixelType, Dimension >  ImageType;
		typedef OrientedImage< PixelType, Dimension >  ImageType;
		typedef typename ImageType::IndexType				IndexType; 
		typedef typename ImageType::PointType				PhysicalPointType;
		typedef typename ImageType::SizeType				SizeType;		
		typedef Vector<float,3U>				VectorType;
		typename ImageType::Pointer image;
		VectorType spacing;
		PhysicalPointType origin;

		//parameters
		float scale;// =1; //defined as radius (not diameter) of vessel
		float beta;//= 0.4;
		float sigma;// =2.5;
		float sFactor;//=0.2;
		float sFactor_;
		float intensityThreshold;

		//other parameters
		int max_iterations;//=140;
		float distance_tolerance;
		float distances[20];
		PhysicalPointType prevPoint;
		VectorType prevVector;
		VectorType prevVector2;
		int trackable; //not to rely on directions till its 0
		float epsilon;// =0.8;
		float gradTolerance;//=0.3;
		int windowSize;// =5;

		//seed model parameters

		int maxSeeds;//=100;
		float angThresh;//= 0.75;
		float distThresh;//= 3.0;
		float axialMean;//=4.0;	//in terms of spacing
		float axialSpread;//=2.5; //in terms of spacing
		float radialSpread;//=1.5;//in terms of spacing
		float radialIncrement;//=0.02; //the circle will move out
		//every iteration
		float axialIncrement;//=0.02; //

		//Seed+Direction
		//typedef typename itkVesselInfo::Pointer			InfoArrayType;

		//typedef TFixedPointSet PointSetType;

	  //typedef itk::VesselInfo<PointSetType,PointSetType> VesselInfoType;
		//typedef typename itkVesselInfo::Pointer			InfoArrayType;
		//typedef VesselInfoType InfoArrayType;

		typedef VesselInfo InfoArrayType;

		std::vector<InfoArrayType> infoVector;
		//PhysicalPointType seedPoint;
		VectorType v_direction; //normalized direction

		//Max
		PhysicalPointType physicalPoint; // two d max

		//Hessian declarations
		typedef HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;
		typedef typename HessianFilterType::OutputImageType	HessianImageType;
		typedef typename HessianImageType::Pointer				HessianImagePointer;
		typename HessianFilterType::Pointer	HessianFilter;// = HessianFilterType::New();

		//Line interator declarations
		typedef LineConstIterator<ImageType>			ImageLineItType; 
		typedef LineConstIterator<HessianImageType>	HessianImageLineItType; 

		//Eigen matrix declarations
		typedef vnl_matrix< float >					MType;
		typedef Matrix< float, 3, 3 >			EigenVMType;
		typedef FixedArray< float, 3 >			EigenVAType;
		typedef SymmetricEigenAnalysis< MType, EigenVAType, EigenVMType > EigenAType;

		//Resampling declarations
		//Filter
		typedef ResampleImageFilter<ImageType,ImageType> FilterType;
		typename FilterType::Pointer filter;// = FilterType::New();
		//Spacing
		VectorType reSpacing;
		VectorType reSpacing_;
		//Interpolator
		typedef LinearInterpolateImageFunction< ImageType, double >  InterpolatorType1;  
		//typedef itk::BSplineInterpolateImageFunction< ImageType, double >  InterpolatorType1;  
		typename InterpolatorType1::Pointer interpolator1;// = InterpolatorType1::New();
		//direction
		typename ImageType::DirectionType direction;
		//Output
		typename ImageType::Pointer resampledImage;

		//Data declarations
		vtkCellArray* checkLines;
		vtkPoints* checkPts;
		vtkIdType nLines;
		float directionLength; //in mm

		//To be Objectified
		int checkInd;
		int checkLin;		
		PhysicalPointType newStart;
		PhysicalPointType prevCheck; //previous checkPoint
		PhysicalPointType prevTracked; //previous point tracked
		int prevTrackedI; //index of the previously tracked point
		float distRemainder;
		float angleThresh;//=0.58;
		float distanceThresh;// =1.0;

		PhysicalPointType nextSeed;
		PhysicalPointType nextDirection;
		VectorType nextVector;
		bool firstTime;
		bool termPremature;

		//Random generator
		itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randomGenerator; 

		//Gradient filter (for radius estimation)
		typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType >  gradientType;
		typename ImageType::Pointer gradImage;

		//interpolator for radius can be linear
		typedef itk::LinearInterpolateImageFunction< ImageType, double >  InterpolatorType2;  
		typename InterpolatorType2::Pointer interpolator2;

		//WRITE
		bool write, write1, write2;

		//TIME
		float time1, time2, time3;

		//FAST SWITCH
		bool fastSwitch;

		// // Tools // //
	public:

		////static functions
		float min(float a, float b);
		float max(float a, float b);
		int smallest_index(VectorType A);
		VectorType compute_vector(PhysicalPointType pointA, PhysicalPointType pointB);
		PhysicalPointType compute_point(PhysicalPointType pointA, VectorType vectorB);
		PhysicalPointType convert_point(PhysicalPointType pointIn);
		void cross_product(VectorType A, VectorType B, VectorType* C);
		void compute_planeV(VectorType C, VectorType* A, VectorType* B);
		float compute_distance(PhysicalPointType curPoint, PhysicalPointType refPoint);
		float compute_angle(PhysicalPointType A, PhysicalPointType B, PhysicalPointType C);
		float compute_angle(PhysicalPointType A, PhysicalPointType B, PhysicalPointType C, PhysicalPointType D);

		////Image based tools
		MType calculate_Gradient(IndexType point);
		int OneD_Maximum(PhysicalPointType seedPoint, PhysicalPointType *oneDMax);
		
		void find_Hessian();
		MType compute_EigenMatrix(HessianImagePointer HessianImage, PhysicalPointType physicalPoint);
		MType make_EigenMatrix(VectorType vA, VectorType vB, VectorType vC);
		
		PhysicalPointType TwoD_Maximum(PhysicalPointType physicalIndex,  MType E_matrix);
		PhysicalPointType TwoD_Maximum_fast(PhysicalPointType physicalIndex,  MType E_matrix);
		
		PhysicalPointType move_Point(PhysicalPointType originalPhysicalPoint, MType E_matrix);
		MType flip_sign(MType E_mtx);
		MType change_sign(MType E_mtx_prev, MType E_mtx);
		MType correct_orientation(MType E_mtx, VectorType direction_);
		MType correct_orientation(PhysicalPointType curPoint, PhysicalPointType refPoint, MType E_mtx);
		//seeding
		float refine_seed(PhysicalPointType unrefinedSeed, PhysicalPointType* refinedSeed, VectorType* seedDirection, PhysicalPointType* testResult,MType E_mtx1);
		int optimize_seed(PhysicalPointType seedPoint, VectorType v_dir, PhysicalPointType expectedSeed, PhysicalPointType* optimalPoint);
		int optimize_seed(PhysicalPointType seedPoint, VectorType v_dir, PhysicalPointType* optimalPoint);
		
		void compute_v_dir(int n);

		void initialize_LOOP();
		
		int get_point(vtkSmartPointer<vtkIdList> pLine, int pInd, PhysicalPointType* returnPoint);
		int get_point(int pLine, int pInd, PhysicalPointType* returnPoint);		
		//checking orientation
		void start_checks(int n);
		PhysicalPointType new_checkPoint(float dist);
		int check_orientation(PhysicalPointType currTracked, float dist);
		int check_check(PhysicalPointType currTracked, int currI);		
		int check_line_orientiation(int n, vtkPoints *inPoints, vtkPoints *outPoints, bool termP);
		int check_line_orientiation(int n, vtkPoints *inPoints, vtkPoints *outPoints);
		void print_Emtx(MType E_mtx);
		//termination
		float track_distances(PhysicalPointType curPoint, PhysicalPointType refPoint);
		void track_directions(PhysicalPointType curPoint);
		bool test_NegEvalue(MType E_matrix);
		float test_GradientZero(IndexType point, MType E_matrix);
		float RidgeSymmetry(MType E_matrix);
		bool termination(PhysicalPointType physicalPoint, MType E_matrix, int n);

		float compute_radius(PhysicalPointType physicalIndex, double prev_rad, double est_rad, MType E_mtx);
		void LOOPER();

		//VTK functions (static)
		void reorder_points(vtkPoints *centerPoints, vtkPoints *orderedPoints, int firstpoints);
		void write_points(vtkPoints *inPoints, vtkPoints *allPoints, int finalInPoint);
		void write_points(vtkPoints *inPoints, vtkPoints *allPoints);
		void trim_equal(vtkPoints *pointsA, vtkPoints *pointsAPermanent, int numPoints);
		void smooth_points(vtkPoints *centerPoints, vtkPoints *fixedPoints, int finalPoint);		
		void smooth_points(vtkPoints *centerPoints, vtkPoints *fixedPoints);
		void bad_points(vtkPoints *centerPoints);
		void write_pointsVTK(const char *filename, vtkPoints *points, vtkCellArray *centerLine, vtkFloatArray* radiusArray);

	};



} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVesselExtractionTools.hxx"
#endif

#endif
