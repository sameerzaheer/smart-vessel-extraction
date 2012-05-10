
#ifndef VESSELSEEDNew_H
#define VESSELSEEDNew_H
#include "itkImage.h"


/** Note:  
  */
#define SPLINECURVENUMPOINTS 300
#define detCZeroThreshold 1e-10
#define ShortestDistThreshold 10

#define SHORTESTDISTTHRESHOLD_ED 1000

//only for test data
#define SPLINECURVENUMPOINTS_TEST SPLINECURVENUMPOINTS




class VesselInfo
{

public:
	static const unsigned int			Dimension = 3;
	typedef float				PixelType;
	typedef itk::Image< PixelType, Dimension >  ImageType;
	typedef ImageType::PointType				PhysicalPointType;
	typedef itk::Vector<float,3U>				VectorType;
	
	unsigned int id;
	float length;
	PhysicalPointType seed; // the true seed
	PhysicalPointType seedU; //the updated seed
	unsigned int seedLine;
	unsigned int seedLineInd;
	PhysicalPointType direction;
	PhysicalPointType direction2;
	float isLeaf; //1 if leaf, -1 if not
	PhysicalPointType endPoint;

	//Immediately computed values
	VectorType v_direction;
	VectorType v_direction2;
	VectorType v_direction3;
	float radiusS;

	//Values computed later
	PhysicalPointType newSeed;
	PhysicalPointType oldNewSeed;
	VectorType startToEnd;
	unsigned int numPoints; //the number of points identified
	unsigned int numTries; 

};

#endif