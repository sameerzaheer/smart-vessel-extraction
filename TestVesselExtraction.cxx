/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: TestVesselExtraction.cxx,v $
 *  Language:  C++
 *  Date:      $Date: 2012-03-31 19:58:58 $
 *  Version:   $Revision: 1.00 $
 *
 *  Copyright (c) CIGITI, HSC. All rights reserved.
 *  See src/CopyrightCIGITI.txt or http://legal.php for details.
 *  
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE. See the above copyright notices for more information.
 *
 *=========================================================================*/



//#include "itkSmartVesselExtraction.h"
#include "itkVesselExtraction.h"


#include "itkPointSet.h"
#include "itkPoint.h"
#include "itkDefaultStaticMeshTraits.h"

//int TestVesselExtraction( int argc, char *argv[] );


//
//****************************************************************************** 
int TestVesselExtraction( int argc, char *argv[] )
{
	std::cout << "argc: " << argc << std::endl;
	std::cout << "argv: " << argv[0] << std::endl;

	//
	//=============================================================================
	const unsigned int NDimensions = 3;
	typedef float TScalarType;
	typedef itk::Point<TScalarType,NDimensions> PointType;
	typedef itk::DefaultStaticMeshTraits<TScalarType, NDimensions, NDimensions, TScalarType, TScalarType> PointSetTraitsType;
	typedef itk::PointSet<PointType, NDimensions, PointSetTraitsType> PointSetType;
  //typedef PointSetType::Pointer                PointSetPointer;
  //typedef PointSetType::PointsContainer        PointsContainer;
	//typedef itk::PointSet<double>  PointSet3DCurveType;

  typedef itk::VesselExtraction<PointSetType,PointSetType> VesselExtractionType;

  typedef VesselExtractionType::Pointer VesselExtractionPointer;
  VesselExtractionPointer Extraction = VesselExtractionType::New();
	//Extraction->Print(std::cout);

	//smartExtraction->Test(argc, argv);
	//int VesselExtraction(int argc, char * argv[]);
	Extraction->Extraction_MAIN(argc, argv);

	//=============================================================================

	//=============================================================================

	return 0;
} // TestVesselExtraction()
//******************************************************************************
