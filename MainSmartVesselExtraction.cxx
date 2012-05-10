/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: MainSmartVesselExtraction.cxx,v $
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
#include<iostream>

int TestVesselExtraction( int argc, char *argv[] );

int main( int argc, char *argv[] )
{

	int errorCode = 0;

	//if(argc<4){
	//	std::cout << "Usage: ComputeACDVTKFile.exe <Extracted Vessels.vtk> <Ground Truth .vtk> <Selected CurveIds .txt>" << std::endl;
	//	return 0;
	//}


	TestVesselExtraction(argc, argv);
	//std::cout << "ACD: " << ACD << std::endl;

	std::cout << "... main() " << std::endl;

	return errorCode;
} // end main()

