#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//#include "VesselTree9.h"

#include "itkPointSet.h"


#include "itkVesselExtractionTools.h"
 
#include <string>


int main( int argc, char * argv[])
{  
  if(argc != 4)
	{
		std:: cout << "Usage: VesselTree9 <image filename> <seed filename> <lines filename>" << std::endl;
		std:: cout << "Example: Build\Debug>VesselTree9 EHPS3R3Crop.mhd EHPS2103BP806Tfm3_2in.txt EHPS2103BP806Tfm3.vtk" << std::endl;
		exit(-1);
	} 	
 
	typedef itk::PointSet<> PointSetType;
	typedef itk::VesselExtractionTools<PointSetType,PointSetType> VesselExtractionToolsType;
	VesselExtractionToolsType::Pointer vesselTools = VesselExtractionToolsType::NEW();
	vesselTools->VesselExtraction(argc,argv);

	//int VesselExtraction(int argc, char * argv[]);


	std::cout << "File written" << std::endl;
  return 0;
}