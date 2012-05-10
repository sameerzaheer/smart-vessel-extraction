#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "VesselTree9.h"
 
#include <string>


int main( int argc, char * argv[])
{  
  if(argc != 4)
	{
		std:: cout << "Usage: VesselTree9 <image filename> <seed filename> <lines filename>" << std::endl;
		std:: cout << "Example: Build\Debug>VesselTree9 EHPS3R3Crop.mhd EHPS2103BP806Tfm3_2in.txt EHPS2103BP806Tfm3.vtk" << std::endl;
		exit(-1);
	} 	
  
  /************Declarations***********
  */
  int j;
	int continue_, continue__;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  vtkPoints *centerPoints = vtkPoints::New(); 
  vtkPoints *seedPoints = vtkPoints::New(); 
  vtkPoints *orderedPoints = vtkPoints::New();
  vtkCellArray *centerLine = vtkCellArray::New();
  vtkFloatArray* radiusArray = vtkFloatArray::New();
  
  MType E_mtx;  
  PhysicalPointType pPoint; 
  	  
  PhysicalPointType physicalEstimation;
  PhysicalPointType displaceme;
  int totalPoints, numPoints;
  float pix_=1;
  float distance_;
  float radius_prev;
	float score=0, score2;
  
  bool end_loop=false;

	std::cout << "Initializing... ";
  /************Initializations************/ 
  //Read image
  ReaderType::Pointer reader = ReaderType::New();
  const char * filename = argv[1];
  reader->SetFileName( filename ); 
  image = reader->GetOutput();
  reader->Update();
	std::cout << "read image ";
  //initialize certain parameters
  initialize();
  radiusArray->SetNumberOfComponents(1);
  radiusArray->SetName("RadiusArray");  
  
  //Get seed and target points
  if (readfile(argv)==-1){
	  exit(-1);
	} 
  if (read_data(argv)==-1){
	
	 exit(-1);
  } 
	std::cout << " read input" << std::endl;
  
  //get_seed_target(argv);   
/*
	start_checks(0);
			exit(-1);*/
  
  //Get Hessian
	std::cout << "Processing Hessian ...";
  find_Hessian();  
  std::cout << "done" << std::endl;

  /************Refine seed***********/  
	
	
  /////////////TRACK THE VESSELS/////////////  

	totalPoints=0;
	for (int n=0; n<infoVector.size();n++){ //infoVector.size();n++){ 
		vtkPoints* linePoints = vtkPoints::New(); 
		vtkPoints *smoothedPoints = vtkPoints::New();
		vtkPoints *trimmedPoints = vtkPoints::New();
		j=0;
		continue_=1; //1 here means yes
		std::cout << "-------N " << n << std::endl;
		
		firstTime=true;

		/////////////Track vessel/////////////
		while (continue_==1 && j<5){
			
			vtkPoints* segmentPoints = vtkPoints::New(); 
			vtkPoints *smoothedSPoints = vtkPoints::New();
			vtkPoints *trimmedSPoints = vtkPoints::New();
			numPoints=MAIN_LOOP(n, segmentPoints,  radiusArray); 			
		//	std::cout << "AT: " << j << " #Pts=" << numPoints << std::endl;
			smooth_points(segmentPoints,smoothedSPoints,numPoints);  //smoothen the points
			continue_=check_line_orientiation(n, smoothedSPoints, trimmedSPoints);
			if (termCircles)
				continue_=1;
			write_points(trimmedSPoints, linePoints); //linePoints are all the points
			j++;
		}
		smooth_points(linePoints,smoothedPoints);  //smoothen the points
		
		infoVector[n].newSeed=infoVector[n].oldNewSeed;
		check_line_orientiation(n, smoothedPoints, trimmedPoints);
		/////////////Creating output/////////////

		//add them to the array of all center points
		write_points(trimmedPoints, centerPoints);
		//find total number of points
		numPoints=(int)linePoints->GetNumberOfPoints();
		//std::cout << "Actual Points =      " << (int)smoothedPoints->GetNumberOfPoints() << std::endl;
		std::cout << " Total Pts=" << numPoints << std::endl;

		vtkIdType* pointIDs = new vtkIdType[numPoints];
		for (int i=0; i<numPoints; i++){
			pointIDs[i]=(vtkIdType)i+totalPoints;
		}  

		//put points in a line
		centerLine->InsertNextCell(numPoints,pointIDs); 
		totalPoints+=numPoints;
	}
  
  //write everything
  std::string outfile, prefix, suffix ="out.vtk";
  prefix=argv[2];
  prefix=prefix.substr(0,prefix.size()-4);
  outfile= prefix+suffix;   
  writePointsVTK(outfile.c_str(), centerPoints, centerLine, radiusArray); 

  std::cout << "File written" << std::endl;
  return 0;
}