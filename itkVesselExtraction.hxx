/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    itkVesselExtraction.hxx
Language:  C++
Date:      $Date$
Version:   $Revision$



=========================================================================*/
#ifndef __itkVesselExtraction_hxx
#define __itkVesselExtraction_hxx

#include "itkVesselExtraction.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

	/** Constructor */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap> 
	VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::VesselExtraction()
	{
		m_DistanceMap = 0;

		maxJ =5;

		// when set to true it will be a bit faster, but it will result in minimizing
		// the sum of distances^4 instead of the sum of distances^2
		m_ComputeSquaredDistance = false; 

		//W2012
		//parameters
		//scale =1; //defined as radius (not diameter) of vessel
		//beta= 0.4;
		//sigma =2.5;
		//sFactor=0.2;

		//other parameters
		//max_iterations=140;
		//epsilon =0.8;
		//gradTolerance=0.3;
		//windowSize =5;

		//seed model parameters
		//maxSeeds=100;
		//angThresh= 0.75;
		//distThresh= 3.0;
		//axialMean=4.0;	//in terms of spacing
		//axialSpread=2.5; //in terms of spacing
		//radialSpread=1.5;//in terms of spacing
		//radialIncrement=0.02; //the circle will move out
		////every iteration
		//axialIncrement=0.02; //

		//Hessian declarations
		/*HessianFilter = HessianFilterType::New();
		filter = FilterType::New();*/

		//Interpolator
		//interpolator1 = InterpolatorType1::New();

		//To be Objectified
		//angleThresh=0.58;
		//distanceThresh=1.0;

		//termPremature=false;

		//Gradient filter (for radius estimation)
		//typename gradientType::Pointer gradientFilter = gradientType::New();

		//interpolator for radius can be linear
		//typename InterpolatorType2::Pointer interpolator2 = InterpolatorType2::New();
		
		//Distances
			//distance_tolerance=5*spacing[0];
			//for (int i=0; i<20;i++)
				//distances[i]=0;	

		//std::cout << "VesselExtraction() ..." << std::endl;
	}

	/** Return the number of values, i.e the number of points in the moving set */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	unsigned int
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>  
		::GetNumberOfValues() const
	{
		MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

		if( !movingPointSet ) 
		{
			itkExceptionMacro( << "Moving point set has not been assigned" );
		}

		return  movingPointSet->GetPoints()->Size();
	}


	/** Get the match Measure */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	typename VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>::MeasureType
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::GetValue( const TransformParametersType & parameters ) const
	{
		MeasureType measure;
		measure.Fill(0.0);

		return measure;

	}

	/** Get the Derivative Measure */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::GetDerivative( const TransformParametersType & itkNotUsed(parameters),
		DerivativeType & itkNotUsed(derivative) ) const
	{

	}

	/** Get both the match Measure and theDerivative Measure  */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	void
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::GetValueAndDerivative(const TransformParametersType & parameters, 
		MeasureType & value, DerivativeType  & derivative) const
	{
		value = this->GetValue(parameters);
		this->GetDerivative(parameters,derivative);
	}



	/** PrintSelf method */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	void
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);

		os << indent << "SmartVesselExtraction:PrintSelf()"<< std::endl;

	}


	//
	////******************************************************************************
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::Test( int argc, char *argv[] )
	{
		std::cout << "VesselExtraction::Test()" << std::endl;

		//testVesselExtraction( );
		//return 0;

		return 0 ;
	}
	//******************************************************************************

	/**** VesselTree ****/


	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::read_file(char * argv[]){

			const char * filename ;//= argv[2];
			std::ifstream myReadFile;	
			std::string temp, temp2, temp3;
			PhysicalPointType directionPoint;	
			float mag=0;	

			myReadFile.open(filename);

			if (myReadFile.is_open() && myReadFile.good()) {		
				int i,j,num=0;
				float tempF[12], mag, mag2, mag3; //PointType is incompatible with sscanf

				VesselInfo iO = VesselInfo(); //information Object for extracting vessel
				//get number of points
				std::getline(myReadFile,temp);
				sscanf(temp.c_str(),"%i",&num);

				//read seeds		
				for (i=0;i<num;i++){
					mag=0;	mag2=0;	mag3=0;
					std::getline(myReadFile,temp);
					sscanf(temp.c_str(),"%i %f %f %f %i %i %f %f %f %f %f %f %f",&(iO.id),&tempF[0],&tempF[1],&tempF[2],&(iO.seedLine),&(iO.seedLineInd),&tempF[3],&tempF[4],&tempF[5],&tempF[6],&tempF[7],&tempF[8],&(iO.length)); 
					for (j=0;j<3;j++){
						if (j==2){
							iO.seed[j]=tempF[j];
							iO.seedU[j]=tempF[j]; //initialize
							iO.direction[j]=tempF[j+3];
							iO.endPoint[j]=tempF[j+6];
						}else{					
							iO.seed[j]=-1*tempF[j];
							iO.seedU[j]=-1*tempF[j]; //initialize
							iO.direction[j]=-1*tempF[j+3];					
							iO.endPoint[j]=-1*tempF[j+6];
						}
						iO.v_direction[j]=iO.direction[j]-iO.seed[j];
						iO.v_direction3[j]=iO.endPoint[j]-iO.seed[j];
						iO.radiusS=1*spacing[0]; //default value, later obtain initial estimate of radius
						mag+=pow(iO.v_direction[j],2);								
						mag3+=pow(iO.v_direction3[j],2);				
					}			
					mag/=sqrt(mag);
					mag3/=sqrt(mag3);
					iO.length=mag3;
					for (j=0;j<3;j++){
						iO.v_direction[j]/=mag;			
						iO.v_direction3[j]/=mag3;
					}						
					iO.isLeaf=-1; //all segments are leaves
					//std::cout << iO.id << "_" << iO.length << std::endl;

					//add info object
					infoVector.push_back(iO);			
				}
				myReadFile.close();

				//std::cout << "File opened: " << infoVector[0].seed<< " " << infoVector[0].direction <<" " << infoVector[0].length<< std::endl;
				return 0;

			} else {
				std::cout << "File not opened" << std::endl;
				return -1;
			}
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::read_data(char * argv[]){

			vtkPolyDataReader *reader = vtkPolyDataReader::New();
			reader->SetFileName( argv[2] ); //.vtk
			reader->Update();
			//Get the necessary data
			vtkPolyData *polyData =  reader->GetOutput(); 
			checkPts = polyData->GetPoints();
			checkLines = polyData->GetLines();
			
			//vtkIdType nLines = lines->GetNumberOfCells();
			//reader->Delete();

			return 0;

	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::generate_info( ){

			int i, j;
			int numPoints;
			vtkIdType pointId; //Id of a point
			float tLength; //"this" length
			float mag, magL;
			vtkSmartPointer<vtkIdList> pLine = vtkSmartPointer<vtkIdList>::New(); //the line gotten
			nLines = checkLines->GetNumberOfCells();
			PhysicalPointType gotPoint;
			
			VesselInfo iO = VesselInfo(); //information Object for extracting vessel
			checkLines->InitTraversal(); 

			for (i=0;i<nLines;i++){
				//set initial info
				iO.id=i;
				iO.seedLine=i+1;
				iO.seedLineInd=3;
				
				checkLines->GetNextCell(pLine); //get next line				
				tLength=0;
				iO.numPoints=pLine->GetNumberOfIds();
				
				//get seed
				get_point(pLine,0,&gotPoint);
				iO.seed=gotPoint;
				iO.seedU=gotPoint;				
				//get last point
				get_point(pLine,iO.numPoints-1,&gotPoint);
				iO.endPoint=gotPoint;				

				for (j=1;j<iO.numPoints;j++){
					
					get_point(pLine,j,&gotPoint);
					iO.direction=gotPoint;
					tLength=compute_distance(gotPoint,iO.seed);					

					if(tLength>=directionLength){
						break;
					}					
				}
				mag=0;
				iO.length=0;
				for (j=0;j<3;j++){
					iO.v_direction[j]=iO.direction[j]-iO.seed[j];
					mag+=pow(iO.v_direction[j],2);					
					iO.length+=pow(iO.endPoint[j]-iO.seed[j],2);					
				}
				mag/=sqrt(mag);
				iO.length=sqrt(iO.length);
				for (j=0;j<3;j++){
					iO.v_direction[j]/=mag;					
				}
				iO.radiusS=1*spacing[0]; //default value, later obtain initial estimate of radius				
				iO.isLeaf=-1; //all segments are leaves

				//add info object
				infoVector.push_back(iO);		
			}			

			return 0;
	}


	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::INITIALIZE(){			
			spacing	=image->GetSpacing();
			spacing[2]=spacing[0]; //normalize spacing

			//interpolator1
			//only used in Spline interpolation
			//interpolator1->SetSplineOrder(3);

			//magnitude gradient
			gradientType::Pointer gradientFilter = gradientType::New();
			gradientFilter->SetInput( image );
			gradientFilter->Update();

			//interpolator2
			gradImage=gradientFilter->GetOutput();
			interpolator2 = InterpolatorType2::New();
			interpolator2->SetInputImage(gradImage);

			//filter	
			filter->SetInput( image );
			filter->SetInterpolator( interpolator1 );
			direction.SetIdentity();  
			filter->SetDefaultPixelValue( 1.0 );  
			reSpacing=spacing*sFactor;

			const VectorType reSpacingConst(reSpacing);
			filter->SetOutputSpacing( reSpacingConst );  

			//filter->SetOutputSpacing( reSpacing );  
			filter->SetOutputDirection( direction );

			//Random generator
			randomGenerator = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
			randomGenerator->Initialize(2861);
			//2861


	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::tracking_LOOP(int n, vtkPoints *centerPoints, vtkFloatArray* radiusArray){
			int i=0,startI=1,j, goodPoints;
			int returnInt=0;
			int optimalSeedFound;
			MType E_mtx_prev, E_mtx;		
			PhysicalPointType physicalEstimation, physicalPoint=infoVector[n].seedU;	
			float radius_, radius_prev=infoVector[n].radiusS;  
			float distance_, distance__;
			float mag=0, cosAngle=0;
			bool endLoop=false;
			bool toContinue=false;
			termPremature=false;

			//***INITIALIZE***
			initialize_LOOP();

			//start_checks(n);
			if (firstTime){
				centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);
				radiusArray->InsertNextValue(radius_prev);
				startI=2;
				firstTime=false;
			}

			///Find a new seed
			optimalSeedFound=optimize_seed(infoVector[n].seedU, infoVector[n].v_direction, &(infoVector[n].newSeed));
			//and check if its good
			if (optimalSeedFound !=1){
				if (!firstTime){
					centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);
					radiusArray->InsertNextValue(radius_prev);
				} else {
					infoVector[n].oldNewSeed=infoVector[n].newSeed;
				}
				physicalPoint=infoVector[n].direction;
				centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);
				radiusArray->InsertNextValue(radius_prev);
				return 2;
			}	

			physicalPoint=infoVector[n].newSeed; //new_seed
			centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);
			radiusArray->InsertNextValue(radius_prev);
			track_distances(physicalPoint,infoVector[n].seedU); //or new_seed
			track_directions(physicalPoint);
			//std::cout << "PP (" << i << ") : " << physicalPoint << std::endl;	

			E_mtx_prev=compute_EigenMatrix(HessianFilter->GetOutput(), physicalPoint);
			E_mtx_prev=correct_orientation(E_mtx_prev,infoVector[n].v_direction);

			//***THE LOOP***
			for (i=startI; i<max_iterations && !endLoop; i++){ //2
				//**Find next cross section**
				E_mtx=compute_EigenMatrix(HessianFilter->GetOutput(), physicalPoint);

				//Ensure the next point has correct orientation
				E_mtx=change_sign(E_mtx_prev,E_mtx); //with respect to previous point
				//with respect to overall direction
				E_mtx=correct_orientation(E_mtx, infoVector[n].v_direction); 

				//**Estimate next point**
				physicalEstimation=move_Point(physicalPoint, E_mtx);
				
				if (write1)
					std::cout << "PE" << i << "," << convert_point(physicalEstimation) << std::endl;	

				//**Find the next point**
				physicalPoint=TwoD_Maximum(physicalEstimation, E_mtx);	
				//Record the point
				centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);

				if (write1)
					std::cout << "PP" << i << "," << convert_point(physicalPoint) << std::endl;	

				//**Find the radius**
				radius_=compute_radius(physicalPoint, radius_prev, spacing[0]*2, E_mtx);
				radiusArray->InsertNextValue(radius_);
				radius_prev=radius_;
				//scale=max(1,min(2,radius_/(2*spacing[0])));
				//std::cout << "Radius " << radius_ << std::endl;					
				//Save the cross-section and direction
				E_mtx_prev=E_mtx; 		

				distance__=compute_distance(physicalPoint,infoVector[n].seed);
				track_distances(physicalPoint,infoVector[n].seedU); //or new_seed
				track_directions(physicalPoint);

				//Determine to end
				if(distance__>infoVector[n].length){
					endLoop=true;
					std::cout << "L " << distance__ << "  " << infoVector[n].length << std::endl;	
				}	else if(termination(physicalPoint, E_mtx)){
					endLoop=true;
					termPremature=true;					
					//i=max(startI+1,i-10);		
					//try i=max(startI+1,i-1);		
				}				

				//if (check_check(physicalPoint,i)){
				//endLoop=true;
				//toContinue=true;
				//std::cout << "O " << std::endl;	
				//}

				//std::cout << "Dist " << distances[0] << std::endl;					
			}

			//if (!toContinue){
			goodPoints=(int)centerPoints->GetNumberOfPoints();
			//} else {
			//goodPoints=prevTrackedI;
			//}

			return i;
	}
	
		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::Vessel_tracking(int n, vtkPoints *smoothedPoints, vtkFloatArray* radiusArray)
	{
		firstTime=true;
		vtkPoints* linePoints = vtkPoints::New(); 
		//tracking_LOOP(n, smoothedPoints,  radiusArray); 		
		tracking_LOOP(n, linePoints,  radiusArray); 		
		smooth_points(linePoints,smoothedPoints);  //smoothen the points

		linePoints->Delete();
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::Vessel_tracking_orientation(int n, vtkPoints *trimmedPoints, vtkFloatArray* radiusArray)
	{
			vtkPoints *segmentPoints;
			vtkPoints *linePointsTemp;
			vtkPoints *linePoints = vtkPoints::New(); //keeps history of points
			vtkPoints *smoothedPoints;
			vtkPoints *trimmedSPoints;

			firstTime=true;
			int j=0;
			int continue_=0; //1 here means yes
			int numSPoints;
			int numAPoints;
			
			while (true){
								
				segmentPoints = vtkPoints::New(); 				
				smoothedPoints = vtkPoints::New();
				trimmedSPoints = vtkPoints::New();	
				linePointsTemp = vtkPoints::New();	
				write_points(linePoints,linePointsTemp); //fill the temp vector
					
				numSPoints=tracking_LOOP(n, segmentPoints, radiusArray); 			
				write_points(segmentPoints,linePointsTemp,numSPoints);
				
				smooth_points(linePointsTemp,smoothedPoints);  //smoothen the points			
				//write_points(linePointsTemp,smoothedPoints);
				if (n!=143){ //43
					continue_=check_line_orientiation(n, smoothedPoints, trimmedSPoints, termPremature);
				} else {
					write=true;
					if(j<2){ 
						continue_=check_line_orientiation(n, smoothedPoints, trimmedSPoints, termPremature);
					} else {
						continue_=check_line_orientiation(n, smoothedPoints, trimmedSPoints, termPremature);
						//write_points(linePointsTemp,trimmedSPoints);
						//write_points(smoothedPoints,trimmedSPoints);
						continue_=0;
						termPremature=false;
						
					}
				}					

				//std::cout << "continue_ " << continue_ << std::endl;
				//std::cout << "tP " << termPremature << std::endl;
				
				if (termPremature)
					continue_=1;
				
				//trim down the linePoints in accordance with check_line_orientation results
				numAPoints=(int)(trimmedSPoints->GetNumberOfPoints()-linePoints->GetNumberOfPoints());
				trim_equal(segmentPoints,linePoints,numAPoints); 				
				
				std::cout << "AT: " << j << " #Pts=" << numSPoints <<  "; " << linePoints->GetNumberOfPoints() << std::endl;

				segmentPoints->Delete();
				smoothedPoints->Delete();		
				linePointsTemp->Delete();
				
				j++; 

				if (continue_!=1 || j>=maxJ){ //maxJ
					write_points(trimmedSPoints, trimmedPoints); //linePoints are all the points
					trimmedSPoints->Delete();
					break;
				}	else{
					trimmedSPoints->Delete();
				}
			}
			//include in any final version but comment out line above
			
			//smooth_points(linePoints,smoothedPoints);  //smoothen the points
			//infoVector[n].newSeed=infoVector[n].oldNewSeed;
			//check_line_orientiation(n, smoothedPoints, trimmedPoints); 
		
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::Extraction_MAIN(int argc, char * argv[])
	{
		std::cout << " VesselExtraction::VesselExtraction()... " << std::endl;
		if(argc != 3)
		{
			std:: cout << "Usage: SmartVesselExtraction <image filename> <lines filename>" << std::endl;
			std:: cout << "Example: SmartVesselExtraction EHPS3R3Crop.mhd EHPS2103BP806Tfm3.vtk" << std::endl;
			exit(-1);
		} 	

		/************Declarations***********
		*/
		typedef itk::ImageFileReader< ImageType >  ReaderType;
		vtkPoints *centerPoints = vtkPoints::New(); 		
		vtkCellArray *centerLine = vtkCellArray::New();
		vtkFloatArray* radiusArray = vtkFloatArray::New();

		typename MType E_mtx;  
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
		INITIALIZE();		
		radiusArray->SetNumberOfComponents(1);
		radiusArray->SetName("RadiusArray");  
		std::cout << "Initialized ";

		if (read_data(argv)==-1){
			exit(-1);
		} 
		//Get seed and target points
		generate_info();

		std::cout << " read input " << infoVector.size() <<  std::endl;		
		
		//get_seed_target(argv);   
		/*
		start_checks(0);
		exit(-1);*/

		//Get Hessian
		std::cout << "Processing Hessian ...";
		find_Hessian();  
		std::cout << "done" << std::endl;

		////EXTRACT VESSEL TREE////

		totalPoints=0;
		for (int n=0; n<infoVector.size();n++){ // infoVector.size();n++){ 
			vtkPoints *trimmedPoints = vtkPoints::New();			
			std::cout << "-------N " << n << std::endl;		

			//if (n==43)
				//write=true;
			
			///////////// Track vessel  /////////////
			
			//Vessel_tracking(n,trimmedPoints,radiusArray);
			Vessel_tracking_orientation(n,trimmedPoints,radiusArray);			

			////Creating output////

			//find total number of points
			numPoints=(int)trimmedPoints->GetNumberOfPoints();
			std::cout << " Total Pts=" << numPoints << std::endl;

			//add them to the array of all center points
			write_points(trimmedPoints, centerPoints);
			
			//std::cout << "Actual Points =      " << (int)smoothedPoints->GetNumberOfPoints() << std::endl;
			
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
		write_pointsVTK(outfile.c_str(), centerPoints, centerLine, radiusArray); 

		std::cout << "File written" << std::endl;

		std::cout << "Times: 1: " << time1 << "s 2: " << time2 << "s " << std::endl;

		return 0;
	}



} // end namespace itk
#endif



