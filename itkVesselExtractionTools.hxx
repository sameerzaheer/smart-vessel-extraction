/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    itkVesselExtractionTools.hxx
Language:  C++
Date:      $Date$
Version:   $Revision$



=========================================================================*/
#ifndef __itkVesselTools_hxx
#define __itkVesselTools_hxx

#include "itkVesselExtractionTools.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

	/** Constructor */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap> 
	VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::VesselExtractionTools()
	{
		m_DistanceMap = 0;

		// when set to true it will be a bit faster, but it will result in minimizing
		// the sum of distances^4 instead of the sum of distances^2
		m_ComputeSquaredDistance = false; 


		//W2012
		
		//parameters
		scale =1; //defined as radius (not diameter) of vessel
		//scale=3;
		float betaSlow= 0.4;
		float betaFast = 0.8;
		sigma =2.5;
		sFactor=0.2;  //Need to increase this
		sFactor_=1;
		intensityThreshold=0.2;

		//other parameters
		max_iterations=140;
		epsilon =0.8;
		gradTolerance=0.3;
		windowSize =5;

		//seed model parameters
		maxSeeds= 50;
		angThresh= 0.75;
		distThresh= 3.0;
		axialMean= 4.0;	//in terms of spacing
		axialSpread= 2.5; //in terms of spacing
		radialSpread= 1.5; //in terms of spacing
		radialIncrement= 0.02; //the circle will move out
		//every iteration
		axialIncrement= 0.02; //

		//Hessian declarations
		HessianFilter = HessianFilterType::New();
		filter = FilterType::New();
		//filter->Print(std::cout);

		//Interpolator
		interpolator1 = InterpolatorType1::New();

		//To be Objectified
		directionLength=3; //in mm

		angleThresh=0.54;
		distanceThresh=1.0;

		termPremature=false;

		//Gradient filter (for radius estimation)
		typename gradientType::Pointer gradientFilter = gradientType::New();

		//interpolator for radius can be linear
		typename InterpolatorType2::Pointer interpolator2 = InterpolatorType2::New();

		//std::cout << "VesselExtractionTools() ..." << std::endl;

		//write
		write=false;
		write1=false;
		write2=false;
		write3=false;
		write4=false;
		write5=false;
		write6=false; //writing cross-section data to file
		write7=false; //writing radius data to file

		//time
		time1=0;
		time2=0;
		time3=0;

		fastSwitch=false; //for old or slow tracking
		//fastSwitch=true; //for fast tracking or new tracking, may be unreliable

		if (fastSwitch){
			beta=betaFast;
		} else {
			beta=betaSlow;
			beta=0.6;
		}
	
	}

	/** Return the number of values, i.e the number of points in the moving set */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	unsigned int
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>  
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
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MeasureType
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::GetValue( const TransformParametersType & parameters ) const
	{
		MeasureType measure;
		measure.Fill(0.0);

		return measure;

	}

	/** Get the Derivative Measure */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::GetDerivative( const TransformParametersType & itkNotUsed(parameters),
		DerivativeType & itkNotUsed(derivative) ) const
	{

	}

	/** Get both the match Measure and theDerivative Measure  */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	void
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::GetValueAndDerivative(const TransformParametersType & parameters, 
		MeasureType & value, DerivativeType  & derivative) const
	{
		value = this->GetValue(parameters);
		this->GetDerivative(parameters,derivative);
	}



	/** PrintSelf method */
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
	void
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);

		os << indent << "SmartVesselExtraction:PrintSelf()"<< std::endl;

	}


	//
	////******************************************************************************
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::Test( int argc, char *argv[] )
	{
		std::cout << "VesselExtractionTools::Test()" << std::endl;

		//testVesselExtraction( );
		//return 0;

		return 0 ;
	}
	//******************************************************************************


	/**** VesselTree ****/
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::min(float a, float b)
	{
		if (a>b){
			return b;
		}else{
			return a;
		}
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::max(float a, float b){
			if (a<b){
				return b;
			}else{
				return a;
			}
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::smallest_index(VectorType A){
			int i, index=0;
			float min=abs(A[index]);

			for (i=index+1;i<3;i++){
				if	(abs(A[i])<min){
					index=i;
					min=abs(A[i]);
				}
			}
			return index;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::VectorType
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_vector(PhysicalPointType pointA, PhysicalPointType pointB){
			VectorType returnVector;
			float mag=0;

			for (int i=0;i<3;i++){
				returnVector[i]=pointB[i]-pointA[i];
				mag+=pow(returnVector[i],2);
			}

			for (int i=0;i<3;i++){
				returnVector[i]/=sqrt(mag);
			}
			return returnVector;
	}

	//return pointA+vectorB
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_point(PhysicalPointType pointA, VectorType vectorB){

			PhysicalPointType returnPoint;
			
			for (int i=0;i<3;i++){
				returnPoint[i]=pointA[i]+vectorB[i];
			}

			return returnPoint;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::convert_point(PhysicalPointType pointIn){

			PhysicalPointType pointOut;
			
			pointOut[0]=pointIn[0]*-1;
			pointOut[1]=pointIn[1]*-1;
			pointOut[2]=pointIn[2];

			return pointOut;
	}
	
		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::VectorType
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::convert_vector(VectorType vectorIn){

			PhysicalPointType vectorOut;
			
			vectorOut[0]=vectorIn[0]*-1;
			vectorOut[1]=vectorIn[1]*-1;
			vectorOut[2]=vectorIn[2];

			return vectorOut;
	}

			template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	double
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::vector_length(VectorType vectorIn){

			double _length=0;
			
			for (int i=0; i<3;i++){
				_length+=vectorIn[i]*vectorIn[i];
			}
			
			_length=sqrt(_length);

			return _length;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::cross_product(VectorType A, VectorType B, VectorType* C){
			(*C)[0]=(*C)[1]=(*C)[2]=0;
			//calculate the i,j,k coefficients of the cross product 
			(*C)[0]=A[1]*B[2]-B[1]*A[2];
			(*C)[1]=-A[0]*B[2]+B[0]*A[2];
			(*C)[2]=A[0]*B[1]-B[0]*A[1];
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_planeV(VectorType C, VectorType* A, VectorType* B){
			int i=0;
			VectorType arb; 

			//create an arbitrary vector
			arb[0]=arb[1]=arb[2]=0; //initialize
			arb[smallest_index(C)]=1; //to catch a singularity

			cross_product(C,arb,A); //compute A
			cross_product(C,*A, B); //compute B
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_distance(PhysicalPointType curPoint, PhysicalPointType refPoint){
			float distance_;
			distance_=sqrt(pow(refPoint[0]-curPoint[0],2)+pow(refPoint[1]-curPoint[1],2)
				+pow(refPoint[2]-curPoint[2],2));
			return distance_;
	}

		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_angle(PhysicalPointType A, PhysicalPointType B, PhysicalPointType C){
			VectorType BA, BC;
			float magBA=0, magBC=0;
			float dotProd=0;

			for (int i=0;i<3;i++){
				BA[i]=A[i]-B[i];
				BC[i]=C[i]-B[i];
				magBA+=pow(BA[i],2);
				magBC+=pow(BC[i],2);
				dotProd+=BA[i]*BC[i];
			}	
			if (magBA*magBC>0)
				return dotProd/sqrt(magBA*magBC);
			else
				return -2; //range would normally be between -1 and 1
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_angle(PhysicalPointType A, PhysicalPointType B, PhysicalPointType C, PhysicalPointType D){
			VectorType BA, CD;
			float magBA=0, magCD=0;
			float dotProd=0;
			for (int i=0;i<3;i++){
				BA[i]=A[i]-B[i];
				CD[i]=D[i]-C[i];
				magBA+=pow(BA[i],2);
				magCD+=pow(CD[i],2);
				dotProd+=BA[i]*CD[i];
			}	
			if (magBA*magCD>0)
				return dotProd/sqrt(magBA*magCD);
			else
				return -2;  //range would normally be between -1 and 1
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::calculate_Gradient(IndexType point){
			ImageType::PixelType pVx, pVy, pVz, qVx, qVy, qVz;	
			MType G(3,1);

			point[0]++; pVx= image->GetPixel(point); point[0]--; point[0]--; qVx= image->GetPixel(point); point[0]++;
			point[1]++; pVy= image->GetPixel(point); point[1]--; point[1]--; qVy= image->GetPixel(point); point[1]++;
			point[2]++; pVz= image->GetPixel(point); point[2]--; point[2]--; qVz= image->GetPixel(point); point[2]++;
			G[0][0]=(pVx-qVx)/(2*spacing[0]);
			G[1][0]=(pVy-qVy)/(2*spacing[1]);
			G[2][0]=(pVz-qVz)/(2*spacing[2]);

			return G;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::OneD_Maximum(PhysicalPointType seedPoint, PhysicalPointType *oneDMax){

			IndexType seedIndex;
			image->TransformPhysicalPointToIndex(seedPoint, seedIndex);
			IndexType imageIndex=seedIndex;

			//Calculate gradient
			ImageType::PixelType mag, pVx, pVy, pVz;
			MType G=calculate_Gradient(seedIndex);
			pVx=G[0][0]; pVy=G[1][0]; pVz=G[2][0];
			//std::cout << "Gradient: " << pVx-pixelValue << ',' << pVy-pixelValue << ',' << pVz-pixelValue <<std::endl; 

			if ((pVx+pVy+pVz)==0){
				//std::cout << "Choose another seed point" <<std::endl;		
				return 0;
			}

			//Determine end point of search
			IndexType endPoint;
			mag=sqrt(pow(pVx,2)+pow(pVy,2)+pow(pVz,2)); 
			endPoint[0]=seedIndex[0]+6*scale*(pVx)/mag; //do a very long search, in case the seed is off
			endPoint[1]=seedIndex[1]+6*scale*(pVy)/mag; //searches along 3 times the diameter of the vessel
			endPoint[2]=seedIndex[2]+6*scale*(pVz)/mag;
			//std::cout << "Far point: " << endPoint[0] << ',' << endPoint[1] << ',' << endPoint[2] <<std::endl; 

			// Find maximum 
			ImageLineItType it(image, seedIndex, endPoint);
			it.GoToBegin();
			PixelType max=-1;
			PixelType temp=0;
			while (!it.IsAtEnd())
			{
				temp=(PixelType)it.Get();
				if (temp>max){
					max=temp;
					imageIndex=it.GetIndex();
				}
				++it;
			}

			image->TransformIndexToPhysicalPoint(imageIndex,*oneDMax);
			return 1;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::find_Hessian(){

			ImageType::RegionType selectRegion;		
						
			HessianFilter->SetInput(image);						
			HessianFilter->SetSigma(sigma);			
			HessianFilter->SetNormalizeAcrossScale(true);   

			HessianFilter->Update();	
			
			//return 0;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_EigenMatrix(HessianImagePointer HessianImage, PhysicalPointType physicalPoint){
			MType T(4,3); //matrix holding both eigenvalues and eigenvectors

			IndexType imageIndex;
			image->TransformPhysicalPointToIndex(physicalPoint, imageIndex);

			HessianImageLineItType Hess_It(HessianImage, imageIndex, imageIndex);
			Hess_It.GoToBegin();
			itk::SymmetricSecondRankTensor< float, 3> HM = Hess_It.Get();
			float HMatrix[9]={	HM[0], HM[1], HM[2],
				HM[1], HM[3], HM[4],
				HM[2], HM[4], HM[5], };

			MType S(HMatrix, 3,3);
			EigenVAType evals; EigenVMType evectors;
			EigenAType ESys(3);
			ESys.SetOrderEigenValues(true);
			ESys.ComputeEigenValuesAndVectors(S, evals, evectors );

			T[0][0]=evectors[0][0]; T[0][1]=evectors[0][1]; T[0][2]=evectors[0][2];
			T[1][0]=evectors[1][0]; T[1][1]=evectors[1][1]; T[1][2]=evectors[1][2];
			T[2][0]=evectors[2][0]; T[2][1]=evectors[2][1]; T[2][2]=evectors[2][2];
			T[3][0]=evals[0]; T[3][1]=evals[1]; T[3][2]=evals[2];

			if (write2){
				std::cout<< "Evecs" << std::endl;
				std::cout<<"["<<T[0][0]<<","<<T[0][1]<<","<<T[0][2]<<"];"<<std::endl;
				std::cout<<"["<<T[1][0]<<","<<T[1][1]<<","<<T[1][2]<<"];"<<std::endl;
				std::cout<<"["<<T[2][0]<<","<<T[2][1]<<","<<T[2][2]<<"];"<<std::endl;
				std::cout<< "Evals" << std::endl;
				std::cout<<"["<<T[3][0]<<","<<T[3][1]<<","<<T[3][2]<<"];"<<std::endl;
			}

			return T;
	}


	//Computes "eigen-matrix" with specified vector as the direction, and computes a
	//plane normal to it
		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_EigenMatrix(VectorType vector){
			MType T(4,3); //matrix holding both eigenvalues and eigenvectors

			T[2][0]=vector[0]; T[2][1]=vector[1]; T[2][2]=vector[2];

			VectorType A, B;
			compute_planeV(vector, &A, &B);

			T[0][0]=A[0]; T[0][1]=A[1]; T[0][2]=A[2];
			T[1][0]=B[0]; T[1][1]=B[1]; T[1][2]=B[2];
			T[3][0]=0; T[3][1]=0; T[3][2]=1;

			return T;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::make_EigenMatrix(VectorType vA, VectorType vB, VectorType vC){
			MType T(4,3); //matrix holding both eigenvalues and eigenvectors

			for (int i=0;i<3;i++){		
				T[0][i]=vC[i]; 
				T[1][i]=vB[i];
				T[2][i]=vA[i];
			}
			T[3][0]=-1; T[3][1]=-1; T[3][2]=1;

			return T;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::TwoD_Maximum(PhysicalPointType physicalIndex,  MType E_matrix){

			//time
			//clock_t start1T, start2T, start3T, final1T, final2T, final3T;

			//start1T=clock();

			//********Resample image********
			IndexType pointIndex;
			PhysicalPointType maxPhysicalIndex;		

			float reScale = scale*spacing[0]; //physical radius of vessel
			origin[0] = physicalIndex[0]-reScale-reSpacing[0]; //reSpacing[0]->extra padding
			origin[1] = physicalIndex[1]-reScale-reSpacing[0];
			origin[2] = physicalIndex[2]-reScale-reSpacing[0];
			filter->SetOutputOrigin( origin );	 

			SizeType   rSize;
			rSize[0]=reScale*2/reSpacing[0]+2; //re_scale/reSpacing[0]=scale in voxels
			rSize[1]=reScale*2/reSpacing[0]+2; //+2 to deal with rounding errors
			rSize[2]=reScale*2/reSpacing[0]+2;
			filter->SetSize( rSize );  

			//final1T=clock();
			//start2T=clock();

			filter->Update();  
			resampledImage=filter->GetOutput();
			resampledImage->TransformPhysicalPointToIndex(physicalIndex, pointIndex);				
			
			//std::cout << "rSize   : "  << rSize<< std::endl;
			//std::cout << "re size   : "  << resampledImage->GetLargestPossibleRegion().GetSize()<< std::endl;

			//final2T=clock();
			//start3T=clock();

			//********Extract eigen vectors********
			IndexType pointA1, pointB1, pointA2, pointB2, currentIndex, maxIndex, AIndex, BInd1, BInd2;
			PhysicalPointType PpointA1, PpointB1, PpointA2, PpointB2;
			float A[3], B[3]; 
			float mag;
			float hSize=reScale/(reSpacing[0]*sqrt(float(2)));			
			//divide by root(2) because that is max sum of two orthogonal vectors

			A[0]=E_matrix[0][0];

			for (int i=0; i<3;i++){
				A[i]=E_matrix[0][i]; 
			}

			mag=sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2)); 
			for (int i=0; i<3;i++){
				A[i]=A[i]/mag; 
			}

			for (int i=0; i<3;i++){
				B[i]=E_matrix[1][i]; 
			}
			mag=sqrt(pow(B[0],2)+pow(B[1],2)+pow(B[2],2)); 
			for (int i=0; i<3;i++){
				B[i]=B[i]/mag; 
			}

			/*if (write2){
				std::cout << "PI : " << convert_point(physicalIndex) << std::endl;	
				std::cout << "A->" << A[0] << "," << A[1] << "," << A[2] << "; B->" << B[0] << "," << B[1] << "," << B[2] << std::endl;
			}*/

			//Iterate in direction A

			for (int i=0;i<3;i++){
				pointA1[i]=min(max(floor(pointIndex[i]-hSize*A[i]+0.5),0),rSize[i]-1); 				
				pointA2[i]=min(max(floor(pointIndex[i]+hSize*A[i]+0.5),0),rSize[i]-1);				
			}

			/*if (write2){
				resampledImage->TransformIndexToPhysicalPoint(pointA1, PpointA1);
				resampledImage->TransformIndexToPhysicalPoint(pointA2, PpointA2);
				std::cout << "hasz " << hSize << std::endl;
				std::cout << "A1 " << convert_point(PpointA1) << std::endl;
				std::cout << "A2 " << convert_point(PpointA2) << std::endl;
			}*/
	
			ImageLineItType itA(resampledImage, pointA1, pointA2); //iterator in directionA
			itA.GoToBegin();
			PixelType currentValue, maxValue=-1;
			currentValue=(PixelType)itA.Get();	
			maxIndex=pointIndex;

			//temp
			bool t=false;
			if (write2){
				t=true;
			}

			while (!itA.IsAtEnd())
			{
				//Determine points of traversing in the B direction
				currentIndex=itA.GetIndex();		

				for (int i=0;i<3;i++){
					pointB1[i]=min(max(floor(currentIndex[i]-hSize*B[i]+0.5),0),rSize[i]-1); 					
					pointB2[i]=min(max(floor(currentIndex[i]+hSize*B[i]+0.5),0),rSize[i]-1);					
				}

				//std::cout << "cI " << currentIndex << std::endl;
				/*if (t){
					t=false;
					resampledImage->TransformIndexToPhysicalPoint(pointB1, PpointB1);
					resampledImage->TransformIndexToPhysicalPoint(pointB2, PpointB2);
					std::cout << "B1 " << convert_point(PpointB1) << std::endl;
					std::cout << "B2 " << convert_point(PpointB2) << std::endl;
				}*/

				ImageLineItType itB(resampledImage, pointB1, pointB2);				

				//Iterate in directionB
				while (!itB.IsAtEnd())
				{
					currentValue=(PixelType)itB.Get();
					if (currentValue>maxValue){
						maxValue=currentValue;
						maxIndex=itB.GetIndex();
						AIndex=currentIndex;
						BInd1=pointB1;
						BInd2=pointB2;
					}			

					++itB;
				}
				++itA;
			}	
			//time
			//final3T=clock();

			//time1+=(double)(final1T-start1T)/CLOCKS_PER_SEC;	
			//time2+=(double)(final2T-start2T)/CLOCKS_PER_SEC;
			//time3+=(double)(final3T-start3T)/CLOCKS_PER_SEC;

			resampledImage->TransformIndexToPhysicalPoint(AIndex, maxPhysicalIndex);			
			//std::cout << "A ind  " << maxPhysicalIndex << std::endl;	

			resampledImage->TransformIndexToPhysicalPoint(maxIndex, maxPhysicalIndex);			

			//std::cout << "| " << physicalIndex << ". Max " << maxPhysicalIndex << "  " << maxIndex << std::endl;	

			return maxPhysicalIndex;
	}

		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
		float 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::Ray_Boundary(PixelType **CS, int *CSsize, int *point, float degree){

			//Ray is the degree position, where 0 is 12'o clock and degree increase in the clockwise direction
			//A is 12'o clock, and B is 3 o'clock

			PixelType currentGrad, maxGrad=0;
			PixelType currentValue, previousValue;
			int maxI=0;
			degree=degree*PI/180;
			int dA=floor(cos(degree)+0.5);//the equivalent of rounding
			int dB=floor(sin(degree)+0.5);
			float dL=sqrt(pow((float)dA,2)+pow((float)dB,2));

			//Iterate in directionC
			//std::cout << "Del " << dA << " " << dB << std::endl;	
			//std::cout << "CSsize " << CSsize[0] << " " << CSsize[1] << std::endl;	
			//std::cout << "Gradients" << std::endl;	

			int a=point[0];
			int b=point[1];
			previousValue=CS[a][b];

			while(a>=0 && a<CSsize[0] && b>=0 && b<CSsize[1]){
				currentValue=CS[a][b];
				currentGrad=(currentValue-previousValue)/dL;
				//std::cout << "a:" << a << " b:" << b << ", " << currentGrad << std::endl;	
				previousValue=currentValue;
				if (currentGrad<maxGrad){
					maxGrad=currentGrad;
					maxI=max(abs(point[0]-a),abs(point[1]-b));
				}
				a+=dA;
				b+=dB;
			}

			//std::cout << "rad: " << maxI << std::endl;	
			//What remains is to set up the test for this

			if (maxGrad<-30){
				return maxI;		
			} else{
				return 0;
			}
	}

		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
		void
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::CS_center(PixelType **CS, int *CSsize, int searchBound, int *Center){
			int center[2]={floor((float)CSsize[0]/2),floor((float)CSsize[1]/2)};
			int searchA1=max(0,center[0]-searchBound);
			int searchA2=min(CSsize[0],center[0]+searchBound);
			int searchB1=max(0,center[1]-searchBound);
			int searchB2=min(CSsize[1],center[1]+searchBound);

			float weight1=0.9, weight2=1-weight1;

			int rayOrigin[2];

			float ray1,ray2,ray3,ray4,ray5,ray6,ray7,ray8,numPairRays;
			float centerLikelihood, maxLikelihood=0;
			Center[0]=center[0]; Center[1]=center[1]; //default

			////Make a Gaussian smoothed cross-section////
			float **sSection= new float*[CSsize[0]];
			for (int i=0;i<CSsize[0];i++){
					sSection[i]= new float[CSsize[1]];
			}

			float gaussAver;
			float weightSum, weights;
			float sigK=1,sigL=1;
			float maxVal=0, minVal=CS[0][0];
			for (int i=0;i<CSsize[0];i++){
				for (int j=0;j<CSsize[1];j++){
					gaussAver=0;
					weightSum=0;
					for (int k=-1;k<2;k++){
						for (int l=-1;l<2;l++){
							if ((i+k)>0 && (j+l)>0 && (i+k)<CSsize[0] && (j+l)<CSsize[1]){
								weights=exp(-1*pow((float)k,2)/pow(sigK,2)-1*pow((float)l,2)/pow(sigL,2));		
								weightSum+=weights;
								gaussAver+=weights*CS[i+k][j+l];
							}
						}
					}	
					sSection[i][j]=gaussAver/weightSum;
					if(maxVal<sSection[i][j]){
						maxVal=sSection[i][j];
					}
					if(minVal>sSection[i][j]){
						minVal=sSection[i][j];
					}

				}
			}

			if (write4){
				for (int i=0;i<CSsize[0];i++){
					for (int j=0;j<CSsize[1];j++){
						std::cout << std::fixed<< setprecision(0)<< setw(7) << floor(1000*CS[i][j]/maxVal) ;
					}
					std::cout << std::endl;
				}
				std::cout << std::endl;

				for (int i=0;i<CSsize[0];i++){
					for (int j=0;j<CSsize[1];j++){
						std::cout << std::fixed<< setprecision(0)<< setw(7) << floor(1000*(sSection[i][j]-minVal)/(maxVal-minVal)) ;
					}
					std::cout << std::endl;
				}
			}
		

			float **cSection= new float*[CSsize[0]];
			for (int i=0;i<CSsize[0];i++){
					cSection[i]= new float[CSsize[1]];
			}

			if (write4)		
				std::cout << "Radii" << std::endl;

			for (int i=0;i<CSsize[0];i++){
				for (int j=0;j<CSsize[1];j++){
					numPairRays=0;
					if (i>=searchA1 && i<searchA2 && j>=searchB1 && j<searchB2){

						rayOrigin[0]=i;
						rayOrigin[1]=j;

						centerLikelihood=0;

						ray1=Ray_Boundary(sSection, CSsize, rayOrigin, 0);	
						ray2=Ray_Boundary(sSection, CSsize, rayOrigin, 180);
						if(max(ray1,ray2)==0) //account for case where both are 0
							ray1=1;
						centerLikelihood+= 1-abs(ray1-ray2)/(ray1+ray2); numPairRays++;
						float cLl1=1-abs(ray1-ray2)/(ray1+ray2);

						ray3=Ray_Boundary(sSection, CSsize, rayOrigin, 90);	
						ray4=Ray_Boundary(sSection, CSsize, rayOrigin, 270);
						if(max(ray3,ray4)==0)
							ray3=1;
						centerLikelihood+= 1-abs(ray3-ray4)/(ray3+ray4); numPairRays++;
						float cLl2=1-abs(ray3-ray4)/(ray3+ray4);

						ray5=Ray_Boundary(sSection, CSsize, rayOrigin, 45);	
						ray6=Ray_Boundary(sSection, CSsize, rayOrigin, 225);
						if(max(ray5,ray6)==0)
							ray5=1;
						centerLikelihood+=1-abs(ray5-ray6)/(ray5+ray6); numPairRays++;
						float cLl3=1-abs(ray5-ray6)/(ray5+ray6);

						ray7=Ray_Boundary(sSection, CSsize, rayOrigin, 135);	
						ray8=Ray_Boundary(sSection, CSsize, rayOrigin, 315);
						if(max(ray7,ray8)==0)
							ray7=1;
						centerLikelihood+=1-abs(ray7-ray8)/(ray7+ray8); numPairRays++;
						float cLl4=1-abs(ray7-ray8)/(ray7+ray8);
						
						//cSection[i][j]= centerLikelihood/2; 
						//the weighted average of the centerlikelihood and distance from the center
						cSection[i][j]= weight1*(centerLikelihood/numPairRays)
							+weight2*(1/(  1+sqrt( pow((float)i-center[0],2)+pow((float)j-center[1],2) )  )); 
						if (maxLikelihood<cSection[i][j]){
							maxLikelihood=cSection[i][j];							
							Center[0]=i-center[0]; Center[1]=j-center[1];
							
						}						
					
						if (write4){
							std::cout << setw(1) <<  setprecision(0) << ray1 <<ray2 <<ray3 <<ray4 <<ray5 <<ray6 <<ray7 <<ray8 << ",";
							//std::cout << setprecision(1) << cLl1 << " " << cLl2 << " " << cLl3 << " " << cLl4;
						}

							//std::cout << setw(9) <<(ray1%10)*10000000 +(ray2%10)*1000000 +(ray3%10)*100000 +ray4%10*10000 +(ray5%10)*1000 +(ray6%10)*100 +(ray7%10)*10 +ray8%10;

					} else{
						cSection[i][j]=0;

						if (write4)
							std::cout << setw(9) <<  0 ;
					}
				}
				if (write4)
					std::cout << std::endl;
			}

			if (write4){
				std::cout << "Likelihood" << std::endl;

				for (int i=0;i<CSsize[0];i++){
					for (int j=0;j<CSsize[1];j++){
						if (cSection[i][j]==0){
							std::cout << setw(7) << setprecision(0) << cSection[i][j] ;
						} else {
							std::cout << setw(7) << setprecision(2) << cSection[i][j] ;
						}
					}
					std::cout << std::endl;
				}
				
			}

			for (int i=0;i<CSsize[0];i++){
					delete [] cSection[i];
			}
			delete [] cSection;			
			//return 0;
		}

		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::TwoD_Recenter(PhysicalPointType physicalIndex,  MType E_matrix){

			//time
			clock_t start1T, final1T;
			start1T=clock();

			//********Sizes********
			IndexType pointIndex;
			PhysicalPointType maxPhysicalIndex, maxPhysicalIIndex;		

			float reScale = scale*spacing[0]; //physical radius of vessel

			SizeType   iSize;
			iSize=image->GetLargestPossibleRegion().GetSize();
			image->TransformPhysicalPointToIndex(physicalIndex, pointIndex);

			int CSsize[2], CSa=0, CSb=0; //sizes of the cross-section matrix
			PixelType **cSection;
			IndexType **cSectionI;
			int Center2D[2];

			//********Extract eigen vectors********
			IndexType pointA1, pointB1, pointA2, pointB2, currentIndex, maxIndex, maxIIndex, AIndex, BInd1, BInd2;
			IndexType pointTemp;
			PhysicalPointType PpointTemp;
			PhysicalPointType PpointA1, PpointB1, PpointA2, PpointB2;
			float A[3], B[3]; 
			float mag;
			//float hSize=reScale/(reSpacing[0]*sqrt(float(2)));
			float hSize=2*reScale/(reSpacing_[0]);
			hSize=hSize*4;
			//std::cout << "hSize: " << hSize << std::endl;
			//divide by root(2) because that is max sum of two orthogonal vectors

			for (int i=0; i<3;i++){
				A[i]=E_matrix[0][i]; 
			}
			mag=sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2)); 
			for (int i=0; i<3;i++){
				A[i]=A[i]/mag; 
			}

			for (int i=0; i<3;i++){
				B[i]=E_matrix[1][i]; 
			}
			mag=sqrt(pow(B[0],2)+pow(B[1],2)+pow(B[2],2)); 
			for (int i=0; i<3;i++){
				B[i]=B[i]/mag; 
			}
			
			//write2=true;
			if (write2){
				std::cout << "PI : " << convert_point(physicalIndex) << std::endl;	
				std::cout << "I : " << pointIndex << std::endl;	
				std::cout << "A->" << A[0] << "," << A[1] << "," << A[2] << "; B->" << B[0] << "," << B[1] << "," << B[2] << std::endl;				
			}
			if (write5){
				std::cout << "A->" << A[0] << "," << A[1] << "," << A[2] << "; B->" << B[0] << "," << B[1] << "," << B[2] << std::endl;				
			}

			//Iterate in direction A

			for (int i=0;i<3;i++){		
				pointA1[i]=min(max(floor(pointIndex[i]-hSize*A[i]+0.5),0),iSize[i]-1);				
				pointA2[i]=min(max(floor(pointIndex[i]+hSize*A[i]+0.5),0),iSize[i]-1);				
			}

			//mock calculation of B
			for (int i=0;i<3;i++){					
					pointB1[i]=min(max(floor(pointA1[i]-hSize*B[i]+0.5),0),iSize[i]-1);
					pointB2[i]=min(max(floor(pointA1[i]+hSize*B[i]+0.5),0),iSize[i]-1);
			}

			//Compute the dimensions of the CS

			CSsize[0]=max(abs(pointA2[0]-pointA1[0])+1,max(abs(pointA2[1]-pointA1[1])+1,abs(pointA2[2]-pointA1[2])+1));
			CSsize[1]=max(abs(pointB2[0]-pointB1[0])+1,max(abs(pointB2[1]-pointB1[1])+1,abs(pointB2[2]-pointB1[2])+1));
			//std::cout << "Sizes: " << CSsize[0] << " " << CSsize[1] << std::endl;
			//Allocate space for matrix CS
			cSection= new float*[CSsize[0]];
			cSectionI= new IndexType*[CSsize[0]];
			for (int i=0;i<CSsize[0];i++){
					cSection[i]= new float[CSsize[1]];
					cSectionI[i]= new IndexType[CSsize[1]];
			}
			
			if (write2){
				image->TransformIndexToPhysicalPoint(pointA1, PpointA1);
				image->TransformIndexToPhysicalPoint(pointA2, PpointA2);
				std::cout << "hasz " << hSize << std::endl;
				std::cout << "A1 " << convert_point(PpointA1) << " " << pointA1 << std::endl;
				std::cout << "A2 " << convert_point(PpointA2) << " " << pointA2 << std::endl;
			}
	
			ImageLineItType itA(image, pointA1, pointA2); //iterator in directionA
			itA.GoToBegin();
			PixelType currentValue, maxValue=-1;
			currentValue=(PixelType)itA.Get();	
			maxIndex=pointIndex;

			//temp
			bool t=false;
			if (write2){
				t=true;
			}

			minIntensity=10000;
			maxIntensity=0;


			while (!itA.IsAtEnd())
			{
				//Determine points of traversing in the B direction
				currentIndex=itA.GetIndex();		
				//std::cout << "i " << currentIndex << std::endl;

				for (int i=0;i<3;i++){					
					pointB1[i]=min(max(floor(currentIndex[i]-hSize*B[i]+0.5),0),iSize[i]-1);
					pointB2[i]=min(max(floor(currentIndex[i]+hSize*B[i]+0.5),0),iSize[i]-1);
				}

				//std::cout << "cI " << currentIndex << std::endl;
				if (t){
					t=false;
					image->TransformIndexToPhysicalPoint(pointB1, PpointB1);
					image->TransformIndexToPhysicalPoint(pointB2, PpointB2);
					std::cout << "B1 " << convert_point(PpointB1) << pointB1 << std::endl;
					std::cout << "B2 " << convert_point(PpointB2) << pointB2 << std::endl;
				}

				//ImageLineItType itB(resampledImage, pointB1, pointB2);
				ImageLineItType itB(image, pointB1, pointB2);

				CSb=0;
				//std::cout << "B " ;
				//Iterate in directionB
				while (!itB.IsAtEnd())
				{
					cSection[CSa][CSb]=(PixelType)itB.Get();	
					cSectionI[CSa][CSb]=itB.GetIndex();

					//std::cout << cSection[CSa][CSb] << " " ;
					if (maxIntensity<cSection[CSa][CSb]){
						maxIntensity=cSection[CSa][CSb];
						maxIIndex=itB.GetIndex();
					}

					if (minIntensity>cSection[CSa][CSb]){
						minIntensity=cSection[CSa][CSb];
					}
					//std::cout << ". ";						
					++itB;
					CSb++;
				}				
				//std::cout << std::endl;
				++itA;
				//std::cout << "CSa " << CSa << std::endl;
				CSa++;
			}	
			//time
			final1T=clock();

			time1+=(double)(final1T-start1T)/CLOCKS_PER_SEC;				

			CS_center(cSection, CSsize, 3, Center2D);
			
			if (write6){
				write_2D_file("newText.txt", cSection, CSsize[0], CSsize[1], thisLabel);
			}

			image->TransformIndexToPhysicalPoint(maxIIndex, maxPhysicalIIndex);	
			
			/*for (int i=0;i<3;i++){		
				maxIndex[i]=min(max(floor(0.5+pointIndex[i]+
					(float)(Center2D[0]/floor((float)CSsize[0]/2))*hSize*A[i]+
					(float)(Center2D[1]/floor((float)CSsize[1]/2))*hSize*B[i]
				),0),iSize[i]-1);				
			}			*/

			/*for (int i=0;i<3;i++){		
				maxIndex[i]=min(max(pointIndex[i]+
					(float)(Center2D[0]/floor((float)CSsize[0]/2))*hSize*A[i]+
					(float)(Center2D[1]/floor((float)CSsize[1]/2))*hSize*B[i]
				,0),iSize[i]-1);				
			}		*/	

			maxIndex=cSectionI[Center2D[0]+(int)floor((float)CSsize[0]/2)]
			[Center2D[1]+(int)floor((float)CSsize[1]/2)];
		
			image->TransformIndexToPhysicalPoint(maxIndex, maxPhysicalIndex);	
			if (write4){
				std::cout << "Center2D  " << Center2D[0] << " " << Center2D[1] << std::endl;
				std::cout << "pI " << pointIndex << " mI " << maxIndex << std::endl;		
				std::cout << "mII " << maxIIndex << " ; " << setprecision(2) << convert_point(maxPhysicalIIndex) << std::endl;
				std::cout << "su " << physicalIndex << " e " << setprecision(2) << convert_point(maxPhysicalIndex) << std::endl;
				write4=false; //ensures this happens only once
			}

			for (int i=0;i<CSsize[0];i++){
					delete [] cSection[i];
			}
			delete [] cSection;

			return maxPhysicalIndex;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::TwoD_Maximum_fast(PhysicalPointType physicalIndex,  MType E_matrix){

			//time
			clock_t start1T, final1T;

			start1T=clock();

			//********Resample image********
			IndexType pointIndex;
			PhysicalPointType maxPhysicalIndex;		

			float reScale = scale*spacing[0]; //physical radius of vessel

			SizeType   iSize;

			iSize=image->GetLargestPossibleRegion().GetSize();

			image->TransformPhysicalPointToIndex(physicalIndex, pointIndex);

			//std::cout << "iSize   : "  << iSize[0] << " " << iSize[1] << " " << iSize[2] << std::endl;
			//std::cout << "re size   : "  << resampledImage->GetLargestPossibleRegion().GetSize()<< std::endl;

			//********Extract eigen vectors********
			IndexType pointA1, pointB1, pointA2, pointB2, currentIndex, maxIndex, AIndex, BInd1, BInd2;
			IndexType pointTemp;
			PhysicalPointType PpointTemp;
			PhysicalPointType PpointA1, PpointB1, PpointA2, PpointB2;
			float A[3], B[3]; 
			float mag;
			//float hSize=reScale/(reSpacing[0]*sqrt(float(2)));
			float hSize=2*reScale/(reSpacing_[0]);
			//divide by root(2) because that is max sum of two orthogonal vectors

			for (int i=0; i<3;i++){
				A[i]=E_matrix[0][i]; 
			}
			mag=sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2)); 
			for (int i=0; i<3;i++){
				A[i]=A[i]/mag; 
			}

			for (int i=0; i<3;i++){
				B[i]=E_matrix[1][i]; 
			}
			mag=sqrt(pow(B[0],2)+pow(B[1],2)+pow(B[2],2)); 
			for (int i=0; i<3;i++){
				B[i]=B[i]/mag; 
			}
			
			/*if (write2){
				std::cout << "PI : " << convert_point(physicalIndex) << std::endl;	
				std::cout << "I : " << pointIndex << std::endl;	
				std::cout << "A->" << A[0] << "," << A[1] << "," << A[2] << "; B->" << B[0] << "," << B[1] << "," << B[2] << std::endl;
			}*/

			//Iterate in direction A

			for (int i=0;i<3;i++){		
				pointA1[i]=min(max(floor(pointIndex[i]-hSize*A[i]+0.5),0),iSize[i]-1);				
				pointA2[i]=min(max(floor(pointIndex[i]+hSize*A[i]+0.5),0),iSize[i]-1);				
			}

			//if (write2){
			//	image->TransformIndexToPhysicalPoint(pointA1, PpointA1);
			//	image->TransformIndexToPhysicalPoint(pointA2, PpointA2);
			//	std::cout << "hasz " << hSize << std::endl;
			//	std::cout << "A1 " << convert_point(PpointA1) << " " << pointA1 << std::endl;
			//	std::cout << "A2 " << convert_point(PpointA2) << " " << pointA2 << std::endl;
			//}
	
			//ImageLineItType itA(resampledImage, pointA1, pointA2); //iterator in directionA
			ImageLineItType itA(image, pointA1, pointA2); //iterator in directionA
			itA.GoToBegin();
			PixelType currentValue, maxValue=-1;
			currentValue=(PixelType)itA.Get();	
			maxIndex=pointIndex;

			//temp
			bool t=false;
			if (write2){
				t=true;
			}

			while (!itA.IsAtEnd())
			{
				//Determine points of traversing in the B direction
				currentIndex=itA.GetIndex();		

				for (int i=0;i<3;i++){					
					pointB1[i]=min(max(floor(currentIndex[i]-hSize*B[i]+0.5),0),iSize[i]-1);
					pointB2[i]=min(max(floor(currentIndex[i]+hSize*B[i]+0.5),0),iSize[i]-1);
				}

				//std::cout << "cI " << currentIndex << std::endl;
				/*if (t){
					t=false;
					image->TransformIndexToPhysicalPoint(pointB1, PpointB1);
					image->TransformIndexToPhysicalPoint(pointB2, PpointB2);
					std::cout << "B1 " << convert_point(PpointB1) << std::endl;
					std::cout << "B2 " << convert_point(PpointB2) << std::endl;
				}*/

				//ImageLineItType itB(resampledImage, pointB1, pointB2);
				ImageLineItType itB(image, pointB1, pointB2);

				//Iterate in directionB
				while (!itB.IsAtEnd())
				{
					currentValue=(PixelType)itB.Get();
					if (currentValue>maxValue){
						maxValue=currentValue;
						maxIndex=itB.GetIndex();
						AIndex=currentIndex;
						BInd1=pointB1;
						BInd2=pointB2;
					}			

					++itB;
				}
				++itA;
			}	
			//time
			final1T=clock();

			time1+=(double)(final1T-start1T)/CLOCKS_PER_SEC;	
			
			//resampledImage->TransformIndexToPhysicalPoint(AIndex, maxPhysicalIndex);
			image->TransformIndexToPhysicalPoint(AIndex, maxPhysicalIndex);
			//std::cout << "A ind  " << maxPhysicalIndex << std::endl;	

			//resampledImage->TransformIndexToPhysicalPoint(maxIndex, maxPhysicalIndex);	
			image->TransformIndexToPhysicalPoint(maxIndex, maxPhysicalIndex);	
			/*if (write2)
				std::cout << "s " << physicalIndex << " e " << maxPhysicalIndex << std::endl;*/

			return maxPhysicalIndex;
	}

	
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::move_Point(PhysicalPointType originalPhysicalPoint, MType E_matrix){

			PhysicalPointType nextPhysicalPoint;
			ImageType::PixelType Ax, Ay, Az, mag;

			//Extract eigen vectors
			Ax=E_matrix[2][0]; Ay=E_matrix[2][1]; Az=E_matrix[2][2];
			mag=sqrt(pow(Ax,2)+pow(Ay,2)+pow(Az,2)); 
			Ax=Ax/mag; Ay=Ay/mag; Az=Az/mag;
			
			//std::cout << "Direction: " << Ax << "," << Ay << "," << Az << std::endl;
			//std::cout << std::endl;
			nextPhysicalPoint[0]=originalPhysicalPoint[0]+beta*spacing[0]*Ax;
			nextPhysicalPoint[1]=originalPhysicalPoint[1]+beta*spacing[0]*Ay;
			nextPhysicalPoint[2]=originalPhysicalPoint[2]+beta*spacing[0]*Az;

			return nextPhysicalPoint;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::flip_sign(MType E_mtx){

			//Extract eigen vectors
			E_mtx[2][0]=E_mtx[2][0]*(-1); 
			E_mtx[2][1]=E_mtx[2][1]*(-1); 
			E_mtx[2][2]=E_mtx[2][2]*(-1); 	

			return E_mtx;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::change_sign(MType E_mtx_prev, MType E_mtx){

			int sign, k, i;
			float dot_product, dot_product_new, magnitude, max, temp;
			ImageType::PixelType A;
			k=2;
			dot_product=E_mtx_prev[2][0]*E_mtx[2][0]+E_mtx_prev[2][1]*E_mtx[2][1]+E_mtx_prev[2][2]*E_mtx[2][2];

			if (dot_product>1)
				dot_product=1;

			//check for switched evectors
			/*
			if (E_mtx[3][2]<=0){ //if the largest evalue is negative
			max=E_mtx[3][2];		
			for (int i=0;i<3;i++){
			//check for the largest dot_product in vectors 1&2
			dot_product_new=E_mtx_prev[2][0]*E_mtx[i][0]+E_mtx_prev[2][1]*E_mtx[i][1]+E_mtx_prev[2][2]*E_mtx[i][2];
			if (dot_product_new>dot_product){
			dot_product=dot_product_new;
			k=i;
			}
			}
			//swap vectors		
			for (int j=0;j<3;j++){
			temp=E_mtx[2][j];
			E_mtx[2][j]=E_mtx[k][j];
			E_mtx[k][j]=temp;
			}
			//std::cout << "Swapped evectors" << std::endl;
			}
			*/
			//Calculate dot product and sign
			//dot_product=E_mtx[1][0]*E_mtx[2][0]+E_mtx[1][1]*E_mtx[2][1]+E_mtx[1][2]*E_mtx[2][2];

			//dot_product=E_mtx[0][0]*E_mtx[2][0]+E_mtx[0][1]*E_mtx[2][1]+E_mtx[0][2]*E_mtx[2][2];
			//std::cout << "dot product 2 " << dot_product << std::endl;

			//check for Flipped direction evector (3rd one)
			if (dot_product<0){
				E_mtx=flip_sign(E_mtx);
				dot_product=dot_product*-1;
			}

			//dot_product=E_mtx_prev[2][0]*E_mtx[2][0]+E_mtx_prev[2][1]*E_mtx[2][1]+E_mtx_prev[2][2]*E_mtx[2][2];
			//std::cout << "dot product " << dot_product << std::endl;

									
				//average vectors
			double weight=0.5; //0.8

			for (i=2;i<3;i++){
				//calculate magnitude of two vectors
				magnitude= sqrt(pow(E_mtx[i][0]*weight+E_mtx_prev[i][0]*(1-weight),2)
					+pow(E_mtx[i][1]*weight+E_mtx_prev[i][1]*(1-weight),2)
					+pow(E_mtx[i][2]*weight+E_mtx_prev[i][2]*(1-weight),2));

				E_mtx[i][0]=(E_mtx[i][0]*weight+E_mtx_prev[i][0]*(1-weight))/magnitude; 
				E_mtx[i][1]=(E_mtx[i][1]*weight+E_mtx_prev[i][1]*(1-weight))/magnitude; 
				E_mtx[i][2]=(E_mtx[i][2]*weight+E_mtx_prev[i][2]*(1-weight))/magnitude;
			}

			//Ensure the cross-section is orthogonal
			VectorType Evec1, Evec2, Evec3;
			Evec3[0]=E_mtx[2][0];	Evec3[1]=E_mtx[2][1];	Evec3[2]=E_mtx[2][2];
			compute_planeV(Evec3, &Evec1, &Evec2);
			E_mtx[0][0]=Evec1[0];	E_mtx[0][1]=Evec1[1];	E_mtx[0][2]=Evec1[2];
			E_mtx[1][0]=Evec2[0];	E_mtx[1][1]=Evec2[1];	E_mtx[1][2]=Evec2[2];			

			
			if (write2){
				std::cout<< "dotProd " << dot_product << std::endl;
				std::cout<< "--New Evecs--" << std::endl;
				std::cout<<"["<<E_mtx[0][0]<<","<<E_mtx[0][1]<<","<<E_mtx[0][2]<<"];"<<std::endl;
				std::cout<<"["<<E_mtx[1][0]<<","<<E_mtx[1][1]<<","<<E_mtx[1][2]<<"];"<<std::endl;
				std::cout<<"["<<E_mtx[2][0]<<","<<E_mtx[2][1]<<","<<E_mtx[2][2]<<"];"<<std::endl;				
			}

			if (write2){
				std::cout<< "prev Evecs" << std::endl;
				std::cout<<"["<<E_mtx_prev[0][0]<<","<<E_mtx_prev[0][1]<<","<<E_mtx_prev[0][2]<<"];"<<std::endl;
				std::cout<<"["<<E_mtx_prev[1][0]<<","<<E_mtx_prev[1][1]<<","<<E_mtx_prev[1][2]<<"];"<<std::endl;
				std::cout<<"["<<E_mtx_prev[2][0]<<","<<E_mtx_prev[2][1]<<","<<E_mtx_prev[2][2]<<"];"<<std::endl;				
			}

			//dot_product=E_mtx_prev[2][0]*E_mtx[2][0]+E_mtx_prev[2][1]*E_mtx[2][1]+E_mtx_prev[2][2]*E_mtx[2][2];

			//if (write3){
			//	std::cout<< "dotProd " << dot_product << std::endl;		
			//}

			return E_mtx;
	}

	//assume correct orientation wrt direction_
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::correct_orientation(MType E_mtx, VectorType direction_){ 	

			float dot_product;
			dot_product=direction_[0]*E_mtx[2][0]+direction_[1]*E_mtx[2][1]+direction_[2]*E_mtx[2][2];

			if (dot_product<0){
				E_mtx=flip_sign(E_mtx);	
				//std::cout << "orientation corrected ******" << std::endl;	
			}

			return E_mtx;
	}
/*
		//assume correct orientation wrt direction_
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::correct_singularity(MType E_mtx){ 	
			
			float thresHold=0.98;
			float dotProduct1,dotProduct2;
			dotProduct1=abs( E_mtx[1][0]*E_mtx[2][0]+E_mtx[1][1]*E_mtx[2][1]+E_mtx[1][2]*E_mtx[2][2] );
			dotProduct2=abs( E_mtx[0][0]*E_mtx[2][0]+E_mtx[0][1]*E_mtx[2][1]+E_mtx[0][2]*E_mtx[2][2] );

			if (dotProduct1 >thresHold && dotProduct2 >thresHold){

			}

			//compute_planeV (VectorType C, VectorType* A, VectorType* B);
			
			dot_product=direction_[0]*E_mtx[2][0]+direction_[1]*E_mtx[2][1]+direction_[2]*E_mtx[2][2];

			if (dot_product<0){
				E_mtx=flip_sign(E_mtx);	
				//std::cout << "orientation corrected ******" << std::endl;	
			}

			return E_mtx;
	}
*/
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::correct_orientation(PhysicalPointType curPoint, PhysicalPointType refPoint, MType E_mtx){ 	

			int i;
			VectorType refTOcurrent;

			for (i=0;i<3;i++)
				refTOcurrent[i]=curPoint[i]-refPoint[i]; 

			return correct_orientation(E_mtx, refTOcurrent);
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::print_Emtx(MType E_mtx){ 	
			
		std::cout<< "---" << std::endl;
		std::cout<< "VD " << E_mtx[2][0] << " " << E_mtx[2][1] << " " << E_mtx[2][2] << std::endl;
		std::cout<< "VC " << E_mtx[0][0] << " " << E_mtx[0][1] << " " << E_mtx[0][2] << std::endl;
		std::cout<< "VC " << E_mtx[1][0] << " " << E_mtx[1][1] << " " << E_mtx[1][2] << std::endl;
		std::cout<< "---" << std::endl;
				
	}

	//seedDirection is both input and output
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::refine_seed(PhysicalPointType unrefinedSeed, PhysicalPointType* refinedSeed, VectorType* seedDirection, PhysicalPointType* testResult,MType E_mtx1){	
			int i=0;
			MType E_mtx,E_mtx_prev; 
			VectorType rDirection;
			float tLength; //"this" length that has been computed thus far
			PhysicalPointType pPoint, pEst;
			float mag=0;
			/*
			//if seed is in the middle of nowhere
			if (OneD_Maximum(unrefinedSeed,refinedSeed)==0){
			return(-1);
			}  	
			//Get eigenvectors  
			E_mtx=compute_EigenMatrix(HessianFilter->GetOutput(), (*refinedSeed));  
			E_mtx=correct_orientation(E_mtx, *seedDirection); //with respect to direction
			//Get 2-D Max  
			*/

			if (fastSwitch){
				*refinedSeed=TwoD_Maximum_fast(unrefinedSeed, E_mtx1); 
			} else{
				*refinedSeed=TwoD_Maximum(*refinedSeed, E_mtx1); 
			}
			
			//Set-up
			E_mtx_prev=compute_EigenMatrix(HessianFilter->GetOutput(), (*refinedSeed));
			E_mtx_prev=correct_orientation(E_mtx_prev, *seedDirection); //with respect to direction

			pPoint=*refinedSeed;	

			//**Estimate next point**
			for (i=0;i<9;i++){
				E_mtx=compute_EigenMatrix(HessianFilter->GetOutput(), pPoint);  
				E_mtx=change_sign(E_mtx_prev,E_mtx);
				E_mtx=correct_orientation(E_mtx, *seedDirection); //with respect to direction
				pEst=move_Point(pPoint, E_mtx);
				
				if (fastSwitch){
					pPoint=TwoD_Maximum_fast(pEst, E_mtx);
				} else {
					pPoint=TwoD_Maximum(pEst, E_mtx);
				}
				//
				E_mtx_prev=E_mtx;
				tLength=compute_distance(pPoint,(*refinedSeed));	
				if(tLength>=directionLength){
						break;
					}
			} //shouldn't this end at the direction length? Why should we compare 
			// how one vessel be 8 mm down the seed to how that vessel will be only 3 mm down

			(*testResult)=pPoint;

			for (i=0;i<3;i++){
				rDirection[i]=pPoint[i]-(*refinedSeed)[i];	
				mag+=pow(rDirection[i],2);
			}

			for (i=0;i<3;i++)
				rDirection[i]/=sqrt(mag);

			(*seedDirection)=rDirection;

			return sqrt(mag);
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::optimize_seed(PhysicalPointType seedPoint, VectorType v_dir, PhysicalPointType expectedSeed, PhysicalPointType* optimalPoint){
			int i,j, output, index=0;
			float variate1,variate2,variate3; //three degrees of freedom in Cyclindrical
			float var2, var3; //two degrees of freedom converted back to Cartesian
			float dist, score, dotProd, maxDotProd=0, maxScoreDist=30, maxScore=0;
			float radialBase=0, axialIncrease=0;	
			MType E_mtx,E_mtx1; 
			PhysicalPointType seedCandidate, bestSeed, temp, bestTemp, testR, bestTest;	
			VectorType directions, maxDirections; //predicts the direction each seed will go to
			PhysicalPointType displacement;
			VectorType vectorC1, vectorC2;	//v_dir->vessel direction; C1,C2 ->cross-section
			bool endFunction=true;
			float mag;

			//Compute the cross-section vectors
			compute_planeV(v_dir, &vectorC1, &vectorC2);
			E_mtx1=make_EigenMatrix(v_dir,vectorC1,vectorC2);

			/*
			std::cout << v_dir << std::endl;
			std::cout << vectorC1 << std::endl;
			std::cout << vectorC2 << std::endl;*/

			//Compute seed estimates (with noise)
			for (i=0;i<maxSeeds;i++){
				variate1= randomGenerator->GetUniformVariate(-1*axialSpread,axialSpread+axialIncrease);				
				//in POLAR coords
				variate2= randomGenerator->GetUniformVariate(0,2); //get angle in terms of PI
				variate3= randomGenerator->GetUniformVariate(radialBase,radialSpread+radialBase); //get radius
				//convert to Cartesian
				var2=cos(variate2*PI)*variate3;
				var3=sin(variate2*PI)*variate3;
				for (j=0;j<3;j++)
					seedCandidate[j]=expectedSeed[j]+v_dir[j]*variate1*spacing[j]+vectorC1[j]*var2*spacing[j]+vectorC2[j]*var3*spacing[j];	
				directions=v_dir;
				temp=seedCandidate;

				mag=refine_seed(seedCandidate,&(seedCandidate),&(directions),&(testR),E_mtx1);
				if (mag!=-1){			
					E_mtx=compute_EigenMatrix(HessianFilter->GetOutput(), seedCandidate);
					score=0;
					dist=0;
					dotProd=0;
					for (j=0;j<3;j++){
						dist+=pow(expectedSeed[j]-seedCandidate[j],2);
						dotProd+=v_dir[j]*directions[j];
					}
					score=dotProd/(sqrt(sqrt(dist)));			
					//			std::cout << mag << "  " << score << " " << std::endl;
					if (score>maxScore){ 	
						endFunction=false;	/*//only continue if score is high enough*/
						bestSeed=seedCandidate;
						bestTemp=temp;
						bestTest=testR;
						maxScore=score;
						maxDotProd=dotProd;
						maxScoreDist=dist;
						maxDirections=directions;				
					}	
					if (maxDotProd>angThresh && maxScoreDist<distThresh ){

						break;
					}
				}
				radialBase+=radialIncrement;
				axialIncrease+=axialIncrement;
			}
			
			//std::cout << "sc: " << maxScore << " dst: " << maxScoreDist << std::endl;	
			//std::cout << " BSd: " << bestSeed << " BDir: " << maxDirections << std::endl;	

			if (write){
				std::cout << "sc: " << maxScore << " BTest: " << bestTest << std::endl;
				std::cout << " SeedP: " << seedPoint  << "   Score: " << maxScore << std::endl;
				std::cout << " ExpSeed: " << expectedSeed << "  VD: " << v_dir << " VDP: " << compute_point(expectedSeed,v_dir) << std::endl;
				std::cout << " BSeed: " << bestSeed << " BDir: " << maxDirections << " BDirP: " << compute_point(bestSeed,maxDirections) << std::endl;
			}
			
			*optimalPoint=bestSeed;
			if (endFunction == true)
				return -1;
			else
				return 1;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::optimize_seed(PhysicalPointType seedPoint, VectorType v_dir, PhysicalPointType* optimalPoint){
			int returnValue;
			PhysicalPointType expectedSeed;
			float axialDist=axialMean;
			//Best estimate for new seed point	
			for (int j=0;j<3;j++)	
				expectedSeed[j]=seedPoint[j]+v_dir[j]*axialDist*spacing[j];

			returnValue=optimize_seed(seedPoint,v_dir,expectedSeed,optimalPoint);
			return returnValue;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_v_dir(int n){
			int j;
			float mag=0;
			for (j=0;j<3;j++){
				infoVector[n].v_direction[j]=infoVector[n].direction[j]-infoVector[n].seedU[j];
				mag+=pow(infoVector[n].v_direction[j],2);		
			}

			for (j=0;j<3;j++){
				infoVector[n].v_direction[j]/=sqrt(mag);					
			}					
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::initialize_LOOP(){

			//tracking distances
			for (int i=0; i<20;i++)
				distances[i]=0;	

			//tracking points
			trackable=3;
			for (int i=0;i<2;i++){
				prevPoint[i]=0;
				prevVector[i]=0;
				prevVector2[i]=0;
			}

	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::get_point(vtkSmartPointer<vtkIdList> pLine, int pInd, PhysicalPointType* returnPoint){
			//PhysicalPointType returnPoint;
			int i;
			int nPtsInLine;
			double gotPoint[3];
			
			//checkLines->GetCell(loc,gotLine); //loc,nPtsInLine,
			if (pLine->GetNumberOfIds()<=pInd){
				//std::cout << "Lim " << gotLine->GetNumberOfIds() << " PI " << pInd << std::endl;
				return -1;
			}

			vtkIdType pointId=pLine->GetId(pInd);	
			//std::cout << "PID " << pointId << std::endl;
			checkPts->GetPoint(pointId, gotPoint);

			(*returnPoint)[0]=-1*gotPoint[0];
			(*returnPoint)[1]=-1*gotPoint[1];
			(*returnPoint)[2]=gotPoint[2];

			//return returnPoint;
			return 1;
	}

	//INEFFICIENT
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::get_point(int pLine, int pInd, PhysicalPointType* returnPoint){
			//PhysicalPointType returnPoint;
			int i, nPtsInLine, returnValue;
			double gotPoint[3];
			vtkSmartPointer<vtkIdList> gotLine = vtkSmartPointer<vtkIdList>::New();
			vtkIdType loc=pLine;

			checkLines->InitTraversal(); 
			for(i=0;i<=pLine;i++){
				checkLines->GetNextCell(gotLine);	
			}

			//checkLines->GetCell(loc,gotLine); //loc,nPtsInLine,
			if (gotLine->GetNumberOfIds()<=pInd){
				//std::cout << "Lim " << gotLine->GetNumberOfIds() << " PI " << pInd << std::endl;
				pInd=gotLine->GetNumberOfIds()-1; //fetch the last point
				returnValue =-1;
			} else {
				returnValue = 1;
			}

			vtkIdType pointId=gotLine->GetId(pInd);	
			//std::cout << "PID " << pointId << std::endl;
			checkPts->GetPoint(pointId, gotPoint);

			(*returnPoint)[0]=-1*gotPoint[0];
			(*returnPoint)[1]=-1*gotPoint[1];
			(*returnPoint)[2]=gotPoint[2];

			//return returnPoint;
			return returnValue;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::start_checks(int n){	
			PhysicalPointType pointA=infoVector[n].seed;
			//PhysicalPointType pointB=infoVector[n].newSeed;
			PhysicalPointType pointB=infoVector[n].point2;
			
			int i=0;
			float cosAng=-1;
			int gotPoint=1;
			checkInd=infoVector[n].seedLineInd-3; 
			//-3 because of MATLAB-C++ indexing differences
			//checkInd++; //to start off with the next point
			checkLin=infoVector[n].seedLine-1;
			//std::cout << "N " << n << std::endl;	
			//	std::cout << "L " << checkLin << " I " << checkInd << std::endl;	
			//Get line
			PhysicalPointType checkPoint;	
			gotPoint=get_point(checkLin,checkInd,&checkPoint);

			while (gotPoint==1){				
				cosAng=compute_angle(pointA,checkPoint,pointB);		
				//std::cout << "P "<< checkPoint << " Ang " << cosAng << std::endl;
				if (cosAng>=0){
					break;
				}
				prevCheck=checkPoint; //prevCheck is one below
				checkInd++;
				gotPoint=get_point(checkLin,checkInd,&checkPoint);
				//i++;
			}

			newStart=checkPoint;			
			distRemainder=0;
			prevTracked=pointB;
			prevTrackedI=2; //this doesn't seem to be used much

			//initialize
			nextSeed=prevCheck;
			nextDirection=checkPoint;
			
			/*
			std::cout << "S " <<pointA << std::endl;
			std::cout << "St " <<checkPoint << std::endl;	// Later on C1
			std::cout << "NS " <<pointB << std::endl; // Later on T1
			std::cout << "CI " <<checkInd << std::endl;
			std::cout << "StartDone " << std::endl;*/
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::PhysicalPointType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::new_checkPoint(float dist){
			//std::cout << "D " << dist;
			dist-=distRemainder;
			float refDist=0;
			int gotPoint=0;
			PhysicalPointType checkPoint;
			//checkInd++; //raise by at least one
			//std::cout << " D " << dist << " rD1 " << refDist << " CI1 " << checkInd ;

			do{		
				gotPoint=get_point(checkLin,checkInd,&checkPoint);
				if (gotPoint!=1){			
					break;
				}
				refDist=compute_distance(prevCheck,checkPoint);
				checkInd++;				
			} while(refDist<=dist);

			distRemainder=refDist-dist;
			//std::cout << " CI2 " << checkInd  << " rD2 " <<  refDist << " dRem " <<distRemainder << std::endl;

			return checkPoint;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::check_orientation(PhysicalPointType currTracked, float dist){

			PhysicalPointType checkPoint=new_checkPoint(dist);
			//determine new checkPoint

			//check orientation
			float cosAngle=compute_angle(checkPoint,prevCheck,prevTracked,currTracked);	

			nextSeed=prevCheck;
			nextDirection=checkPoint;
			if (cosAngle>-2){
				prevCheck=checkPoint;	
				if (cosAngle<angleThresh){
					/*
					std::cout << "cA----------" << cosAngle << std::endl;

					std::cout << "C1, " <<convert_point(prevCheck) << ",1,1" << std::endl;
					std::cout << "C2, " <<convert_point(checkPoint) << ",1,1" << std::endl;
					std::cout << "T1, " <<convert_point(prevTracked) << ",1,1" << std::endl;	
					std::cout << "T2, " <<convert_point(currTracked) << ",1,1" << std::endl;
					std::cout << "TDist, " <<compute_distance(prevTracked,currTracked)
						<< "  CDist " << compute_distance(prevCheck,checkPoint) << std::endl;*/				
					return -1;
				}else{
					/*
					std::cout << "cA     " << cosAngle << std::endl;
					std::cout << "C1, " <<convert_point(prevCheck) << ",1,1" << std::endl;
					std::cout << "C2, " <<convert_point(checkPoint) << ",1,1" << std::endl;
					std::cout << "T1, " <<convert_point(prevTracked) << ",1,1" << std::endl;	
					std::cout << "T2, " <<convert_point(currTracked) << ",1,1" << std::endl;	
					
					std::cout << "TDist, " <<compute_distance(prevTracked,currTracked)
						<< "  CDist " << compute_distance(prevCheck,checkPoint) << std::endl;				
				*/
					return 1;
				}
			}else{
				return 0;
			}
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::check_check(PhysicalPointType currTracked, int currI){
			int returnValue;
			float dist=compute_distance(currTracked,prevTracked);
			//std::cout << "CT " << currTracked << std::endl;
			//std::cout << "PT " << prevTracked << std::endl;

			if (dist>distanceThresh){
				returnValue=check_orientation(currTracked,dist);
				if(returnValue>0){
					//std::cout << "1" << std::endl;
					prevTracked=currTracked;			
				} 
				prevTrackedI=currI;
				//std::cout << "-1" << std::endl;
				//std::cout << "EP " << currTracked << std::endl;
				//std::cout << "ECP " << prevCheck << std::endl;
				//prevTracked=currTracked;  Be careful when uncommenting this
				return returnValue;		
			} else {
				//	std::cout << "0" << std::endl;
				return 2;	
			}
			//2 -> not enough information
			//1  -> everything ok 
			//0  -> something wrong with angle calculation
			//-1 -> tracking deviates significantly from guidance

	}

	
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::check_line_orientiation(int n, vtkPoints *inPoints, vtkPoints *outPoints, bool termP){
			int i,j;
			int returnValue=0;
			int decide;
			double getP[3];	
			PhysicalPointType thePoint;
			int total_points=(int)inPoints->GetNumberOfPoints();

			if (n==61){
				i =0;
			}

			/*
			//initialize new_seed
			inPoints->GetPoint(1,getP);	
			thePoint[0]=getP[0]*(-1);
			thePoint[1]=getP[1]*(-1);
			thePoint[2]=getP[2];
			infoVector[n].newSeed=thePoint;*/
			
			if (total_points>2){
				inPoints->GetPoint(2,getP);	
			} else {
				inPoints->GetPoint(1,getP);	
			}
			thePoint[0]=getP[0]*(-1);
			thePoint[1]=getP[1]*(-1);
			thePoint[2]=getP[2];
			infoVector[n].point2=thePoint;
			
			start_checks(n);

			for (i=3; i<total_points;i++){			
				inPoints->GetPoint(i,getP);	
				thePoint[0]=getP[0]*(-1);
				thePoint[1]=getP[1]*(-1);
				thePoint[2]=getP[2];
				//std::cout<< "thePoint " << thePoint << std::endl;

				decide=check_check(thePoint,i);				

				if (decide==0 && termP==true){
					decide=-1;
				}
				if (decide==-1){
					total_points=prevTrackedI; 
					std::cout << " Ptracked " << prevTrackedI << std::endl;
					
					returnValue=1;
					break;
				} else if (decide==0){
					break;
				}
			}

			if (returnValue==1 || termP==true){
					infoVector[n].seedU=nextSeed;
					infoVector[n].direction=nextDirection;
					compute_v_dir(n);


					/*std::cout << " decide " << decide << std::endl;
					std::cout << " nextSeed " << convert_point(nextSeed) << std::endl;
					std::cout << " nextDirection " << convert_point(nextDirection) << std::endl;*/
			}

			//reorder the first points
			for (i=0;i<total_points;i++){

				inPoints->GetPoint(i,getP);
				outPoints->InsertNextPoint (getP[0], getP[1],getP[2]);
			}	

			return returnValue;
	}

	//Default: not terminated by circles
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::check_line_orientiation(int n, vtkPoints *inPoints, vtkPoints *outPoints){

			return check_line_orientiation(n,inPoints,outPoints,false);

	}

		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::track_distances(PhysicalPointType curPoint, PhysicalPointType refPoint){
			float distance_;
			distance_=compute_distance(curPoint,refPoint);
			//distance_=sqrt(pow(refPoint[0]-curPoint[0],2)+pow(refPoint[1]-curPoint[1],2)+pow(refPoint[2]-curPoint[2],2));

			//store last 20 distances
			for (int i=19; i>0;i--){
				distances[i]=distances[i-1];
			}

			distances[0]=distance_;

			return distance_;
	}

			template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::track_directions(PhysicalPointType curPoint){
			
			prevVector2=prevVector;
			prevVector=compute_vector(curPoint, prevPoint);
			
			prevPoint=curPoint;

			if (trackable>0)
					trackable--;

			//std::cout << "trackable " << trackable << std::endl;
			
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	bool  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::test_NegEvalue(MType E_matrix){	
			bool result;

			if (E_matrix[3][0]<0 && E_matrix[3][1]<0)
				result = false;
			else //terminate
				result = true;

			return result;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float   
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::test_GradientZero(IndexType point, MType E_matrix){
			float pVx, pVy, pVz;
			float dotProduct1, dotProduct2, mag;

			//Calculate gradient
			MType G=calculate_Gradient(point);
			pVx=G[0][0]*1.0; pVy=G[1][0]*1.0; pVz=G[2][0]*1.0;

			//dot product squared
			dotProduct1=pow((pVx*E_matrix[0][0]+pVy*E_matrix[0][1]+pVz*E_matrix[0][2]),2); 
			dotProduct2=pow((pVx*E_matrix[1][0]+pVy*E_matrix[1][1]+pVz*E_matrix[1][2]),2); 

			//std::cout << "T " << dotProduct1+dotProduct2<< std::endl;

			return dotProduct1+dotProduct2;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float   
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::RidgeSymmetry(MType E_matrix){	
			float symmetry= E_matrix[3][1]/E_matrix[3][0] - 1 ;

			if (symmetry<0)
				symmetry*=-1;
			//std::cout << "Symmetry: " << symmetry  << std::endl;

			return symmetry;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	bool  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::termination(PhysicalPointType physicalPoint, MType E_matrix, int n){
			IndexType point;
			float backCount=0;
			float dotProduct=0;
			bool terminate_=false;
			float distance__;
			int i;
			
			//do some basic calculations
			image->TransformPhysicalPointToIndex(physicalPoint, point);
			track_distances(physicalPoint,infoVector[n].seedU); //or new_seed
			track_directions(physicalPoint);
			distance__=compute_distance(physicalPoint,infoVector[n].seed);

			//.............End now, repeat not.............//

			//terminate by darkness
			float value = image->GetPixel( point );
			
			if(value< intensityThreshold){
					//this is to simply end for this vessel, without retry
					terminate_=true;
					std::cout << "Int " << value << std::endl;
					return (terminate_);
			}

			//termination by length
			if(distance__>infoVector[n].length){
					//this is to simply end for this vessel, without retry
					terminate_=true;
					//std::cout << "L " << distance__ << "  " << infoVector[n].length << std::endl;
					return (terminate_);
			}

			/*
			//termination by gradient
			if ( test_GradientZero(point, E_matrix)>= gradTolerance){
			terminate_=true;
			std::cout << "Termination by improper gradient" << std::endl;
			}

			//termination by positive eigenvalues
			if (test_NegEvalue(E_matrix)){
			terminate_=true;
			std::cout << "Termination by positive eigenvalues" << std::endl;
			}
			*/

			//std::cout << test_NegEvalue(E_matrix) << "  T_Grad: " << test_GradientZero(point, E_matrix) << std::endl;

			//termination by assymetry	
			//if (RidgeSymmetry(E_matrix) >=epsilon){
			//terminate_=true;
			//std::cout << "Termination by assymetry" << std::endl;
			//}

			//.............End now, repeat again.............//

			//termination by going in circles
			for (i=1; i<20;i++){
				if (distances[i]>distances[i-1]) //these are distances from seed
					backCount++;
			}

			if (backCount>8){
				terminate_=true;
				termPremature=true;		
				std::cout << "Term via Circles" << std::endl;
			}
		
			//termination by reversal
			for (i=0;i<3;i++){
				dotProduct+=prevVector[i]*prevVector2[i];
			}
			if (trackable==0){ //check if enough values available
				if (dotProduct<0) {
					//terminate_=true;

					//std::cout << " Term via reversal " << dotProduct<< std::endl;					
				}			
				
			}
			
			//Return conclusion
			return terminate_;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_radius(PhysicalPointType physicalIndex, double prev_rad, double est_rad, MType E_mtx){	

			//declarations	
			//temporary index at which the gradient will be tested
			PhysicalPointType testIndex;
			float radius, final_radius=0;
			float incremen = est_rad/20; // increment
			//coefficients to traverse circle, based on sin(theta)
			float Script1[8]={	1, 0.707, 0, -0.707, 
				-1,-0.707, 0, 0.707 };
			float Script2[8]={	0, 0.707, 1, 0.707, 
				0,-0.707, -1, -0.707 };
			//set the vectors that form the cross-section of the vessel
			MType T(4,3);
			T[0][0]=E_mtx[0][0]; T[0][1]=E_mtx[0][1]; T[0][2]=E_mtx[0][2];
			T[1][0]=E_mtx[1][0]; T[1][1]=E_mtx[1][1]; T[1][2]=E_mtx[1][2];
			float max_intens=0, curr_intens;	

			//***Main Loop***
			for (radius=min(0.75*prev_rad,0.5*est_rad); radius<min(1.25*prev_rad,2*est_rad);radius=radius+incremen){
				curr_intens=0;
				for (int i=0;i<8;i++){// 8 points around

					for (int j=0;j<3;j++) // get the x, y, z of the point to be tested
						testIndex[j]=physicalIndex[j]+radius*(T[0][j]*Script1[i]+T[1][j]*Script2[i]);			

					curr_intens=curr_intens+interpolator2->Evaluate(testIndex);
				}				
				//if the function is maximized
				if (curr_intens>max_intens){
					max_intens=curr_intens;		
					// set the radius with max grad intensity 				
					final_radius=radius*1.41421; //multiply by sqrt(2) because of paper
				}			
			}
			return final_radius;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::LOOPER(){			
			//first_points=MAIN_LOOP(physicalPoint, E_mtx, centerPoints,  radiusArray, radius_prev);
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::reorder_points(vtkPoints *centerPoints, vtkPoints *orderedPoints, int firstpoints){
			double getP[3];	
			int total_points=(int)centerPoints->GetNumberOfPoints();

			//reorder the first points
			for (int i=0; i<firstpoints;i++){			
				centerPoints->GetPoint(firstpoints-1-i,getP);		
				orderedPoints->InsertNextPoint (getP[0], getP[1],getP[2]);
			}
			//the rest are same
			for (int i=firstpoints;i<total_points;i++){
				centerPoints->GetPoint(i,getP);		
				orderedPoints->InsertNextPoint (getP[0], getP[1],getP[2]);
			}

	}

	//add inPoints to the end of allPoints
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_points(vtkPoints *inPoints, vtkPoints *allPoints, int finalInPoint){
			double getP[3];				

			//reorder the first points
			for (int i=0; i<finalInPoint;i++){			
				inPoints->GetPoint(i,getP);	
				allPoints->InsertNextPoint (getP[0], getP[1],getP[2]);
			}

	}

		//add inPoints to the end of allPoints
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_array_elements(vtkFloatArray *inPoints, vtkFloatArray *allPoints, int finalInPoint){
			float getV;				

			//reorder the first points
			for (int i=0; i<finalInPoint;i++){			
				getV=inPoints->GetValue(i);	
				allPoints->InsertNextValue(getV);
			}

	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_points(vtkPoints *inPoints, vtkPoints *allPoints){

			int finalInPoint=(int)inPoints->GetNumberOfPoints();
			write_points(inPoints, allPoints, finalInPoint);
	}
	
		//trim A in accordance with B
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::trim_equal(vtkPoints *pointsA, vtkPoints *pointsAPermanent, int numPoints){
			
			int pointsLengthA=(int)pointsA->GetNumberOfPoints();
			int newPointsLength= min(numPoints,pointsLengthA);
			
			write_points(pointsA, pointsAPermanent, newPointsLength);	
			
			//std::cout << "leng " << newPointsLength << " " << pointsA->GetNumberOfPoints()<< std::endl;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::trim_equal_array(vtkFloatArray *arrayA, vtkFloatArray *arrayAPermanent, int numPoints){
			
			double getP[3];	
			int pointsLengthA=(int)arrayA->GetNumberOfTuples();;
			int newPointsLength= min(numPoints,pointsLengthA);
			
			write_array_elements(arrayA, arrayAPermanent, newPointsLength);
			
			//std::cout << "leng " << newPointsLength << " " << pointsA->GetNumberOfPoints()<< std::endl;
	}


	//Should make smoothing function more accurate by taking into account the distance between points
	//Currently, the computation assumes a distance =1
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::smooth_points(vtkPoints *centerPoints, vtkPoints *fixedPoints, int finalPoint){
			double a[3], getP[3], av[3];
			double weights[11]={1,0.956,0.835,0.667,0.487,0.325,0.198,0.110,0.0563,0.0262,0.0112}; //Gaussian weighted
			//c=3.335
			double weight_sum;	

			for (int i=0; i<finalPoint; i++){	
				weight_sum=0;
				a[0]=a[1]=a[2]=0;
				//Look at the point itself
				centerPoints->GetPoint(i,getP);
				for (int k=0;k<3;k++)
					a[k]+=weights[0]*getP[k];
				weight_sum+=weights[0];
				//Look at its neighbours
				for (int j=-1*windowSize; j<0;j++){
					if((i+j)>=0 && (i-j)<finalPoint){
						centerPoints->GetPoint(i+j,getP);
						for (int k=0;k<3;k++)
							a[k]+=weights[abs(j)]*getP[k];
						weight_sum+=weights[abs(j)];

						centerPoints->GetPoint(i-j,getP);
						for (int k=0;k<3;k++)
							a[k]+=weights[abs(j)]*getP[k];
						weight_sum+=weights[abs(j)];

					}
				}
				//the new points
				fixedPoints->InsertNextPoint (a[0]/weight_sum, a[1]/weight_sum,a[2]/weight_sum);
			}

	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::smooth_points(vtkPoints *centerPoints, vtkPoints *fixedPoints){

			int finalPoint=(int)centerPoints->GetNumberOfPoints();
			smooth_points(centerPoints,fixedPoints,finalPoint);
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::bad_points(vtkPoints *centerPoints){
			int i,j;
			int numPoints=(int)centerPoints->GetNumberOfPoints();
			double getP[3];

			for (i=0;i<numPoints;i++){
				centerPoints->GetPoint(i,getP);
				for (j=0;j<3;j++){
					if (getP[j]>5000 || getP[j]<-5000){
						centerPoints->SetPoint(i,0,0,0);
						break;
					}
				}
			}

	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_pointsVTK(const char *filename, vtkPoints *points, vtkCellArray *centerLine, vtkFloatArray* radiusArray){

			vtkPolyData *polydata = vtkPolyData::New();
			polydata->SetPoints(points);
			polydata->SetLines(centerLine);
			polydata->GetPointData()->AddArray(radiusArray);

			vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
			writer->SetFileName(filename);
			writer->SetInput(polydata);
			writer->Write();
	}

		template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_1D_file(char *filename, vtkFloatArray *points, int length, std::string label){

			ofstream fileToWrite(filename, std::ios_base::app);
			if (fileToWrite.is_open())
			{
				fileToWrite<< label << " = [";
				for (int i=0; i<length;i++){					
					fileToWrite<< points->GetValue(i) << " ";					
					fileToWrite << "\n";
				}
				fileToWrite<< "]; \n";
				fileToWrite.close();
			} else cout << "Unable to open file";
			
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_1D_file(char *filename, float *points, int length, std::string label){

			ofstream fileToWrite(filename, std::ios_base::app);
			if (fileToWrite.is_open())
			{
				fileToWrite<< label << " = [";
				for (int i=0; i<length;i++){					
					fileToWrite<< points[i] << " ";					
					fileToWrite << "\n";
				}
				fileToWrite<< "]; \n";
				fileToWrite.close();
			} else cout << "Unable to open file";
			
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_2D_file(char *filename, float **points, int hieght, int width, std::string label){

			ofstream fileToWrite(filename, std::ios_base::app);
			if (fileToWrite.is_open())
			{
				fileToWrite<< label << " = [";
				for (int i=0; i<hieght;i++){
					for (int j=0; j< width; j++){
						fileToWrite<< points[i][j] << " ";
					}
					fileToWrite << "\n";
				}
				fileToWrite<< "]; \n";
				fileToWrite.close();
			} else cout << "Unable to open file";
			
	}


} // end namespace itk
#endif



