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
		beta= 0.4;
		sigma =2.5;
		sFactor=0.2;

		//other parameters
		max_iterations=140;
		epsilon =0.8;
		gradTolerance=0.3;
		windowSize =5;

		//seed model parameters
		maxSeeds=50;
		angThresh= 0.75;
		distThresh= 3.0;
		axialMean=4.0;	//in terms of spacing
		axialSpread=2.5; //in terms of spacing
		radialSpread=1.5;//in terms of spacing
		radialIncrement=0.02; //the circle will move out
		//every iteration
		axialIncrement=0.02; //

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

		//time
		time1=0;
		time2=0;
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
			clock_t start1T, start2T, final1T, final2T;

			start1T=clock();

			//********Resample image********
			IndexType pointIndex;
			PhysicalPointType maxPhysicalIndex;		

			float reScale = scale*spacing[0]; //physical radius of vessel
			origin[0] = physicalIndex[0]-reScale-reSpacing[0]; //reSpacing[0]->extra padding
			origin[1] = physicalIndex[1]-reScale-reSpacing[0];
			origin[2] = physicalIndex[2]-reScale-reSpacing[0];
			filter->SetOutputOrigin( origin );	 

			ImageType::SizeType   rSize;
			rSize[0]=reScale*2/reSpacing[0]+2; //re_scale/reSpacing[0]=scale in voxels
			rSize[1]=reScale*2/reSpacing[0]+2; //+2 to deal with rounding errors
			rSize[2]=reScale*2/reSpacing[0]+2;
			filter->SetSize( rSize );  
			filter->Update();  
			resampledImage=filter->GetOutput();
			resampledImage->TransformPhysicalPointToIndex(physicalIndex, pointIndex);	

			final1T=clock();

			start2T=clock();	
			//std::cout << "rSize   : "  << rSize<< std::endl;
			//std::cout << "re size   : "  << resampledImage->GetLargestPossibleRegion().GetSize()<< std::endl;

			//********Extract eigen vectors********
			IndexType pointA1, pointB1, pointA2, pointB2, currentIndex, maxIndex, AIndex, BInd1, BInd2;
			PhysicalPointType PpointA1, PpointB1, PpointA2, PpointB2;
			float A[3], B[3]; 
			float mag;
			float hSize=reScale/(reSpacing[0]*sqrt(float(2)));
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

			//std::cout << "A->" << A[0] << "," << A[1] << "," << A[2] << "; B->" << B[0] << "," << B[1] << "," << B[2] << std::endl;

			//Iterate in direction A


			for (int i=0;i<3;i++){
				pointA1[i]=min(max(floor(pointIndex[i]-hSize*A[i]+0.5),0),rSize[i]-1); 
				pointA2[i]=min(max(floor(pointIndex[i]+hSize*A[i]+0.5),0),rSize[i]-1);
			}

			//resampledImage->TransformIndexToPhysicalPoint(pointA1, PpointA1);
			//resampledImage->TransformIndexToPhysicalPoint(pointA2, PpointA2);
			//std::cout << "hasz " << h_size << std::endl;
			//std::cout << "A1 " << convert_point(PpointA1) << std::endl;
			//std::cout << "A2 " << convert_point(PpointA2) << std::endl;

			ImageLineItType itA(resampledImage, pointA1, pointA2); //iterator in directionA
			itA.GoToBegin();
			PixelType currentValue, maxValue=-1;
			currentValue=(PixelType)itA.Get();	
			maxIndex=pointIndex;

			//temp
			bool t=true;

			while (!itA.IsAtEnd())
			{
				//Determine points of traversing in the B direction
				currentIndex=itA.GetIndex();		

				for (int i=0;i<3;i++){
					pointB1[i]=min(max(floor(currentIndex[i]-hSize*B[i]+0.5),0),rSize[i]-1); 
					pointB2[i]=min(max(floor(currentIndex[i]+hSize*B[i]+0.5),0),rSize[i]-1);
				}

				//std::cout << "cI " << currentIndex << std::endl;
				//if (t){
				//	t=false;
				//	resampledImage->TransformIndexToPhysicalPoint(pointB1, PpointB1);
				//	resampledImage->TransformIndexToPhysicalPoint(pointB2, PpointB2);
				//	std::cout << "B1 " << convert_point(PpointB1) << std::endl;
				//	std::cout << "B2 " << convert_point(PpointB2) << std::endl;
				//}

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
			final2T=clock();

			time1+=(double)(final1T-start1T)/CLOCKS_PER_SEC;	
			time2+=(double)(final2T-start2T)/CLOCKS_PER_SEC;

			resampledImage->TransformIndexToPhysicalPoint(AIndex, maxPhysicalIndex);
			//std::cout << "A ind  " << maxPhysicalIndex << std::endl;	

			resampledImage->TransformIndexToPhysicalPoint(maxIndex, maxPhysicalIndex);	

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
			}

			//std::cout << "dot product " << dot_product << std::endl;

			//average vectors
			for (i=0;i<3;i++){
				//calculate magnitude of two vectors
				magnitude= sqrt(pow(E_mtx[i][0]+E_mtx_prev[i][0],2)+pow(E_mtx[i][1]+E_mtx_prev[i][1],2)
					+pow(E_mtx[i][2]+E_mtx_prev[i][2],2));

				E_mtx[i][0]=(E_mtx[i][0]+E_mtx_prev[i][0])/magnitude; 
				E_mtx[i][1]=(E_mtx[i][1]+E_mtx_prev[i][1])/magnitude; 
				E_mtx[i][2]=(E_mtx[i][2]+E_mtx_prev[i][2])/magnitude;
			}

			dot_product=E_mtx_prev[2][0]*E_mtx[2][0]+E_mtx_prev[2][1]*E_mtx[2][1]+E_mtx_prev[2][2]*E_mtx[2][2];

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

	//seedDirection is both input and output
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::refine_seed(PhysicalPointType unrefinedSeed, PhysicalPointType* refinedSeed, VectorType* seedDirection, PhysicalPointType* testResult,MType E_mtx1){	
			int i=0;
			MType E_mtx,E_mtx_prev; 
			VectorType rDirection;
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
			*refinedSeed=TwoD_Maximum(*refinedSeed, E_mtx); 
			*/
			*refinedSeed=TwoD_Maximum(unrefinedSeed, E_mtx1); 
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
				pPoint=TwoD_Maximum(pEst, E_mtx);
				E_mtx_prev=E_mtx;
			}

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
			int i;
			int nPtsInLine;
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
				return -1;
			}

			vtkIdType pointId=gotLine->GetId(pInd);	
			//std::cout << "PID " << pointId << std::endl;
			checkPts->GetPoint(pointId, gotPoint);

			(*returnPoint)[0]=-1*gotPoint[0];
			(*returnPoint)[1]=-1*gotPoint[1];
			(*returnPoint)[2]=gotPoint[2];

			//return returnPoint;
			return 1;
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
				checkInd++;
				gotPoint=get_point(checkLin,checkInd,&checkPoint);
				//i++;
			}

			newStart=checkPoint;
			prevCheck=checkPoint;
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

			while(refDist<dist){		
				gotPoint=get_point(checkLin,checkInd,&checkPoint);
				if (gotPoint!=1)
					break;
				refDist=compute_distance(prevCheck,checkPoint);
				checkInd++;				
			}

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
				return 1;	
			}
			//1  -> everything ok or not enough information								
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
					/*
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
		::termination(PhysicalPointType physicalPoint, MType E_matrix){
			IndexType point;
			float backCount=0;
			float dotProduct=0;
			bool terminate_=false;
			int i;

			image->TransformPhysicalPointToIndex(physicalPoint, point);
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

			float value = image->GetPixel( point );

			//termination by going in circles
			for (i=1; i<20;i++){
				if (distances[i]>distances[i-1]) //these are distances from seed
					backCount++;
			}

			if (backCount>8){
				terminate_=true;
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
			
			double getP[3];	
			int pointsLengthA=(int)pointsA->GetNumberOfPoints();
			int newPointsLength= min(numPoints,pointsLengthA);
			
			write_points(pointsA, pointsAPermanent, newPointsLength);	
			
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



} // end namespace itk
#endif



