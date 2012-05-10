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
		beta= 0.4;
		sigma =2.5;
		sFactor=0.2;

		//other parameters
		max_iterations=140;
		epsilon =0.8;
		gradTolerance=0.3;
		windowSize =5;

		//seed model parameters
		maxSeeds=100;
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

		//Interpolator
		interpolator1 = InterpolatorType1::New();

		//To be Objectified
		angleThresh=0.58;
		distanceThresh=1.0;

		termCircles=false;

		//Gradient filter (for radius estimation)
		typename gradientType::Pointer gradientFilter = gradientType::New();

		//interpolator for radius can be linear
		typename InterpolatorType2::Pointer interpolator2 = InterpolatorType2::New();


		//std::cout << "VesselExtractionTools() ..." << std::endl;
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

			//********Resample image********
			IndexType pointIndex;
			PhysicalPointType max_physicalIndex;		

			float re_scale = scale*spacing[0]; //physical radius of vessel
			origin[0] = physicalIndex[0]-re_scale-re_spacing[0]; //re_spacing[0]->extra padding
			origin[1] = physicalIndex[1]-re_scale-re_spacing[0];
			origin[2] = physicalIndex[2]-re_scale-re_spacing[0];
			filter->SetOutputOrigin( origin );	 

			ImageType::SizeType   r_size;
			r_size[0]=re_scale*2/re_spacing[0]+2; //re_scale/re_spacing[0]=scale in voxels
			r_size[1]=re_scale*2/re_spacing[0]+2; //+2 to deal with rounding errors
			r_size[2]=re_scale*2/re_spacing[0]+2;
			filter->SetSize( r_size );  
			filter->Update();  
			resampled_image=filter->GetOutput();
			resampled_image->TransformPhysicalPointToIndex(physicalIndex, pointIndex);	

			//std::cout << "r_size   : "  << r_size<< std::endl;
			//std::cout << "re size   : "  << resampled_image->GetLargestPossibleRegion().GetSize()<< std::endl;

			//********Extract eigen vectors********
			IndexType pointA1, pointB1, pointA2, pointB2, currentIndex, maxIndex, AIndex, BInd1, BInd2;
			float A[3], B[3]; 
			float mag;
			float h_size=re_scale/(re_spacing[0]*sqrt(float(2)));
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
				pointA1[i]=min(max(floor(pointIndex[i]-h_size*A[i]+0.5),0),r_size[i]-1); 
				pointA2[i]=min(max(floor(pointIndex[i]+h_size*A[i]+0.5),0),r_size[i]-1);
			}

			//std::cout << "hasz " << h_size << std::endl;
			//std::cout << "A1 " << pointA1 << std::endl;
			//std::cout << "A2 " << pointA2 << std::endl;

			ImageLineItType itA(resampled_image, pointA1, pointA2); //iterator in directionA
			itA.GoToBegin();
			PixelType currentValue, maxValue=-1;
			currentValue=(PixelType)itA.Get();	
			maxIndex=pointIndex;

			while (!itA.IsAtEnd())
			{
				//Determine points of traversing in the B direction
				currentIndex=itA.GetIndex();		

				for (int i=0;i<3;i++){
					pointB1[i]=min(max(floor(currentIndex[i]-h_size*B[i]+0.5),0),r_size[i]-1); 
					pointB2[i]=min(max(floor(currentIndex[i]+h_size*B[i]+0.5),0),r_size[i]-1);
				}

				//std::cout << "cI " << currentIndex << std::endl;
				//std::cout << "B1 " << pointB1 << std::endl;
				//std::cout << "B2 " << pointB2 << std::endl;

				ImageLineItType itB(resampled_image, pointB1, pointB2);

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
			resampled_image->TransformIndexToPhysicalPoint(AIndex, max_physicalIndex);
			//std::cout << "A ind  " << max_physicalIndex << std::endl;	

			resampled_image->TransformIndexToPhysicalPoint(maxIndex, max_physicalIndex);	

			return max_physicalIndex;
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
		::flipSign(MType E_mtx){

			//Extract eigen vectors
			E_mtx[2][0]=E_mtx[2][0]*(-1); 
			E_mtx[2][1]=E_mtx[2][1]*(-1); 
			E_mtx[2][2]=E_mtx[2][2]*(-1); 	

			return E_mtx;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::changeSign(MType E_mtx_prev, MType E_mtx){

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
				E_mtx=flipSign(E_mtx);
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
		::correctOrientation(MType E_mtx, VectorType direction_){ 	

			float dot_product;
			dot_product=direction_[0]*E_mtx[2][0]+direction_[1]*E_mtx[2][1]+direction_[2]*E_mtx[2][2];

			if (dot_product<0){
				E_mtx=flipSign(E_mtx);	
				//std::cout << "orientation corrected ******" << std::endl;	
			}

			return E_mtx;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	typename VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>::MType 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::correctOrientation(PhysicalPointType curPoint, PhysicalPointType refPoint, MType E_mtx){ 	

			int i;
			VectorType refTOcurrent;

			for (i=0;i<3;i++)
				refTOcurrent[i]=curPoint[i]-refPoint[i]; 

			return correctOrientation(E_mtx, refTOcurrent);
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
			E_mtx=correctOrientation(E_mtx, *seedDirection); //with respect to direction
			//Get 2-D Max  
			*refinedSeed=TwoD_Maximum(*refinedSeed, E_mtx); 
			*/
			*refinedSeed=TwoD_Maximum(unrefinedSeed, E_mtx1); 
			//Set-up
			E_mtx_prev=compute_EigenMatrix(HessianFilter->GetOutput(), (*refinedSeed));
			E_mtx_prev=correctOrientation(E_mtx_prev, *seedDirection); //with respect to direction

			pPoint=*refinedSeed;	

			//**Estimate next point**
			for (i=0;i<11;i++){
				E_mtx=compute_EigenMatrix(HessianFilter->GetOutput(), pPoint);  
				E_mtx=changeSign(E_mtx_prev,E_mtx);
				E_mtx=correctOrientation(E_mtx, *seedDirection); //with respect to direction
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
			/*
			std::cout << "sc: " << maxScore << " dst: " << maxScoreDist << std::endl;	
			std::cout << " BSd: " << bestSeed << " BDir: " << maxDirections << std::endl;	

			std::cout << "sc: " << maxScore << " BTest: " << bestTest << std::endl;
			std::cout << " SeedP: " << seedPoint  << "   Score: " << maxScore << std::endl;
			std::cout << " ExpSeed: " << expectedSeed << "  VD: " << v_dir <<std::endl;
			std::cout << " BSeed: " << bestSeed << " BDir: " << maxDirections << std::endl;
			*/
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
		::initialize_distances(){
			for (int i=0; i<20;i++)
				distances[i]=0;	
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
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::start_checks(int n){	
			PhysicalPointType pointA=infoVector[n].seed;
			PhysicalPointType pointB=infoVector[n].newSeed;
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
				//std::cout << "P " << checkPoint << " Ang " << cosAng << std::endl;
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
			prevTrackedI=2;
			/*
			std::cout << "S " <<pointA << std::endl;
			std::cout << "St " <<checkPoint << std::endl;	
			std::cout << "NS " <<pointB << std::endl;
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
			PhysicalPointType checkPoint;
			//checkInd++; //raise by at least one
			//	std::cout << " D " << dist << " rD1 " << refDist << " CI1 " << checkInd ;

			while(get_point(checkLin,checkInd,&checkPoint)==1 && refDist<dist){		
				refDist=compute_distance(prevCheck,checkPoint);
				checkInd++;
			}
			distRemainder=refDist-dist;
			//	std::cout << " rD2 " <<  refDist << " CI2 " <<checkInd << std::endl;

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

			if (cosAngle>-2){
				if (cosAngle<angleThresh){
					nextSeed=prevCheck;
					nextDirection=checkPoint;
					/*std::cout << "cA        " << cosAngle << std::endl;

					std::cout << "C1 " <<prevCheck<< std::endl;
					std::cout << "C2 " <<checkPoint<< std::endl;
					std::cout << "T1 " <<prevTracked << std::endl;	
					std::cout << "T2 " <<currTracked << std::endl;*/
					prevCheck=checkPoint;	
					return -1;
				}else{
					/*	std::cout << "cA        " << cosAngle << std::endl;
					std::cout << "C1 " <<prevCheck<< std::endl;
					std::cout << "C2 " <<checkPoint<< std::endl;
					std::cout << "T1 " <<prevTracked << std::endl;	
					std::cout << "T2 " <<currTracked << std::endl;	
					*/
					prevCheck=checkPoint;	
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
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	float  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::compute_distances(PhysicalPointType curPoint, PhysicalPointType refPoint){
			float distance_;
			distance_=compute_distance(curPoint,refPoint);
			//distance_=sqrt(pow(refPoint[0]-curPoint[0],2)+pow(refPoint[1]-curPoint[1],2)+pow(refPoint[2]-curPoint[2],2));
			distances[0]=distance_;

			//store last 10 distances
			for (int i=19; i>0;i--){
				distances[i]=distances[i-1];
			}
			return distance_;
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	bool  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::termination(PhysicalPointType physicalPoint, MType E_matrix){
			IndexType point;
			float back_count=0;
			bool terminate_=false;

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
			for (int i=1; i<20;i++){
				if (distances[i]>distances[i-1]) //these are distances from seed
					back_count++;
			}
			if (back_count>8){
				terminate_=true;
				//	std::cout << "Term via Circles" << std::endl;
			}
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
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::MAIN_LOOP(int n, vtkPoints *centerPoints, vtkFloatArray* radiusArray){
			int i,startI=1,j, goodPoints;
			int returnInt=0;
			MType E_mtx_prev, E_mtx;		
			PhysicalPointType physicalEstimation, physicalPoint=infoVector[n].seedU;	
			float radius_, radius_prev=infoVector[n].radiusS;  
			float distance_, distance__;
			float mag=0, cosAngle=0;
			initialize_distances();
			bool endLoop=false;
			bool toContinue=false;
			termCircles=false;

			//***INITIALIZE***
			//start_checks(n);
			if (firstTime){
				centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);
				radiusArray->InsertNextValue(radius_prev);
				startI=2;
				firstTime=false;
			}

			///Find a new seed
			//and check if its good
			if (optimize_seed(infoVector[n].seedU, infoVector[n].v_direction, &(infoVector[n].newSeed)) !=1){
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

			E_mtx_prev=compute_EigenMatrix(HessianFilter->GetOutput(), physicalPoint);
			E_mtx_prev=correctOrientation(E_mtx_prev,infoVector[n].v_direction);

			//***THE LOOP***
			for (i=startI; i<max_iterations && !endLoop; i++){ //2
				//**Find next cross section**
				E_mtx=compute_EigenMatrix(HessianFilter->GetOutput(), physicalPoint);

				//Ensure the next point has correct orientation
				E_mtx=changeSign(E_mtx_prev,E_mtx); //with respect to previous point
				//with respect to overall direction
				E_mtx=correctOrientation(E_mtx, infoVector[n].v_direction); 

				//**Estimate next point**
				physicalEstimation=move_Point(physicalPoint, E_mtx);

				//**Find the next point**
				physicalPoint=TwoD_Maximum(physicalEstimation, E_mtx);	
				//Record the point
				centerPoints->InsertNextPoint (physicalPoint[0]*(-1), physicalPoint[1]*(-1), physicalPoint[2]);
				//std::cout << "PP (" << i << ") : " << physicalPoint << std::endl;	
				//**Find the radius**
				radius_=compute_radius(physicalPoint, radius_prev, spacing[0]*2, E_mtx);
				radiusArray->InsertNextValue(radius_);
				radius_prev=radius_;
				//scale=max(1,min(2,radius_/(2*spacing[0])));
				//std::cout << "Radius " << radius_ << std::endl;					
				//Save the cross-section and direction
				E_mtx_prev=E_mtx; 		

				distance__=compute_distance(physicalPoint,infoVector[n].seed);
				distance_=compute_distances(physicalPoint,infoVector[n].direction); //or new_seed

				//Determine to end
				if(termination(physicalPoint, E_mtx)){
					endLoop=true;
					termCircles=true;
					i=max(startI+1,i-10);
				}		
				if(distance__>infoVector[n].length){
					endLoop=true;
					//std::cout << "L " << distance__ << "  " << infoVector[n].length << std::endl;	
				}	
				/*
				if (check_check(physicalPoint,i)){
				endLoop=true;
				toContinue=true;
				std::cout << "O " << std::endl;	
				}*/

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
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::check_line_orientiation(int n, vtkPoints *inPoints, vtkPoints *outPoints, int termCircle){
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

			start_checks(n);

			//reorder the first points
			for (i=2; i<total_points;i++){			
				inPoints->GetPoint(i,getP);	
				thePoint[0]=getP[0]*(-1);
				thePoint[1]=getP[1]*(-1);
				thePoint[2]=getP[2];
				//std::cout<< "thePoint " << thePoint << std::endl;

				decide=check_check(thePoint,i);

				if (decide==0 && termCircle==1){
					decide=-1;
				}
				if (decide==-1){
					total_points=prevTrackedI;
					//std::cout << " Ptracked " << total_points << std::endl;
					infoVector[n].seedU=nextSeed;
					infoVector[n].direction=nextDirection;
					compute_v_dir(n);
					returnValue=1;
					break;
				} else if (decide==0){
					break;
				}
			}

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

			return check_line_orientiation(n,inPoints,outPoints,0);

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

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::write_points(vtkPoints *inPoints, vtkPoints *outPoints){
			double getP[3];	
			int total_points=(int)inPoints->GetNumberOfPoints();

			//reorder the first points
			for (int i=0; i<total_points;i++){			
				inPoints->GetPoint(i,getP);	
				outPoints->InsertNextPoint (getP[0], getP[1],getP[2]);
			}

	}
	//Should make smoothing function more accurate by taking into account the distance between points
	//Currently, the computation assumes a distance =1
	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::smooth_points(vtkPoints *centerPoints, vtkPoints *fixedPoints, int total_points){
			double a[3], getP[3], av[3];
			double weights[11]={1,0.956,0.835,0.667,0.487,0.325,0.198,0.110,0.0563,0.0262,0.0112}; //Gaussian weighted
			//c=3.335
			double weight_sum;	

			for (int i=0; i<total_points;i++){	
				weight_sum=0;
				a[0]=a[1]=a[2]=0;
				//Look at the point itself
				centerPoints->GetPoint(i,getP);
				for (int k=0;k<3;k++)
					a[k]+=weights[0]*getP[k];
				weight_sum+=weights[0];
				//Look at its neighbours
				for (int j=-1*windowSize; j<0;j++){
					if((i+j)>=0 && (i-j)<total_points){
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

			int total_points=(int)centerPoints->GetNumberOfPoints();
			smooth_points(centerPoints,fixedPoints,total_points);
	}

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::writePointsVTK(const char *filename, vtkPoints *points, vtkCellArray *centerLine, vtkFloatArray* radiusArray){

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
	int
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::read_file(char * argv[]){

			const char * filename = argv[2];
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
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::read_data(char * argv[]){

			vtkPolyDataReader *reader = vtkPolyDataReader::New();
			reader->SetFileName( argv[3] ); //.vtk
			reader->Update();
			//Get the necessary data
			vtkPolyData *polyData =  reader->GetOutput(); 
			checkPts = polyData->GetPoints();
			checkLines = polyData->GetLines();
			//vtkIdType nLines = lines->GetNumberOfCells();

			return 0;

	}


	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	void 
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
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

			//Distances
			distance_tolerance=5*spacing[0];
			for (int i=0; i<20;i++)
				distances[i]=0;	

			//filter	
			filter->SetInput( image );
			filter->SetInterpolator( interpolator1 );
			direction.SetIdentity();  
			filter->SetDefaultPixelValue( 1.0 );  
			re_spacing=spacing*sFactor;

			const VectorType re_spacingConst(re_spacing);
			filter->SetOutputSpacing( re_spacingConst );  

			//filter->SetOutputSpacing( re_spacing );  
			filter->SetOutputDirection( direction );

			//Random generator
			randomGenerator = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
			randomGenerator->Initialize(2861);
			//2861


	}

	//int VesselExtraction(int argc, char * argv[]);

	template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
	int  
		VesselExtractionTools<TFixedPointSet,TMovingPointSet,TDistanceMap>
		::VesselExtraction(int argc, char * argv[])
	{
		std::cout << " VesselExtractionTools::VesselExtraction()... "<< std::endl;
		if(argc != 4)
		{
			std:: cout << "Usage: SmartVesselExtraction <image filename> <seed filename> <lines filename>" << std::endl;
			std:: cout << "Example: SmartVesselExtraction EHPS3R3Crop.mhd EHPS2103BP806Tfm3_2in.txt EHPS2103BP806Tfm3.vtk" << std::endl;
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
		INITIALIZE();
		radiusArray->SetNumberOfComponents(1);
		radiusArray->SetName("RadiusArray");  

		//Get seed and target points
		if (read_file(argv)==-1){
			exit(-1);
		} 

		if (read_data(argv)==-1){

			exit(-1);
		} 
		std::cout << " read input2" << std::endl;
		
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



} // end namespace itk
#endif



