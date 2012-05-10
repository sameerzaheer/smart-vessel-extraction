/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkSmartVesselExtraction.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSmartVesselExtraction_hxx
#define __itkSmartVesselExtraction_hxx

#include "itkSmartVesselExtraction.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

/** Constructor */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap> 
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
::SmartVesselExtraction()
{
  m_DistanceMap = 0;
  
  // when set to true it will be a bit faster, but it will result in minimizing
  // the sum of distances^4 instead of the sum of distances^2
  m_ComputeSquaredDistance = false; 

	//std::cout << "SmartVesselExtraction() ..." << std::endl;
}

/** Return the number of values, i.e the number of points in the moving set */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
unsigned int
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>  
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
typename SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>::MeasureType
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
::GetValue( const TransformParametersType & parameters ) const
{
	MeasureType measure;
	measure.Fill(0.0);

  return measure;

}

/** Get the Derivative Measure */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
void
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
::GetDerivative( const TransformParametersType & itkNotUsed(parameters),
                 DerivativeType & itkNotUsed(derivative) ) const
{

}

/** Get both the match Measure and theDerivative Measure  */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
void
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
::GetValueAndDerivative(const TransformParametersType & parameters, 
                        MeasureType & value, DerivativeType  & derivative) const
{
  value = this->GetValue(parameters);
  this->GetDerivative(parameters,derivative);
}



/** PrintSelf method */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
void
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "SmartVesselExtraction:PrintSelf()"<< std::endl;

  ////os << indent << "DistanceMap: " << m_DistanceMap << std::endl;
  //if(m_ComputeSquaredDistance)
  //  {
  //  os << indent << "m_ComputeSquaredDistance: True"<< std::endl;
  //  }
  //else
  //  {
  //  os << indent << "m_ComputeSquaredDistance: False"<< std::endl;
  //  }
}


//
////******************************************************************************
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
int
SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
::Test( int argc, char *argv[] )
{
	std::cout << "SmartVesselExtraction::Test()" << std::endl;

	//TestSmartVesselExtraction( );
	//return 0;

	return 0 ;
}
//******************************************************************************



////
//////******************************************************************************
//template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
//int
//SmartVesselExtraction<TFixedPointSet,TMovingPointSet,TDistanceMap>
//::MainLoop( int argc, char *argv[] )
//{
//	std::cout << "SmartVesselExtraction::MainLoop()" << std::endl;
//
//	//TestSmartVesselExtraction( );
//	//return 0;
//
//	return 0 ;
//}
////******************************************************************************
//
//





} // end namespace itk
#endif



