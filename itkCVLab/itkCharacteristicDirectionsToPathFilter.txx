//**********************************************************
//Copyright 2011 Fethallah Benmansour
//
//Licensed under the Apache License, Version 2.0 (the "License"); 
//you may not use this file except in compliance with the License. 
//You may obtain a copy of the License at
//
//http://www.apache.org/licenses/LICENSE-2.0 
//
//Unless required by applicable law or agreed to in writing, software 
//distributed under the License is distributed on an "AS IS" BASIS, 
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
//See the License for the specific language governing permissions and 
//limitations under the License.
//**********************************************************


#ifndef __itkCharacteristicDirectionsToPathFilter_txx
#define __itkCharacteristicDirectionsToPathFilter_txx

#include "itkCharacteristicDirectionsToPathFilter.h"

namespace itk
{

template <class TInputImage, class TOutputPath>
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::CharacteristicDirectionsToPathFilter()
{
  m_TerminationDistance = 0.5;
  m_CurrentOutput = 0;
	m_Interpolator = InterpolatorType::New();
	m_Step = 1;
	m_NbMaxIter = 5000;
}


/**
 *
 */
template <class TInputImage, class TOutputPath>
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::~CharacteristicDirectionsToPathFilter()
{
}


template <class TInputImage, class TOutputPath>
void
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    InputImagePointer image =
      const_cast< InputImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

/**
 *
 */
template<class TInputImage, class TOutputPath>
unsigned int
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::GetNumberOfPathsToExtract() const
{
  return m_EndPointList.size();
}


/**
 *
 */
template<class TInputImage, class TOutputPath>
const typename CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>::PointType &
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::GetNextEndPoint()
{
  return m_EndPointList[m_CurrentOutput];
}


/**
 *
 */
template <class TInputImage, class TOutputPath>
void
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::SetStep( double step )
{
	
	SpacingType spacing;
	if ( this->GetInput() )
	{
		spacing = this->GetInput()->GetSpacing();
		double minSpacing  = spacing[0];
		for (unsigned int i = 1; i < SetDimension; i++)
		{
			if (minSpacing > spacing[i])
			{
				minSpacing = spacing[i];
			}
		}
		
		m_Step = step*minSpacing;
	}
	else {
		std::cerr << "set the input image first !!" << std::endl;
		exit(-2);
	}

}
	
/**
 *
 */
template <class TInputImage, class TOutputPath>
void
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::SetTerminationDistanceFactor( double factor )
{
	
	SpacingType spacing;
	if ( this->GetInput() )
	{
		spacing = this->GetInput()->GetSpacing();
		double maxSpacing  = spacing[0];
		for (unsigned int i = 1; i < SetDimension; i++)
		{
			maxSpacing = vnl_math_max(maxSpacing, spacing[i]);
		}
		
		m_TerminationDistance = factor*maxSpacing;
	}
	else {
		std::cerr << "set the input image first !!" << std::endl;
		exit(-2);
	}
	
}
	
	
/**
 *
 */
template <class TInputImage, class TOutputPath>
void
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::GenerateData( void )
{
  // Get the input
  InputImagePointer input = 
    const_cast< InputImageType * >( this->GetInput() );
	
	// cache some buffered region information
	m_StartIndex = input->GetRequestedRegion().GetIndex();
	m_LastIndex = m_StartIndex + input->GetRequestedRegion().GetSize();
	typename InputImageType::OffsetType offset;
	offset.Fill( 1 );
	m_LastIndex -= offset;
	
	
	
  if ( input.IsNull() )
	{
    itkExceptionMacro( "Input image must be provided" );
    return;
  }

  // Check the number of paths is not none
  unsigned int numberOfOutputs = GetNumberOfPathsToExtract();
  if ( numberOfOutputs == 0 )
    {
    itkExceptionMacro( "At least one path must be specified for extraction" );
    }
  this->ProcessObject::SetNumberOfRequiredOutputs( numberOfOutputs );
	
	// Set the interpolator input
	m_Interpolator->SetInputImage( this->GetInput() );
	
  // Do for each output
  for ( unsigned int n=0; n < numberOfOutputs; n++ )
    {
			// Set the output	index
			// NOTE: m_CurrentOutput is used in Execute() and GetNextEndPoint()
			m_CurrentOutput = n;
        
			// Make the output
			OutputPathPointer output
				= static_cast<TOutputPath*>( this->MakeOutput(n).GetPointer() ); 
			this->ProcessObject::SetNthOutput( n, output.GetPointer() );

			// Get the end point to back propagate from
			PointType currentPoint = this->GetNextEndPoint();
			double dist2source = 0;
			for(unsigned int d = 0; d < SetDimension; d++)
			{
				dist2source += (currentPoint[d] - m_StartPoint[d])*(currentPoint[d] - m_StartPoint[d]);
			}
			dist2source = sqrt(dist2source);
			

			unsigned int count = 0;
			
			//TODO : Add Oscialltion checker and getr rid of it
			
			while (dist2source > m_TerminationDistance && count < m_NbMaxIter) {
				
				// Convert currentPoint to continuous index
				ContinuousIndexType cindex;
				input->TransformPhysicalPointToContinuousIndex( currentPoint, cindex );
				
				for (unsigned int d = 0; d < SetDimension; d++) {
					if (cindex[d] < m_StartIndex[d]) 
					{
						cindex[d] = m_StartIndex[d];
					}
					
					if (cindex[d] > m_LastIndex[d]) 
					{
						cindex[d] = m_LastIndex[d];
					}
				}
				// Check that the new position is inside the image domain
				bool isCurrentPointInImageDomain = m_Interpolator->IsInsideBuffer( cindex );
				if ( !isCurrentPointInImageDomain ) 
				{
					itkExceptionMacro("Point out of image domain: should not happen");
				}
				
				// Add point as vertex in path
				output->AddVertex( cindex );
				
				//Get the steepest descent direction
				VectorType grad = m_Interpolator->EvaluateAtContinuousIndex(cindex);
				grad.Normalize();
				//Compute next point
				for (unsigned int d = 0; d < SetDimension; d++) {
					currentPoint[d] = currentPoint[d] - m_Step*grad[d];
					// chekc if current point is not out of requested region bounds
										
				}
				
				// Update the fucking distance to the source
				dist2source = 0;
				for(unsigned int d = 0; d < SetDimension; d++)
				{
					dist2source += (currentPoint[d] - m_StartPoint[d])*(currentPoint[d] - m_StartPoint[d]);
				}
				dist2source = sqrt(dist2source);
				count++;
			}
			
			if( count >=  m_NbMaxIter )
			{
				itkWarningMacro("Start point not reached, increase number of iterations");
			}
			
			// Convert currentPoint to continuous index
			ContinuousIndexType cindex;
			input->TransformPhysicalPointToContinuousIndex( m_StartPoint, cindex );
			
			
			for (unsigned int d = 0; d < SetDimension; d++) {
				if (cindex[d] < m_StartIndex[d]) 
				{
					cindex[d] = m_StartIndex[d];
				}
				
				if (cindex[d] > m_LastIndex[d]) 
				{
					cindex[d] = m_LastIndex[d];
				}
			}
			// Check that the new position is inside the image domain
			bool isCurrentPointInImageDomain = m_Interpolator->IsInsideBuffer( cindex );
			if ( !isCurrentPointInImageDomain ) 
			{
				itkExceptionMacro("Point out of image domain: should not happen");
			}
			// Add point as vertex in path
			output->AddVertex( cindex );
			
    }
}

/**
 *
 */
template<class TInputImage, class TOutputPath>
void 
CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "TerminationDistance: " << m_TerminationDistance << std::endl;
  os << indent << "NumberOfEndPoints: "		<< m_EndPointList.size() << std::endl;
	os << indent << "TerminationDistance"		<< m_TerminationDistance << std::endl;
	os << indent << "CurrentOutput"					<< m_CurrentOutput << std::endl;
	os << indent << "Interpolator"					<< m_Interpolator << std::endl;
	os << indent << "Step"									<< m_Step << std::endl;
	os << indent << "NbMaxIter"							<< m_NbMaxIter << std::endl;
	
}
	
}

#endif