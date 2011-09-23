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


#ifndef __itkIterateNeighborhoodCharacteristicDirectionsToPathFilter_txx
#define __itkIterateNeighborhoodCharacteristicDirectionsToPathFilter_txx

#include "itkIterateNeighborhoodCharacteristicDirectionsToPathFilter.h"

namespace itk
{
	
	template <class TInputImage, class TOutputPath>
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::IterateNeighborhoodCharacteristicDirectionsToPathFilter()
	{
		m_NbMaxIter = 5000;
	}
	
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::~IterateNeighborhoodCharacteristicDirectionsToPathFilter()
	{
	}
	
	
	template <class TInputImage, class TOutputPath>
	void
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
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
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GetNumberOfPathsToExtract() const
	{
		return m_EndPointList.size();
	}
	
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	const typename IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>::PointType &
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GetNextEndPoint()
	{
		return m_EndPointList[m_CurrentOutput];
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
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
			IndexType newPosition;
			IndexType currentIndex;
			input->TransformPhysicalPointToIndex( currentPoint, currentIndex);
			double bestValue = input->GetPixel( currentIndex );
			
			double dist2source = 0;
			for(unsigned int d = 0; d < SetDimension; d++)
			{
				dist2source += (currentPoint[d] - m_StartPoint[d])*(currentPoint[d] - m_StartPoint[d]);
			}
			dist2source = sqrt(dist2source);
			
			unsigned int count = 0;
			
			//Get the neighboring point with the smallest value
			typedef NeighborhoodIterator< InputImageType > NeighborhoodIterator;
			typename NeighborhoodIterator::RadiusType radius;
			radius.Fill(1);
				
			NeighborhoodIterator it(radius, input, input->GetRequestedRegion());
			
			while (dist2source > 0 && count < m_NbMaxIter) {
				
				// Add point as vertex in path
				output->AddVertex( currentIndex );

				it.SetLocation( currentIndex );
				for(unsigned int i = 0; i < it.Size(); ++i)
				{

					IndexType neighborPosition = it.GetIndex(i);
					
					bool isNeighborInsideRegion = true;
					
					for (unsigned int d = 0; d < SetDimension; d++) {
						if (neighborPosition[d] < m_StartIndex[d]) 
						{
							isNeighborInsideRegion = false;
							continue;
						}
						
						if (neighborPosition[d] > m_LastIndex[d]) 
						{
							isNeighborInsideRegion = false;
							continue;
						}
					}
					
					if ( isNeighborInsideRegion ) 
					{
						double neighborValue = input->GetPixel( neighborPosition );
						if( neighborValue < bestValue )
						{
							bestValue = neighborValue;
							newPosition = neighborPosition;
						}
					}
				}
				
				currentIndex = newPosition;
				
				input->TransformIndexToPhysicalPoint( currentIndex, currentPoint );
				
				// Update the fucking distance to the source
				dist2source = 0;
				for(unsigned int d = 0; d < SetDimension; d++)
				{
					dist2source += (currentPoint[d] - m_StartPoint[d])*(currentPoint[d] - m_StartPoint[d]);
				}
				dist2source = sqrt(dist2source);
				count++;
			}

			// Add point as vertex in path
			output->AddVertex( newPosition );
			
    }
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	void 
	IterateNeighborhoodCharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "NumberOfEndPoints: "		<< m_EndPointList.size() << std::endl;
		os << indent << "CurrentOutput"					<< m_CurrentOutput << std::endl;
		os << indent << "NbMaxIter"							<< m_NbMaxIter << std::endl;
		
	}
	
}

#endif