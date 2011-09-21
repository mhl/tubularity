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


#ifndef __itkHessianToOrientedFluxFilter_txx
#define __itkHessianToOrientedFluxFilter_txx

#include "itkHessianToOrientedFluxFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
	
	
	/**
	 * Constructor
	 */
	template <typename TInputImage, typename TOutputImage >
	HessianToOrientedFluxFilter<TInputImage,TOutputImage>
	::HessianToOrientedFluxFilter()
	{
		m_MinRadius = 0.0;
		m_MaxRadius = 1.0;
		m_EffectiveRadius = 0.0;
	}
	
	/**
	 * set min Radius
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxFilter<TInputImage,TOutputImage>
	::SetMinRadius( RealType minRadius )
	{
		m_MinRadius = minRadius;
		
		this->Modified();
	}
	
	/**
	 * set max Radius
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxFilter<TInputImage,TOutputImage>
	::SetMaxRadius( RealType maxRadius )
	{
		m_MaxRadius = maxRadius;
		
		this->Modified();
	}
	
	//
	//
	//
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxFilter<TInputImage,TOutputImage>
	::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
	{
		// call the superclass' implementation of this method. this should
		// copy the output requested region to the input requested region
		Superclass::GenerateInputRequestedRegion();
		
		// This filter needs all of the input
		typename HessianToOrientedFluxFilter<TInputImage,TOutputImage>::InputImagePointer image = const_cast<InputImageType *>( this->GetInput() );
		if (image)
    {
			image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
    }
	}
	
	
	//
	//
	//
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxFilter<TInputImage,TOutputImage>
	::EnlargeOutputRequestedRegion(DataObject *output)
	{
		TOutputImage *out = dynamic_cast<TOutputImage*>(output);
		
		if (out)
    {
			out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
	}
	
	/**
	 * Compute filter for Gaussian kernel
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	HessianToOrientedFluxFilter<TInputImage,TOutputImage >
	::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
												 int threadId)
	{
		//Get Input and Output
		InputImageConstPointer input  = this->GetInput();
		OutputImagePointer     output = this->GetOutput();
		
		if (m_MinRadius >= m_MaxRadius) 
		{
			itkExceptionMacro("min radius larger than max radius !!?");
		}
		
		if ( input.IsNull() )
		{
			itkExceptionMacro( "Input image must be provided" );
		}
		
		//Define Boundary conditions
		ZeroFluxNeumannBoundaryCondition<InputImageType>      nbc;
		
		// Radius and spacing
		RadiusType radius;
		SpacingType spacing = input->GetSpacing();

		/** Define the radius of the shaped iterator for the current scale */
		for (unsigned int d = 0; d < ImageDimension; d++)
		{
			radius[d] = floor(m_MaxRadius / spacing[d]) + 1;
		}
		
		m_EffectiveRadius = 0.0;
		// Find the data-set boundary "faces"
		typename NeighborhoodAlgorithm::
		ImageBoundaryFacesCalculator<InputImageType>::
		FaceListType																					faceList;
		NeighborhoodAlgorithm::
		ImageBoundaryFacesCalculator<InputImageType>					bC;		
		faceList = bC(input, outputRegionForThread, radius);
		
		typename NeighborhoodAlgorithm::
		ImageBoundaryFacesCalculator<InputImageType>::
		FaceListType::iterator																fit;
		
		// support progress methods/callbacks
		ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
		
		// Define output iterator
		ImageRegionIterator<OutputImageType>                  outputIt;
		
		// define a neighborhood iterator useful for the shaped one 
		NeighborhoodIteratorType	 it( radius, input, input->GetRequestedRegion() );
		it.NeedToUseBoundaryConditionOff();
		
		// Process each of the boundary faces.  These are N-d regions which border
		// the edge of the buffer.
		for (fit=faceList.begin(); fit != faceList.end(); ++fit)
    {
			outputIt = ImageRegionIterator<OutputImageType>(output, *fit);
			ShapedNeighborhoodIteratorType bit( radius, input, *fit );
			bit.OverrideBoundaryCondition(&nbc);
			bit.ClearActiveList();
			
			for (unsigned int j = 0; j < it.Size(); j++) 
			{
				OffsetType offset = it.GetOffset( j );
				RealType distance = 0.0;
				//compute the distance to the center
				for (unsigned int d = 0; d < ImageDimension; d++) 
				{
					distance += vnl_math_sqr( static_cast<RealType>(offset[d]) * spacing[d] );
				}
				distance = sqrt(distance);
				if (distance <= m_MaxRadius && distance > m_MinRadius)
				{
					m_EffectiveRadius = vnl_math_max(m_EffectiveRadius, distance);
					bit.ActivateOffset(offset);
				}	
			}
			
			bit.GoToBegin();
			outputIt.GoToBegin();
			while ( ! bit.IsAtEnd() )
      {
				ConstIterator ci;
				PixelType pixel;
				pixel.SetIdentity();
				pixel *= 0.0;
				for (ci = bit.Begin(); ci != bit.End(); ci++)
        {
					pixel += ci.Get();
				}
				outputIt.Set(pixel);
				
				++bit;
				++outputIt;
				progress.CompletedPixel();
			}
		}
	}
	
	
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxFilter<TInputImage,TOutputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);
		os << indent << "Minimal Radius: " << std::endl
		<< this->m_MinRadius << std::endl;
		os << indent << "Maximal Radius: " << std::endl
		<< this->m_MaxRadius << std::endl;
	}
	
	
} // end namespace itk

#endif
