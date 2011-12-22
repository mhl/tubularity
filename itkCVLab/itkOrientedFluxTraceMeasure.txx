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


#ifndef __itkOrientedFluxTraceMeasureFilter_txx
#define __itkOrientedFluxTraceMeasureFilter_txx

#include "itkOrientedFluxTraceMeasure.h"

namespace itk
{
	
	
	/**
	 * Constructor
	 */
	template <typename TInputImage, typename TOutputImage >
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::OrientedFluxTraceMeasureFilter()
	{
		m_IsBright = true;
	}
	
	/**
	 * Set Bright Object
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::SetBrightObject( bool brightObj )
	{
		m_IsBright = brightObj;
		
		this->Modified();
	}
	
	/**
	 * Set Bright Object
	 */
	template <typename TInputImage, typename TOutputImage>
	bool
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::GetBrightObject( ) const
	{
		return m_IsBright;
	}
	
	//
	//
	//
	template <typename TInputImage, typename TOutputImage>
	void
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
	{
		// call the superclass' implementation of this method. this should
		// copy the output requested region to the input requested region
		Superclass::GenerateInputRequestedRegion();
		
		// This filter needs all of the input
		typename OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>::InputImagePointer image = const_cast<InputImageType *>( this->GetInput() );
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
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
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
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage >
	::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
												 int threadId)
	{
		
		// support progress methods/callbacks
		ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
		
		
		ImageRegionIterator<OutputImageType>				 outputIt;
		typedef ImageRegionConstIterator< InputImageType> OrientedFluxIteratorType;
		OrientedFluxIteratorType orientedFluxIt;
		
		outputIt = ImageRegionIterator<OutputImageType>(this->GetOutput(), outputRegionForThread );
		orientedFluxIt = OrientedFluxIteratorType(this->GetInput(), outputRegionForThread );
		
		outputIt.GoToBegin();
		orientedFluxIt.GoToBegin();
		while ( !outputIt.IsAtEnd() ) 
		{
			RealType value = 0.0;
			for(unsigned int i = 0; i < ImageDimension; i++)
			{
				value += orientedFluxIt.Get()(i, i);
			}
			
			if ( m_IsBright )
			{
				value = -value;
			}
			
			// Allow negative responses and a higher dynamic range.	
			// The following line is commented by eturetken on 27.05.2011.
			//			value = vnl_math_max((double)value, (double)0.0);
			outputIt.Set( value );
			++outputIt;
			++orientedFluxIt;
			progress.CompletedPixel();
		}
		
		
	}
	
	
	template <typename TInputImage, typename TOutputImage>
	void
	OrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);
		os << indent << "is Bright: " << std::endl
		<< this->m_IsBright << std::endl;
	}
	
	
} // end namespace itk

#endif
