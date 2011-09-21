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


#ifndef __itkHessianToOrientedFluxTraceMeasureFilter_txx
#define __itkHessianToOrientedFluxTraceMeasureFilter_txx

#include "itkHessianToOrientedFluxTraceMeasureImageFilter.h"

namespace itk
{
	
	
	/**
	 * Constructor
	 */
	template <typename TInputImage, typename TOutputImage >
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::HessianToOrientedFluxTraceMeasureFilter()
	{
		m_Sigma = 1.0;
		m_IsBright = true;
	}
	
	/**
	 * Set Sigma
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::SetSigma( RealType sigma )
	{
		m_Sigma = sigma;
		
		this->Modified();
	}
	
	/**
	 * Get Sigma
	 */
	template <typename TInputImage, typename TOutputImage>
	typename HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>::RealType
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::GetSigma( ) const
	{
		return m_Sigma;
	}
	
	/**
	 * Set Bright Object
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
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
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::GetBrightObject( ) const
	{
		return m_IsBright;
	}
		
	//
	//
	//
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
	{
		// call the superclass' implementation of this method. this should
		// copy the output requested region to the input requested region
		Superclass::GenerateInputRequestedRegion();
		
		// This filter needs all of the input
		typename HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>::InputImagePointer image = const_cast<InputImageType *>( this->GetInput() );
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
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::EnlargeOutputRequestedRegion(DataObject *output)
	{
		TOutputImage *out = dynamic_cast<TOutputImage*>(output);
		
		if (out)
    {
			out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
	}
	
	/**
	 * Before Threaded Generate Data:
	 * Computes the Oriented Flux matrix and allocate the output 
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage >
	::BeforeThreadedGenerateData()
	{
		//Get Input and Output
		InputImageConstPointer input  = this->GetInput();
		
		if ( input.IsNull() )
		{
			itkExceptionMacro( "Input image must be provided" );
		}
		
		/** set the oriented flux filter */
		m_OFFilter = HessianToOFluxFilterType::New();
		m_OFFilter->SetInput( input );
		m_OFFilter->SetMinRadius(0.0);
		m_OFFilter->SetMaxRadius(m_Sigma);
		m_OFFilter->ReleaseDataFlagOn();
		
		m_OFFilter->Update();
		
		OutputImagePointer     output = this->GetOutput();
		output->CopyInformation(m_OFFilter->GetOutput());
		output->SetBufferedRegion(m_OFFilter->GetOutput()->GetBufferedRegion());
		output->Allocate();
		
	}
	
	/**
	 * Compute filter for Gaussian kernel
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage >
	::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
												 int threadId)
	{
		
		// support progress methods/callbacks
		ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
		
		
		ImageRegionIterator<OutputImageType>				 outputIt;
		typedef ImageRegionIterator< InputImageType> OrientedFluxIteratorType;
		OrientedFluxIteratorType orientedFluxIt;
		
		outputIt = ImageRegionIterator<OutputImageType>(this->GetOutput(), outputRegionForThread );
		orientedFluxIt = OrientedFluxIteratorType(m_OFFilter->GetOutput(), outputRegionForThread );
		
		RealType normalizationFactor = pow(m_OFFilter->GetEffectiveRadius(), static_cast<RealType>(ImageDimension-1));
		
		outputIt.GoToBegin();
		orientedFluxIt.GoToBegin();
		while ( !outputIt.IsAtEnd() ) 
		{
			RealType value = 0.0;
			for(unsigned int i = 0; i < ImageDimension; i++)
			{
				value += orientedFluxIt.Get()(i, i);
			}
			
			value /= normalizationFactor;
			
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
	HessianToOrientedFluxTraceMeasureFilter<TInputImage,TOutputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);
		os << indent << "Sigma: " << std::endl
		<< this->m_Sigma << std::endl;
	}
	
	
} // end namespace itk

#endif
