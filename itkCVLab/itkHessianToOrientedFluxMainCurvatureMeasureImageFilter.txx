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


#ifndef __itkHessianToOrientedFluxMainCurvatureMeasureFilter_txx
#define __itkHessianToOrientedFluxMainCurvatureMeasureFilter_txx

#include "itkHessianToOrientedFluxMainCurvatureMeasureImageFilter.h"

namespace itk
{
	
	
	/**
	 * Constructor
	 */
	template <typename TInputImage, typename TOutputImage >
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
	::HessianToOrientedFluxMainCurvatureMeasureFilter()
	{
		m_Sigma = 1.0;
		m_IsBright = true;
	}
	
	/**
	 * Set Sigma
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
	::SetSigma( RealType sigma )
	{
		m_Sigma = sigma;
		
		this->Modified();
	}
	
	/**
	 * Get Sigma
	 */
	template <typename TInputImage, typename TOutputImage>
	typename HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>::RealType
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
	::GetSigma( ) const
	{
		return m_Sigma;
	}
	
	/**
	 * Set Bright Object
	 */
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
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
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
	::GetBrightObject( ) const
	{
		return m_IsBright;
	}
	
	//
	//
	//
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
	::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
	{
		// call the superclass' implementation of this method. this should
		// copy the output requested region to the input requested region
		Superclass::GenerateInputRequestedRegion();
		
		// This filter needs all of the input
		typename HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>::InputImagePointer image = const_cast<InputImageType *>( this->GetInput() );
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
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
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
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage >
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
		
		/**Eigen Analisys Filter */
		m_eigenAnalysisFilter = EigenAnalysisImageFilterType::New();
		m_eigenAnalysisFilter->SetInput( m_OFFilter->GetOutput() );
		m_eigenAnalysisFilter->SetGenerateEigenVectorImage(false);
		m_eigenAnalysisFilter->SetGenerateFirstEigenVectorOrientImage(false);
		m_eigenAnalysisFilter->Update();
		
		OutputImagePointer     output = this->GetOutput();
		output->CopyInformation(m_eigenAnalysisFilter->GetOutput());
		output->SetBufferedRegion(m_eigenAnalysisFilter->GetOutput()->GetBufferedRegion());
		output->Allocate();
		
	}
	
	/**
	 * Compute filter for Gaussian kernel
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage >
	::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
												 int threadId)
	{
		// support progress methods/callbacks
		ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
		
		ImageRegionIterator<OutputImageType>				 outputIt;
		typedef ImageRegionIterator< EigenValueImageType> EigenValueIteratorType;
		EigenValueIteratorType eigenIt;
		
		outputIt = ImageRegionIterator< OutputImageType >( this->GetOutput(), outputRegionForThread );
		eigenIt = EigenValueIteratorType( const_cast<EigenValueImageType*>(m_eigenAnalysisFilter->GetEigenValueImage()),
																		 outputRegionForThread );
		
		RealType normalizationFactor = pow(m_OFFilter->GetEffectiveRadius(), static_cast<RealType>(ImageDimension-1));
		
		outputIt.GoToBegin();
		eigenIt.GoToBegin();
		while ( !outputIt.IsAtEnd() ) 
		{
			RealType value = 0.0;
			if (m_IsBright) 
			{
				for (unsigned int i = 0; i < ImageDimension-1; i++) 
				{
					value -= eigenIt.Get()[i];
				}
			}
			else
			{
				for (unsigned int i = 1; i < ImageDimension; i++) 
				{
					value += eigenIt.Get()[i];
				}
			}
			
			value /= normalizationFactor;
	
// Allow negative responses and a higher dynamic range.	
// The following line is commented by eturetken on 27.05.2011.
			//value = vnl_math_max((double)value, (double)0.0);
			
			outputIt.Set( value );
			++outputIt;
			++eigenIt;
		}
	}
	
	
	template <typename TInputImage, typename TOutputImage>
	void
	HessianToOrientedFluxMainCurvatureMeasureFilter<TInputImage,TOutputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);
		os << indent << "Sigma: " << std::endl
		<< this->m_Sigma << std::endl;
	}
	
	
} // end namespace itk

#endif
