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


#ifndef __itkFFTOrientedFluxMatrixImageFilter_txx
#define __itkFFTOrientedFluxMatrixImageFilter_txx

#include "itkFFTOrientedFluxMatrixImageFilter.h"

#include <itkImageFileWriter.h>

// todo get rid of this header for the insight submission
#include <omp.h>

namespace itk
{
	template <typename TInputImage, typename TOutputImage > SimpleFastMutexLock FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>::m_Mutex;
	
	/**
	 * Constructor
	 */
	template <typename TInputImage, typename TOutputImage >
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>
	::FFTOrientedFluxMatrixImageFilter()
	{
		m_Sigma0 = 1.0;
		m_Radius = 1.0;
		
		m_GreatestPrimeFactor = 13;
		
		m_ImageAdaptor = OutputImageAdaptorType::New();
	}
	
	/**
	 * Set Sigma0
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>
	::SetSigma0( RealType sigma0 )
	{
		if(m_Sigma0 != sigma0)
		{
			m_Sigma0 = sigma0;
			this->Modified();
		}
	}
	
	/**
	 * Get Sigma0
	 */
	template <typename TInputImage, typename TOutputImage >
	typename FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>::RealType
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>
	::GetSigma0( )
	{
		return	m_Sigma0;
	}
	
	/**
	 * Set Radius
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>
	::SetRadius( RealType radius )
	{
		if(m_Radius != radius)
		{
			m_Radius = radius;
			this->Modified();
		}
	}

	/**
	 * Get Radius
	 */
	template <typename TInputImage, typename TOutputImage >
	typename FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>::RealType
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>
	::GetRadius( )
	{
		return m_Radius;
	}
	
	/***************************************************************************************
	 *  For 2 given directions, Generates the oriented flux matix kernel in the fourier
	 *  domain as Given by Eq.8 in:
	 *  Max W. K. Law and Albert C. S. Chung, 
	 *	“Three Dimensional Curvilinear Structure Detection using Optimally Oriented Flux”
	 *  The Tenth European Conference on Computer Vision, (ECCV’ 2008)
	 *
	 * \author Fethallah Benmansour
	 ***************************************************************************************/
	template <typename TInputImage, typename TOutputImage >
	void
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage >
	::GenerateOrientedFluxMatrixElementKernel(ComplexImagePointerType &kernel,
																					const ComplexImageType* input, 
																					unsigned int derivA, unsigned int derivB, 
																					float radius, float sigma0)
	{
		// i and j should both be smaller than the dimension
		if (derivA >= ImageDimension || derivB >= ImageDimension)
		{
			itkGenericExceptionMacro("Derivatives along the dimensions, these indices should be less than the dimension");
		}
		
		RegionType region = input->GetLargestPossibleRegion();
		
		typedef  typename InputImageType::SpacingType		SpacingType;
		
		SpacingType imageSpacing = input->GetSpacing();
		SizeType    size    = region.GetSize();
		IndexType  start    = region.GetIndex();
		
		SpacingType freqSpacing;
		PointType  origin;
		for(unsigned int i = 0; i < ImageDimension; i++)
		{
			origin[i] = -0.5/ imageSpacing[i];
			
			if(i != 0)
			{
				freqSpacing[i] = (1.0 / static_cast<double>(size[i]-1)) / imageSpacing[i];
			}
			else
			/**
			 * we assume that fftw flags are activated.
			 */
			{
				freqSpacing[i] = (0.5 / static_cast<double>(size[i]-1)) / imageSpacing[i];
			}
		}
		
		kernel = ComplexImageType::New();
		kernel->CopyInformation( input );
		kernel->SetSpacing( freqSpacing );
		kernel->SetOrigin(origin);
		kernel->SetBufferedRegion( input->GetBufferedRegion() );
		kernel->Allocate();
		/**
		 * 2 corners per dimension except for the first one.
		 */
		unsigned int nbCorners = pow(float(2),float( ImageDimension-1));
		
		
		typedef std::vector<IndexType> listOfIndexCorners;
		typedef std::vector<PointType> listOfPointCorners;
		listOfIndexCorners listOfIdxC;
		listOfPointCorners listOfPtsC;
		listOfIdxC.clear();
		listOfPtsC.clear();
		// Get the corners (2D) *
		//*--------------/
		//---------------/
		//---------------/
		//---------------/
		//---------------/
		//---------------/
		//*--------------/
		
		for(unsigned int i = 0; i < nbCorners; i++)
		{
			IndexType corner;
			unsigned int idx = 0;
			corner[idx] = start[idx]; 
			idx++;
			unsigned int k = i;
			for(unsigned int j = 0; j < ImageDimension-1; j++)
			{
				if( k%2 )
				{
					corner[idx] = start[idx];
				}
				else
				{
					corner[idx] = start[idx]+size[idx]-1;
				}
				idx++;
				k /= 2;
			}
			listOfIdxC.push_back(corner);
		}
			
			

		
		// Get the physical location of these corners
		PointType pt;
		for(unsigned int i = 0; i < nbCorners; i++)
		{
			kernel->TransformIndexToPhysicalPoint(listOfIdxC[i], pt);
			listOfPtsC.push_back(pt);
		}
		
		// Start computation of the filter after some initializations 
		double eps = itk::NumericTraits<float>::epsilon();	
		double Ui = 0.0, Uj = 0.0;
		double normU;
		
		typedef itk::ImageRegionIterator< ComplexImageType > ImageIterator;
		ImageIterator  it( kernel, kernel->GetLargestPossibleRegion() );
		it.GoToBegin();
		
		while(!it.IsAtEnd())
		{
			//DO the compoutation of the Kernel 
			IndexType index  = it.GetIndex();
			PointType point;
			kernel->TransformIndexToPhysicalPoint(index, point);
			//compute norm of U
			
			normU = 1e9;
			for(unsigned int i = 0; i < nbCorners; i++)
			{
				double distToCorner = point.EuclideanDistanceTo(listOfPtsC[i]);
				if (normU > distToCorner)
				{
					normU = distToCorner;
					Ui = point[derivA] - (listOfPtsC[i])[derivA];
					Uj = point[derivB] - (listOfPtsC[i])[derivB];
				}
			}
			
			ComplexPixelType value(0.0, 0.0);
			if(normU >= eps)
			{
				double phase = 2.0 * vnl_math::pi * radius * normU;
				
				if(ImageDimension == 2)
				{// scale normalized response
					value  += 2.0 * vnl_math::pi *exp( -2.0* (vnl_math::pi * normU * sigma0)*(vnl_math::pi * normU * sigma0) )*
					Ui * Uj * j1(phase) / normU;
				}
				else if(ImageDimension == 3)
				{
					value  += 4.0 *  vnl_math::pi * exp( -2.0* (vnl_math::pi * normU * sigma0)*(vnl_math::pi * normU * sigma0) ) *
					Ui * Uj * (cos(phase) - sin(phase)/phase) / (radius* normU*normU);
				}
				else
				{
					itkGenericExceptionMacro("Oriented Flux filter in the Fourier domain not imlemented for other dimensions than 2 and 3");
				}
			}
			
			it.Set(value);
			++it;
		}
		//
		typedef ImageFileWriter<ComplexImageType> WriterType;
		typename WriterType::Pointer writer  = WriterType::New();
		if(derivA == 2 && derivB == 2)
		{
			writer->SetInput( kernel );
			writer->SetFileName("kernel.nrrd");
			writer->Update();
		}
	}
	
	/**
	 * Generate data:
	 * Conpute the convolution componnent per componnent in the spectral domain 
	 *
	 * \author: F. Benmansour
	 */
	template <typename TInputImage, typename TOutputImage >
	void
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage >
	::GenerateData( )
	{
		InputImageConstPointer input  = this->GetInput();
		OutputImagePointer		 output = this->GetOutput();
		
		if( input.IsNull() )
		{
			itkExceptionMacro("Input image must be provided");
		}

		// Prepare Image adaptor
		m_ImageAdaptor->SetImage( this->GetOutput() );
		m_ImageAdaptor->SetLargestPossibleRegion( this->GetInput()->GetLargestPossibleRegion() );
		m_ImageAdaptor->SetBufferedRegion( this->GetInput()->GetBufferedRegion() );
		m_ImageAdaptor->SetRequestedRegion( this->GetInput()->GetRequestedRegion() );
		m_ImageAdaptor->Allocate();
		
		// Fake a kernel(wrt to the input scale) and pad the input image accordingly
		// TODO: this one might be replaced by an other padding function that does not require a kernel but just the input image and the pdding size
		float minSpacing   = input->GetSpacing().GetVnlVector().min_value();
		int kernelPixelSize = Math::Round<int>( this->GetRadius() / minSpacing ) + 1;
		IndexType index; 
		index.Fill(-kernelPixelSize);
		SizeType size;
		size.Fill(2*kernelPixelSize+1);
		RegionType region;
		region.SetIndex(index); 
		region.SetSize(size);
		
		InputImagePointer fakeKernel = InputImageType::New();
		fakeKernel->SetBufferedRegion(region);
		fakeKernel->Allocate();
		fakeKernel->FillBuffer(0.0);
		typedef FFTPadImageFilter2< InputImageType, InputImageType, InternalImageType, InternalImageType > PadType;
		typename PadType::Pointer pad = PadType::New();
		pad->SetInput( input );
		pad->SetInputKernel( fakeKernel );
		pad->SetReleaseDataFlag( true );
		pad->SetGreatestPrimeFactor( this->GetGreatestPrimeFactor() );
		
		typename FFTFilterType::Pointer fft = FFTFilterType::New();
		fft->SetInput( pad->GetOutput() );
		
		// todo uncomment for the insight submission
//		fft->SetNumberOfThreads( this->GetNumberOfThreads() );
		
		// todo comment for the insight submission
		fft->SetNumberOfThreads( omp_get_max_threads() );
		
		// This call to the FFT library has to be locked for 
		// threaded (private) access to different instances of this object.
		m_Mutex.Lock();
		fft->Update();
		m_Mutex.Unlock();
		
		unsigned int element = 0;
		for(unsigned int i = 0; i < ImageDimension; i++)
		{
			for(unsigned int j = i; j < ImageDimension; j++)
			{
				
				ComplexImagePointerType kernel = NULL;
				
				GenerateOrientedFluxMatrixElementKernel( kernel, fft->GetOutput(), i, j, this->GetRadius(), this->GetSigma0() );

				typedef itk::MultiplyImageFilter< ComplexImageType,
				ComplexImageType,
				ComplexImageType > MultType;
				typename MultType::Pointer mult = MultType::New();
				mult->SetInput( 0, fft->GetOutput() );
				mult->SetInput( 1, kernel );
				mult->SetNumberOfThreads( this->GetNumberOfThreads() );
				mult->SetReleaseDataFlag( true );
				mult->SetInPlace( false );
				mult->Update();
				
				typename IFFTFilterType::Pointer ifft = IFFTFilterType::New();
				ifft->SetInput( mult->GetOutput() );
				ifft->SetActualXDimensionIsOdd( false );

				// todo uncomment for the insight submission
//				ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
				
				// todo comment for the insight submission
				ifft->SetNumberOfThreads( omp_get_max_threads() );
				
				ifft->SetReleaseDataFlag( true );
				
				// This call to the FFT library has to be locked for 
				// threaded (private) access to different instances of this object.
				m_Mutex.Lock();
				ifft->Update();
				m_Mutex.Unlock();
								
				typedef ClipImageFilter< InternalImageType, InternalImageType > ClipType;
				typename ClipType::Pointer clip = ClipType::New();
				clip->SetInput( ifft->GetOutput() );
				clip->SetNumberOfThreads( this->GetNumberOfThreads() );
				clip->SetReleaseDataFlag( true );
				clip->SetInPlace( true );
				
				typedef itk::RegionFromReferenceImageFilter< InternalImageType, InternalImageType > CropType;
				typename CropType::Pointer crop = CropType::New();
				crop->SetInput( clip->GetOutput() );
				crop->SetReferenceImage( input );
				crop->SetNumberOfThreads( this->GetNumberOfThreads() );
				crop->Update();
				
				ImageRegionIteratorWithIndex< InternalImageType > it(crop->GetOutput(), crop->GetOutput()->GetRequestedRegion());
				m_ImageAdaptor->SelectNthElement( element++ );
				ImageRegionIteratorWithIndex< OutputImageAdaptorType > ot( m_ImageAdaptor, m_ImageAdaptor->GetRequestedRegion());
				
				it.GoToBegin();
				ot.GoToBegin();
				while( !it.IsAtEnd() )
				{
					ot.Set( it.Get() );
					++it;
					++ot;
				}
			}
		}
	}
	
	template <typename TInputImage, typename TOutputImage>
	void
	FFTOrientedFluxMatrixImageFilter<TInputImage,TOutputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);
		os << indent << "Sigma0 (for smoothing): " << std::endl
		<< this->m_Sigma0 << std::endl;
		os << indent << "Radius: " << std::endl
		<< this->m_Radius << std::endl;
		os << indent << "ImageAdaptor: " << std::endl
		<< this->m_ImageAdaptor << std::endl;
	}
} // end namespace itk

#endif
