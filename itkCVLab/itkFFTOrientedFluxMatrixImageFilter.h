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

#ifndef __itkFFTOrientedFluxMatrixImageFilter_h
#define __itkFFTOrientedFluxMatrixImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkPixelTraits.h>
#include <itkNthElementImageAdaptor.h>
#include <itkFFTPadImageFilter.h>
#include <itkFFTPadImageFilter2.h>
#include <itkFFTRealToComplexConjugateImageFilter.h>
#include <itkFFTComplexConjugateToRealImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkClipImageFilter.h>
#include <itkRegionFromReferenceImageFilter.h>


namespace itk
{
	
	/** \class FFTOrientedFluxMatrixImageFilter
	 * \brief This filter takes as input an image
	 * and convolve it with the oriented flux kernel.
	 * The convolution being done in the spectral domain.
	 * PixelType of the input image is supposed to be scalar.
	 * 
	 * PixelType of the output image type is SymmetricSecondRankTensor.
	 *
	 * \author : Fethallah Benmansour
	 */
	template <typename TInputImage, 
	typename TOutputImage= Image< SymmetricSecondRankTensor< 
  ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType>::RealType,
  ::itk::GetImageDimension<TInputImage>::ImageDimension >,
	::itk::GetImageDimension<TInputImage>::ImageDimension >  >
	class ITK_EXPORT FFTOrientedFluxMatrixImageFilter:
	public ImageToImageFilter<TInputImage,TOutputImage>
	{
	public:
		/** Standard class typedefs. */
		typedef FFTOrientedFluxMatrixImageFilter									Self;
		typedef ImageToImageFilter<TInputImage,TOutputImage>			Superclass;
		typedef SmartPointer<Self>																Pointer;
		typedef SmartPointer<const Self>													ConstPointer;
		
		/** Pixel Type of the input image */
		typedef TInputImage																				InputImageType;
		typedef typename InputImageType::Pointer									InputImagePointer;
		typedef typename InputImageType::ConstPointer							InputImageConstPointer;
		typedef typename InputImageType::PixelType								PixelType;
		typedef typename NumericTraits<PixelType>::RealType				RealType;
		
		/** Define the image type for internal computations 
		 RealType is usually 'double' in NumericTraits. 
		 Here we prefer float in order to save memory.  */
		
		typedef float                                             InternalRealType;
		typedef Image<
    InternalRealType, 
    ::itk::GetImageDimension<TInputImage>::ImageDimension >		InternalImageType;
		
		/** Utilities for the input image */
		typedef typename InputImageType::IndexType								IndexType;
		typedef typename InputImageType::RegionType								RegionType;
		typedef typename InputImageType::SizeType									SizeType;
		typedef typename InputImageType::PointType								PointType;

		/** Image dimension. */
		itkStaticConstMacro(ImageDimension, unsigned int,
												::itk::GetImageDimension<TInputImage>::ImageDimension);
		
		/** FFT filters */
		typedef FFTRealToComplexConjugateImageFilter
		< InternalRealType, ImageDimension >											FFTFilterType;
		typedef FFTComplexConjugateToRealImageFilter
		< InternalRealType, ImageDimension >											IFFTFilterType;
		
		typedef typename FFTFilterType::TOutputImageType					ComplexImageType;
		typedef	typename ComplexImageType::PixelType							ComplexPixelType;
		typedef	typename ComplexImageType::Pointer								ComplexImagePointerType;
		
		
		/** Type of the output Image */
		typedef TOutputImage                                      OutputImageType;
		typedef typename OutputImageType::Pointer									OutputImagePointer;
		typedef typename OutputImageType::PixelType								OutputPixelType;
		typedef typename PixelTraits<OutputPixelType>::ValueType  OutputComponentType;
			
		/** Image adaptor */
		typedef NthElementImageAdaptor
		<OutputImageType, OutputComponentType>										OutputImageAdaptorType;
		typedef typename OutputImageAdaptorType::Pointer					OutputImageAdaptorPointer;
		
		
		/**declare types for regions */
		typedef typename InputImageType::RegionType               InputImageRegionType;
    typedef typename OutputImageType::RegionType              OutputImageRegionType;
		
		/** Run-time type information (and related methods).   */
		itkTypeMacro( FFTOrientedFluxMatrixImageFilter, ImageToImageFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Set/Get for the smoothing parameter \Sigma0 and for the scale Radius */
		void SetSigma0( RealType sigma0 );
		RealType GetSigma0( );
		void SetRadius( RealType radius);
		RealType GetRadius( );

		/**
		 * Set/Get the greatest prime factor allowed on the size of the padded image.
		 */
		itkGetConstMacro(GreatestPrimeFactor, int);
		itkSetMacro(GreatestPrimeFactor, int);
		
#ifdef ITK_USE_CONCEPT_CHECKING
		/** Begin concept checking */
		itkConceptMacro(InputHasNumericTraitsCheck,
										(Concept::HasNumericTraits<PixelType>));
		itkConceptMacro(OutputHasPixelTraitsCheck,
										(Concept::HasPixelTraits<OutputPixelType>));
		/** End concept checking */
#endif		
		
	protected:
		
		FFTOrientedFluxMatrixImageFilter();
		virtual ~FFTOrientedFluxMatrixImageFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;
		
		void GenerateOrientedFluxMatrixElementKernel(ComplexImagePointerType &kernel,
																								 const ComplexImageType* input, 
																								 unsigned int derivA, unsigned int derivB, 
																								 float radius, float sigma0);
		
		/** Generate Data */
		void GenerateData( );
		
	private:
		
		FFTOrientedFluxMatrixImageFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		RealType										m_Sigma0;
		RealType										m_Radius;
		
		int													m_GreatestPrimeFactor;
		
		OutputImageAdaptorPointer		m_ImageAdaptor;
		
	};
	
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTOrientedFluxMatrixImageFilter.txx"
#endif

#endif
