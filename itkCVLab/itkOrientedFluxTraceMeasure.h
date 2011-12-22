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

#ifndef __itkOrientedFluxTraceMeasureFilter_h
#define __itkOrientedFluxTraceMeasureFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
	
	/** \class OrientedFluxTraceMeasureFilter
	 * \brief This filter takes as input the oriented flux response of an image
	 * and computes the trace.
	 * PixelType of the input image is supposed to be SymmetricSecondRankTensor
	 * TODO: give this filter a more generic name, 
	 * since this filter is doing nothing but computing traces
	 *
	 * \author : Fethallah Benmansour
	 */
	template <typename TInputImage, 
	typename TOutputImage = itk::Image<typename NumericTraits<
	typename TInputImage::PixelType::ValueType>::RealValueType, 
	TInputImage::ImageDimension > >
	class ITK_EXPORT OrientedFluxTraceMeasureFilter:
	public ImageToImageFilter<TInputImage,TOutputImage>
	{
	public:
		/** Standard class typedefs. */
		typedef OrientedFluxTraceMeasureFilter										Self;
		typedef ImageToImageFilter<TInputImage,TOutputImage>			Superclass;
		typedef SmartPointer<Self>																Pointer;
		typedef SmartPointer<const Self>													ConstPointer;
		
		/** Pixel Type of the input image */
		typedef TInputImage																				InputImageType;
		typedef typename InputImageType::ConstPointer							InputImageConstPointer;
		typedef typename InputImageType::Pointer									InputImagePointer;
		typedef typename InputImageType::PixelType								PixelType;
		typedef typename NumericTraits<PixelType>::ValueType			RealType;
		typedef typename InputImageType::SpacingType							SpacingType;
		
		/** Image dimension. */
		itkStaticConstMacro(ImageDimension, unsigned int,
												::itk::GetImageDimension<TInputImage>::ImageDimension);
		
		/** Type of the output Image */
		typedef TOutputImage                                      OutputImageType;
		typedef typename OutputImageType::Pointer									OutputImagePointer;
		typedef typename OutputImageType::PixelType								OutputPixelType;
		typedef typename PixelTraits<OutputPixelType>::ValueType  OutputComponentType;
		
		/**declare types for regions */
		typedef typename InputImageType::RegionType               InputImageRegionType;
    typedef typename OutputImageType::RegionType              OutputImageRegionType;
		
		/** Run-time type information (and related methods).   */
		itkTypeMacro( OrientedFluxTraceMeasureFilter, ImageToImageFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		void SetBrightObject( bool bIsBrightObject );
		bool GetBrightObject() const;
		
		/** OrientedFluxTraceMeasureFilter needs all of the input to produce an
		 * output. Therefore, OrientedFluxTraceMeasureFilter needs to provide
		 * an implementation for GenerateInputRequestedRegion in order to inform
		 * the pipeline execution model.
		 * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
		virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);
		
	protected:
		
		OrientedFluxTraceMeasureFilter();
		virtual ~OrientedFluxTraceMeasureFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;
		
		/** Before Threaded Generate Data */
		//void BeforeThreadedGenerateData( );
		
		/** Threaded Generate Data */
		void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
															int threadId );
		
		// Override since the filter produces the entire dataset
		void EnlargeOutputRequestedRegion(DataObject *output);
		
	private:
		
		OrientedFluxTraceMeasureFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		bool																			m_IsBright;
		
	};
	
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOrientedFluxTraceMeasure.txx"
#endif

#endif
