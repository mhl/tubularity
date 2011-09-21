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

#ifndef __itkHessianToOrientedFluxMainCurvatureMeasureFilter_h
#define __itkHessianToOrientedFluxMainCurvatureMeasureFilter_h

#include "itkImageToImageFilter.h"
#include "itkHessianToOrientedFluxFilter.h"
#include "itkSymmetricEigenAnalysisImageFilter2.h"
#include "itkImage.h"

namespace itk
{
	
	/** \class HessianToOrientedFluxMainCurvatureMeasureFilter
	 * \brief This filter takes as input a hessian of an image
	 * and convolve it with the oriented flux kernel and outputs the most inportant eigen value.
	 * PixelType of the input image is supposed to be SymmetricSecondRankTensor
	 *
	 * \author : Fethallah Benmansour
	 */
	template <typename TInputImage, 
	typename TOutputImage = itk::Image<typename NumericTraits<
	typename TInputImage::PixelType::ValueType>::RealType, 
	TInputImage::ImageDimension > >
	class ITK_EXPORT HessianToOrientedFluxMainCurvatureMeasureFilter:
	public ImageToImageFilter<TInputImage,TOutputImage>
	{
	public:
		/** Standard class typedefs. */
		typedef HessianToOrientedFluxMainCurvatureMeasureFilter												Self;
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
		
		typedef HessianToOrientedFluxFilter<InputImageType>				HessianToOFluxFilterType;
		typedef typename HessianToOFluxFilterType::Pointer				HessianToOFluxFilterPointer;
		
		typedef typename HessianToOFluxFilterType::OutputImageType OrientedFluxImageType;
		
		typedef Image< FixedArray<RealType, ImageDimension>, ImageDimension > EigenValueImageType;
		
		typedef SymmetricEigenAnalysisImageFilter2<OrientedFluxImageType, EigenValueImageType> EigenAnalysisImageFilterType;
		typedef typename EigenAnalysisImageFilterType::Pointer						EigenAnalysisImageFilterPointer;
		
		/** Run-time type information (and related methods).   */
		itkTypeMacro( HessianToOrientedFluxMainCurvatureMeasureFilter, ImageToImageFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		void SetSigma( RealType );
		RealType GetSigma( ) const;
		
		void SetBrightObject( bool bIsBrightObject );
		bool GetBrightObject() const;
		
		/** HessianToOrientedFluxMainCurvatureMeasureFilter needs all of the input to produce an
		 * output. Therefore, HessianToOrientedFluxMainCurvatureMeasureFilter needs to provide
		 * an implementation for GenerateInputRequestedRegion in order to inform
		 * the pipeline execution model.
		 * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
		virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);
		
	protected:
		
		HessianToOrientedFluxMainCurvatureMeasureFilter();
		virtual ~HessianToOrientedFluxMainCurvatureMeasureFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;
		
		/** Generate the shaped neighborhood iterator*/
		//ShapedNeighborhoodIteratorType GenerateShapedNeighborhood();
		
		/** Before Threaded Generate Data */
		void BeforeThreadedGenerateData( );
		
		/** Threaded Generate Data */
		void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
															int threadId );
		
		// Override since the filter produces the entire dataset
		void EnlargeOutputRequestedRegion(DataObject *output);
		
	private:
		
		HessianToOrientedFluxMainCurvatureMeasureFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		bool																			m_IsBright;
		RealType																	m_Sigma;
		HessianToOFluxFilterPointer								m_OFFilter;
		EigenAnalysisImageFilterPointer						m_eigenAnalysisFilter;
		
	};
	
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHessianToOrientedFluxMainCurvatureMeasureImageFilter.txx"
#endif

#endif
