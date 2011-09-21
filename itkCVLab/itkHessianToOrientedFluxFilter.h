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

#ifndef __itkHessianToOrientedFluxFilter_h
#define __itkHessianToOrientedFluxFilter_h

#include "itkImageToImageFilter.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImage.h"
#include "itkProgressAccumulator.h"
#include "itkNeighborhoodAlgorithm.h"

namespace itk
{
	
	/** \class HessianToOrientedFluxFilter
	 * \brief This filter takes as input a hessian of an image
	 * and convolve it with the oriented flux kernel.
	 * PixelType of the input image is supposed to be SymmetricSecondRankTensor
	 * the output image type is the same.
	 *
	 * \author : Fethallah Benmansour
	 */
	template <typename TInputImage, 
	typename TOutputImage= TInputImage >
	class ITK_EXPORT HessianToOrientedFluxFilter:
	public ImageToImageFilter<TInputImage,TOutputImage>
	{
	public:
		/** Standard class typedefs. */
		typedef HessianToOrientedFluxFilter												Self;
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
		
		/**Declare types for neighborhoods */
		typedef ConstShapedNeighborhoodIterator
		< InputImageType >																				ShapedNeighborhoodIteratorType;
		typedef typename ShapedNeighborhoodIteratorType
		::ConstIterator																						ConstIterator;
		typedef ConstNeighborhoodIterator
		< InputImageType >																				NeighborhoodIteratorType;
		typedef typename 
		ShapedNeighborhoodIteratorType::RadiusType								RadiusType;
		typedef typename 
		ShapedNeighborhoodIteratorType::OffsetType								OffsetType;
		
		/**declare types for regions */
		typedef typename InputImageType::RegionType               InputImageRegionType;
    typedef typename OutputImageType::RegionType              OutputImageRegionType;
		
		/** Run-time type information (and related methods).   */
		itkTypeMacro( HessianToOrientedFluxFilter, ImageToImageFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Set/Get min and max radius */
		void SetMinRadius( RealType );
		itkGetConstMacro(MinRadius, RealType);
		void SetMaxRadius( RealType );
		itkGetConstMacro(MaxRadius, RealType);
		
		/** Get effective maximal radius */
		itkGetConstMacro(EffectiveRadius, RealType);
		
		/** HessianToOrientedFluxFilter needs all of the input to produce an
		 * output. Therefore, HessianToOrientedFluxFilter needs to provide
		 * an implementation for GenerateInputRequestedRegion in order to inform
		 * the pipeline execution model.
		 * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
		virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

	protected:
		
		HessianToOrientedFluxFilter();
		virtual ~HessianToOrientedFluxFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;
		
		/** Generate the shaped neighborhood iterator*/
		//ShapedNeighborhoodIteratorType GenerateShapedNeighborhood();
		
		/** Generate Data */
		void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
															int threadId);
		
		// Override since the filter produces the entire dataset
		void EnlargeOutputRequestedRegion(DataObject *output);
		
		
	private:
		
		HessianToOrientedFluxFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		RealType																	m_MinRadius;
		RealType																	m_MaxRadius;
		RealType																	m_EffectiveRadius;
		
	};
	
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHessianToOrientedFluxFilter.txx"
#endif

#endif
