/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkMultiScaleOrientedFluxBasedMeasureFFTImageFilter.h,v $
 Language:  C++
 Date:      $Date: 2010-02-04 12:16:43 $
 Version:   $Revision: 1.10 $
 
 This class has been obtained by modifiying the 
 itkMultiScaleHessianBasedMeasureImageFilter.h file of the ITK library 
 distributed and coyrighted by Insight Software Consortium.  
 The reason for modifications is to enable calling the SetSigma method of the 
 OrientedFluxToMeasureFilter (if it exists) for each scale level. Another reason is 
 to output an (N+1)-D output hessian-based measure image instead of an N-D one.
 Modified By:	Engin Turetken
 Date:			17.10.2010
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __itkMultiScaleOrientedFluxBasedMeasureFFTImageFilter_h
#define __itkMultiScaleOrientedFluxBasedMeasureFFTImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkDivideByConstantImageFilter.h>
#include <itkFFTOrientedFluxMatrixImageFilter.h>
#include <itkTimeProbe.h>

namespace itk
{
	/**\class MultiScaleOrientedFluxBasedMeasureFFTImageFilter
	 * \brief A filter to enhance structures using Hessian eigensystem-based 
	 * measures in a multiscale framework
	 * 
	 * The filter evaluates a Hessian-based enhancement measure, such as vesselness 
	 * or objectness, at different scale levels. The Hessian-based measure is computed 
	 * from the Hessian image at each scale level and the best response is selected. 
	 *
	 * Minimum and maximum sigma value can be set using SetMinSigma and SetMaxSigma
	 * methods respectively. The number of scale levels is set using 
	 * SetNumberOfSigmaSteps method. Exponentially distributed scale levels are 
	 * computed within the bound set by the minimum and maximum sigma values 
	 * 
	 * The filter computes a second output image (accessed by the GetScaleOutput method)
	 * containing the scales at which each pixel gave the best reponse. 
	 *
	 * \author Luca Antiga Ph.D.  Medical Imaging Unit,
	 *                            Bioengineering Deparment, Mario Negri Institute, Italy.
	 *
	 * \sa HessianToObjectnessMeasureImageFilter 
	 * \sa Hessian3DToVesselnessMeasureImageFilter 
	 * \sa HessianSmoothed3DToVesselnessMeasureImageFilter 
	 * \sa HessianRecursiveGaussianImageFilter 
	 * \sa SymmetricEigenAnalysisImageFilter
	 * \sa SymmetricSecondRankTensor
	 * 
	 * \ingroup IntensityImageFilters TensorObjects
	 *
	 */
	template <typename TInputImage,
	typename THessianImage, 
	typename TScaleImage,
	typename TOrientedFluxToMeasureFilter,
	typename TOutputNDImage = Image<typename NumericTraits<typename TInputImage::PixelType>::ScalarRealType, 
	::itk::GetImageDimension<TInputImage>::ImageDimension> >
	class ITK_EXPORT MultiScaleOrientedFluxBasedMeasureFFTImageFilter 
	: public ImageToImageFilter< TInputImage, TOutputNDImage > 
	{
	public:
		/** Standard class typedefs. */
		typedef MultiScaleOrientedFluxBasedMeasureFFTImageFilter									Self;
		typedef ImageToImageFilter<TInputImage, TOutputNDImage>										Superclass;
		typedef SmartPointer<Self>																								Pointer;
		typedef SmartPointer<const Self>																					ConstPointer;
		
		typedef typename NumericTraits
		<typename TInputImage::PixelType>::ScalarRealType													RealType;
		typedef TInputImage																												InputImageType;
		typedef TOutputNDImage																										OutputNDImageType;
		typedef THessianImage																											HessianImageType;
		typedef FFTOrientedFluxMatrixImageFilter<InputImageType>									OrientedFluxFilterType;
		typedef TScaleImage																												ScaleImageType;
		typedef TOrientedFluxToMeasureFilter																			OrientedFluxToMeasureFilterType;
		
		typedef typename OrientedFluxFilterType::OutputImageType									OrientedFluxImageType;
		typedef typename OrientedFluxImageType::Pointer														OrientedFluxImagePointer;
		
		/** Declare Addition filter for OrientedFluxImageType */
		typedef AddImageFilter<OrientedFluxImageType, OrientedFluxImageType>			AddImageFilterType;
		typedef typename AddImageFilterType::Pointer															AddImageFilterPointer;
		
		/** Declare DivideByConstant filter for OrientedFluxImageType */
		typedef DivideByConstantImageFilter
		<OrientedFluxImageType, double, OrientedFluxImageType>										DivideByCstImageFilterType;
		typedef typename DivideByCstImageFilterType::Pointer											DivideByCstImageFilterPointer;
		
		/** Declare outputs types */
		typedef	Image<typename OutputNDImageType::PixelType, 
		::itk::GetImageDimension<OutputNDImageType>::ImageDimension + 1>					OutputNPlus1DImageType;
		typedef Image<typename HessianImageType::PixelType, 
		::itk::GetImageDimension<HessianImageType>::ImageDimension + 1>						NPlus1DHessianImageType;
		
		typedef typename TInputImage::PixelType																		InputPixelType;
		typedef typename TInputImage::RegionType																	InputRegionType;
		typedef typename TOutputNDImage::PixelType																OutputNDPixelType;
		typedef typename TOutputNDImage::RegionType																OutputNDRegionType;
		typedef typename OutputNPlus1DImageType::RegionType												OutputNPlus1DRegionType;
		
		typedef ImageToImageFilterDetail::ImageRegionCopier<OutputNPlus1DImageType::ImageDimension,
		InputImageType::ImageDimension>																						InputToOutputRegionCopierType;
		
		/** Image dimension. */
		itkStaticConstMacro(ImageDimension, unsigned int, 
												::itk::GetImageDimension<InputImageType>::ImageDimension);
		
		/** Types for Scale image */
		typedef typename ScaleImageType::PixelType																ScalePixelType;
		
		/** Hessian computation filter. */
		typedef HessianRecursiveGaussianImageFilter
		< InputImageType, HessianImageType>																				HessianFilterType;
		
		/** Update image buffer that holds the best objectness response. This is not redundant from
		 the output image because the latter may not be of float type, which is required for the comparisons 
		 between responses at different scales. */ 
		typedef Image< double, itkGetStaticConstMacro(ImageDimension) >						UpdateBufferType;
		typedef typename UpdateBufferType::ValueType															BufferValueType;
    
		typedef typename Superclass::DataObjectPointer														DataObjectPointer;
		
		typedef FFTOrientedFluxMatrixImageFilter< InputImageType, HessianImageType > FFTOrientedFluxType;
		typedef typename OutputNDImageType::Pointer															 OutputNDImagePointer;
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Runtime information support. */
		itkTypeMacro(MultiScaleOrientedFluxBasedMeasureFFTImageFilter, 
                 ImageToImageFilter);
		
		/** Set/Get macros for SigmaMin */
		itkSetMacro(SigmaMinimum, double);
		itkGetConstMacro(SigmaMinimum, double);
		
		/** Set/Get macros for SigmaMax */
		itkSetMacro(SigmaMaximum, double);
		itkGetConstMacro(SigmaMaximum, double);
		
		/** Set/Get macros for Number of Scale */
		void SetNumberOfSigmaSteps(unsigned int );
		//itkSetMacro(NumberOfSigmaSteps, unsigned int);
		itkGetConstMacro(NumberOfSigmaSteps, unsigned int);
		
		/** 
		 * Methods to turn on/off flag to use a fixed sigma value for computing 
		 * the Hessian images. By default, it is false, which means we use the 
		 * given set of sigma values.
		 */
		/* Not needed anymore since it would be always the case */
		/*
		 itkSetMacro(UseAFixedSigmaForComputingHessianImage,bool);
		 itkGetConstMacro(UseAFixedSigmaForComputingHessianImage,bool);
		 itkBooleanMacro(UseAFixedSigmaForComputingHessianImage);
		 */
		
		/** 
		 * Set/Get macros for the fixed sigma value of the Hessian image.
		 * This parameter is used only if the 
		 * UseAFixedSigmaForComputingHessianImage flag is true.
		 */
		itkSetMacro(FixedSigmaForHessianImage, double);
		itkGetConstMacro(FixedSigmaForHessianImage, double);
		
		/** Get the image containing the Hessian computed at the best
		 * response scale */
		HessianImageType* GetHessianOutput();
		
		/** Get the (N+1)-D image containing the Hessian computed at all the
		 * scales */
		NPlus1DHessianImageType* GetNPlus1DHessianOutput();
		
		/** Get the image containing the scales at which each pixel gave the
		 * best response */
		ScaleImageType* GetScaleOutput();
		
		/** Get the (N+1)-D image containing the hessian based measure
		 * responses at all scales. */
		OutputNPlus1DImageType* GetNPlus1DImageOutput();
		
		void EnlargeOutputRequestedRegion (DataObject *);
		
		/** Methods to turn on/off flag to generate an image with scale values at
		 *  each pixel for the best vesselness response */
		itkSetMacro(GenerateScaleOutput,bool);
		itkGetConstMacro(GenerateScaleOutput,bool);
		itkBooleanMacro(GenerateScaleOutput);
		
		/** 
		 * Methods to turn on/off flag to generate an image with hessian 
		 * matrices at each pixel for the best vesselness response 
		 */
		itkSetMacro(GenerateHessianOutput,bool);
		itkGetConstMacro(GenerateHessianOutput,bool);
		itkBooleanMacro(GenerateHessianOutput);
		
		/** 
		 * Methods to turn on/off flag to generate the (N+1)-D image with 
		 * hessian matrices at each pixel for all possible scales. 
		 */
		itkSetMacro(GenerateNPlus1DHessianOutput,bool);
		itkGetConstMacro(GenerateNPlus1DHessianOutput,bool);
		itkBooleanMacro(GenerateNPlus1DHessianOutput);
		
		/** 
		 * Methods to turn on/off flag to treat the structures as bright or dark. 
		 * Its value is true by default.
		 */
		itkSetMacro(BrightObject,bool);
		itkGetConstMacro(BrightObject,bool);
		itkBooleanMacro(BrightObject);
		
		/** Methods to turn on/off flag to generate an image with hessian-based objectness 
		 * measure values at each pixel. */
		itkSetMacro(GenerateNPlus1DHessianMeasureOutput,bool);
		itkGetConstMacro(GenerateNPlus1DHessianMeasureOutput,bool);
		itkBooleanMacro(GenerateNPlus1DHessianMeasureOutput);
		
		/** This is overloaded to create the Scale and Hessian output images */
		virtual DataObjectPointer MakeOutput(unsigned int idx);
		
	protected:
		MultiScaleOrientedFluxBasedMeasureFFTImageFilter();
		~MultiScaleOrientedFluxBasedMeasureFFTImageFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;
		
		
		virtual void CallCopyInputRegionToOutputRegion(OutputNPlus1DRegionType &destRegion,
																									 const InputRegionType &srcRegion);
		virtual void GenerateOutputInformation();
		
		/** Generate Data */
		void GenerateData( void );
		
	private:
		void UpdateMaximumResponse(double sigma, unsigned int scaleLevel);
		double ComputeSigmaValue(int scaleLevel);
		
		void AllocateUpdateBuffer();
		
		//purposely not implemented
		MultiScaleOrientedFluxBasedMeasureFFTImageFilter(const Self&); 
		void operator=(const Self&); //purposely not implemented
		
		double																						m_SigmaMinimum;
		double																						m_SigmaMaximum;
		unsigned int																			m_NumberOfSigmaSteps;
		std::vector< RealType >														m_Sigmas;
		
		double																						m_FixedSigmaForHessianImage;
		//typename OrientedFluxToMeasureFilterType::Pointer	m_OrientedFluxToMeasureFilter;
		std::vector<typename OrientedFluxToMeasureFilterType::Pointer>		m_OrientedFluxToMeasureFilterList;
		typename UpdateBufferType::Pointer								m_UpdateBuffer;
		
		bool																							m_GenerateScaleOutput;
		bool																							m_GenerateHessianOutput;
		bool																							m_GenerateNPlus1DHessianOutput;	
		bool																							m_GenerateNPlus1DHessianMeasureOutput;
		
		bool																							m_BrightObject;
		
	};
	
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleOrientedFluxBasedMeasureFFTImageFilter.txx"
#endif

#endif
