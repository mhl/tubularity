/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianBasedMeasureImageFilter2.txx,v $
  Language:  C++
  Date:      $Date: 2009-08-26 19:09:35 $
  Version:   $Revision: 1.13 $
  
  This class has been obtained by modifiying the 
  itkMultiScaleHessianBasedMeasureImageFilter.txx file of the ITK library 
  distributed and coyrighted by Insight Software Consortium.  
  The reason for modifications is to enable calling the SetSigma method of the 
	HessianToMeasureFilter (if it exists) for each scale level.  Another reason is 
	to output an (N+1)-D output hessian-based measure image instead of an N-D one.
  Modified By:	Engin Turetken
  Date:			17.10.2010

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiScaleHessianBasedMeasureImageFilter2_txx
#define __itkMultiScaleHessianBasedMeasureImageFilter2_txx

#include "itkMultiScaleHessianBasedMeasureImageFilter2.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

template <int I>
struct Int2Type
{
  enum { value = I };
};

// Check if the class THessianToMeasureFilter has a member function with one of 
// the following  prototypes:
// void SetSigma(double);
// void SetSigma(const double&);
// void SetSigma(float);
// void SetSigma(const float&);
template<typename THessianToMeasureFilter>
struct DoesHessianToMeasureFilterHaveSetSigmaMethod
{
		template<typename U, void (U::*)(double)> struct SFINAE_double {};
		template<typename U, void (U::*)(const double&)> struct SFINAE_const_double {};
		template<typename U, void (U::*)(float)> struct SFINAE_float {};
		template<typename U, void (U::*)(const float&)> struct SFINAE_const_float {};
        template<typename U> static char Test(SFINAE_double<U, &U::SetSigma>*);
		template<typename U> static char Test(SFINAE_const_double<U, &U::SetSigma>*);
		template<typename U> static char Test(SFINAE_float<U, &U::SetSigma>*);
		template<typename U> static char Test(SFINAE_const_float<U, &U::SetSigma>*);		
        template<typename U> static int Test(...);
        static const bool Has = sizeof(Test<THessianToMeasureFilter>(0)) == sizeof(char);
};

template <typename THessianToMeasureFilterPtr>
void CallHessianToMeasureFilterSigmaMethod(THessianToMeasureFilterPtr filter, double sigma, Int2Type<true>)
{
	// Call Set sigma method of the HessianToMeasureFilter
	filter->SetSigma( sigma );
}

template <typename THessianToMeasureFilterPtr>
void CallHessianToMeasureFilterSigmaMethod(THessianToMeasureFilterPtr filter, double sigma, Int2Type<false>)
{
	// Do nothing.
}

namespace itk
{

/**
 * Constructor
 */
template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::MultiScaleHessianBasedMeasureImageFilter2()
{
  m_SigmaMinimum = 0.2;
  m_SigmaMaximum = 2.0;

  m_NumberOfSigmaSteps = 10;

  m_HessianFilter = HessianFilterType::New();
  m_HessianToMeasureFilter = NULL;

  //Instantiate Update buffer
  m_UpdateBuffer = UpdateBufferType::New();

	m_UseAFixedSigmaForComputingHessianImage = false;
	m_FixedSigmaForHessianImage = 1.0;
	
  m_GenerateScaleOutput = false;
  m_GenerateHessianOutput = false;
	m_GenerateNPlus1DHessianOutput = false;
	m_GenerateNPlus1DHessianMeasureOutput = false;

  this->ProcessObject::SetNumberOfRequiredOutputs(5);
  this->ProcessObject::SetNthOutput(1,this->MakeOutput(1));
  this->ProcessObject::SetNthOutput(2,this->MakeOutput(2));
  this->ProcessObject::SetNthOutput(3,this->MakeOutput(3));	
	this->ProcessObject::SetNthOutput(4,this->MakeOutput(4));	
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::EnlargeOutputRequestedRegion (DataObject *output)
{
  // currently this filter can not stream so we must set the outputs
  // requested region to the largest
  output->SetRequestedRegionToLargestPossibleRegion();
	
	// Also enlarge the requested region of the (N+1)-D image output. Note that, the 
	// requested regions of the remaining outputs will be set when GenerateOutputRequestedRegion() 
	// function of the ProcessObject is called. Unfortunately, for the (N+1)-D image output, this 
	// doesn't work since the dynamic cast in the SetRequestedRegion() of the ImageBase<N> fails.
	// That is why, we have to do it manually here.
	typename OutputNPlus1DImageType::Pointer  outputPtr = 
	dynamic_cast<OutputNPlus1DImageType*>(this->ProcessObject::GetOutput(3));
	outputPtr->SetRequestedRegionToLargestPossibleRegion();
	
	typename NPlus1DHessianImageType::Pointer  nPlus1DHessianPtr = 
	dynamic_cast<NPlus1DHessianImageType*>(this->ProcessObject::GetOutput(4));
	nPlus1DHessianPtr->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
typename MultiScaleHessianBasedMeasureImageFilter2 
  <TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>::DataObjectPointer
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::MakeOutput(unsigned int idx)
{
  if (idx == 1)
	{
    return static_cast<DataObject*>(ScaleImageType::New().GetPointer());
	}
  else if (idx == 2)
	{
    return static_cast<DataObject*>(HessianImageType::New().GetPointer());
	}
	else if (idx == 3)
	{
    return static_cast<DataObject*>(OutputNPlus1DImageType::New().GetPointer());
	}
	else if (idx == 4)
	{
    return static_cast<DataObject*>(NPlus1DHessianImageType::New().GetPointer());
	}	
	else // (idx == 0)
	{
    return static_cast<DataObject*>(OutputNDImageType::New().GetPointer());		
	}
}


template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::CallCopyInputRegionToOutputRegion(OutputNPlus1DRegionType &destRegion,
                                    const InputRegionType &srcRegion)
{
  InputToOutputRegionCopierType regionCopier;
  regionCopier(destRegion, srcRegion);
}

/** Copy information from the input image to output images except 
	* the (N+1)-D hessian-based image. For that, we compute the scales 
	* and create a scale-space image. */
template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::GenerateOutputInformation()
{
	// Copy information for the first three output images first.
  DataObject::ConstPointer input = this->ProcessObject::GetInput(0);
	DataObject::Pointer output;
  if( input )
	{
		for(unsigned int idx = 0; idx < 3; ++idx)
		{
			output = this->ProcessObject::GetOutput(idx);
			if( output )
			{
				output->CopyInformation(input);
			}
		}
	}
	
	// Now, copy information for the (N+1)-D output image.
	typename OutputNPlus1DImageType::Pointer  outputPtr = 
	dynamic_cast<OutputNPlus1DImageType*>(this->ProcessObject::GetOutput(3));
	typename InputImageType::ConstPointer  inputPtr  = this->GetInput();
	
	if ( !outputPtr || !inputPtr)
	{
		return;
	}
	
	// Set the output image largest possible region.  Use a RegionCopier
	// so that the input and output images can be different dimensions.
	OutputNPlus1DRegionType outputLargestPossibleRegion;
  this->CallCopyInputRegionToOutputRegion(outputLargestPossibleRegion,
                                          inputPtr->GetLargestPossibleRegion());
	typename OutputNPlus1DImageType::SizeType regionSize = outputLargestPossibleRegion.GetSize();
	regionSize[OutputNPlus1DImageType::ImageDimension-1] = m_NumberOfSigmaSteps;
	outputLargestPossibleRegion.SetSize( regionSize );
	outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
	
	// Set the output spacing and origin	
	if( inputPtr.IsNotNull() )
	{
		// Copy what we can from the image from spacing and origin of the input.
		unsigned int i, j;
		const typename InputImageType::SpacingType& inputSpacing = inputPtr->GetSpacing();
		const typename InputImageType::PointType& inputOrigin = inputPtr->GetOrigin();
		const typename InputImageType::DirectionType& inputDirection = inputPtr->GetDirection();
		
		typename OutputNPlus1DImageType::SpacingType outputSpacing;
		typename OutputNPlus1DImageType::PointType outputOrigin;
		typename OutputNPlus1DImageType::DirectionType outputDirection;
		
		// copy the input to the output and fill the rest of the
		// output with zeros.
		for( i=0; i < InputImageType::ImageDimension; ++i )
		{
			outputSpacing[i] = inputSpacing[i];
			outputOrigin[i] = inputOrigin[i];
			for( j=0; j < OutputNPlus1DImageType::ImageDimension; j++ )
			{
				if( j < InputImageType::ImageDimension )
				{
					outputDirection[j][i] = inputDirection[j][i];
				}
				else
				{
					outputDirection[j][i] = 0.0;
				}
			}
		}
		for( ; i < OutputNPlus1DImageType::ImageDimension; ++i )
		{
			outputSpacing[i] = vnl_math_max(NumericTraits<double>::epsilon(), 
																			(m_SigmaMaximum - m_SigmaMinimum ) / (m_NumberOfSigmaSteps - 1));
			outputOrigin[i] = m_SigmaMinimum;
			
			for( j=0; j < OutputNPlus1DImageType::ImageDimension; j++ )
			{
				if (j == i)
				{
					outputDirection[j][i] = 1.0;
				}
				else
				{
					outputDirection[j][i] = 0.0;
				}
			}
		}
		
		// set the spacing and origin
		outputPtr->SetSpacing( outputSpacing );
		outputPtr->SetOrigin( outputOrigin );
		outputPtr->SetDirection( outputDirection );
		outputPtr->SetNumberOfComponentsPerPixel( // propagate vector length info
																						 inputPtr->GetNumberOfComponentsPerPixel());
																						 
		// Now, do the same thing for the (N+1)-D Hessian image.
		typename NPlus1DHessianImageType::Pointer  outputNPlus1DHessianPtr = 
		dynamic_cast<NPlus1DHessianImageType*>(this->ProcessObject::GetOutput(4));
		if ( !outputNPlus1DHessianPtr )
		{
			return;
		}
		outputNPlus1DHessianPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
		outputNPlus1DHessianPtr->SetSpacing( outputSpacing );
		outputNPlus1DHessianPtr->SetOrigin( outputOrigin );
		outputNPlus1DHessianPtr->SetDirection( outputDirection );
		outputNPlus1DHessianPtr->SetNumberOfComponentsPerPixel( // propagate vector length info
																						 inputPtr->GetNumberOfComponentsPerPixel());
	}
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::AllocateUpdateBuffer()
{
  /* The update buffer looks just like the output and holds the best response
     in the  objectness measure */

  typename TOutputNDImage::Pointer output = this->GetOutput();

  // this copies meta data describing the output such as origin,
  // spacing and the largest region
  m_UpdateBuffer->CopyInformation(output);

  m_UpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_UpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_UpdateBuffer->Allocate();

	m_UpdateBuffer->FillBuffer( itk::NumericTraits< BufferValueType >::NonpositiveMin() );
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::GenerateData()
{
	if( m_HessianToMeasureFilter.IsNull() )
	{
		itkExceptionMacro( " HessianToMeasure filter is not set. Use SetHessianToMeasureFilter() " );
	}

  // TODO: Move the allocation to a derived AllocateOutputs method
  // Allocate the output
  this->GetOutput()->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
  this->GetOutput()->Allocate();

  if (m_GenerateScaleOutput)
	{
    typename ScaleImageType::Pointer scaleImage = 
		dynamic_cast<ScaleImageType*>(this->ProcessObject::GetOutput(1));
		
    scaleImage->SetBufferedRegion(scaleImage->GetRequestedRegion());
    scaleImage->Allocate();
    scaleImage->FillBuffer(0);
	}
	
  if (m_GenerateHessianOutput)
	{
    typename HessianImageType::Pointer hessianImage = 
		dynamic_cast<HessianImageType*>(this->ProcessObject::GetOutput(2));
		
    hessianImage->SetBufferedRegion(hessianImage->GetRequestedRegion());
    hessianImage->Allocate();
    // SymmetricSecondRankTensor is already filled with zero elements at construction. 
    // No strict need of filling the buffer, but we do it explicitly here to make sure.
    typename HessianImageType::PixelType zeroTensor(0.0);
    hessianImage->FillBuffer(zeroTensor);
	}
	
	if (m_GenerateNPlus1DHessianOutput)
	{
    typename NPlus1DHessianImageType::Pointer nPlus1DHessianImage = 
		dynamic_cast<NPlus1DHessianImageType*>(this->ProcessObject::GetOutput(4));
		
    nPlus1DHessianImage->SetBufferedRegion(nPlus1DHessianImage->GetRequestedRegion());
    nPlus1DHessianImage->Allocate();
	}
	
	if( m_GenerateNPlus1DHessianMeasureOutput )
	{
		typename OutputNPlus1DImageType::Pointer outputNPlus1DImage = 
		dynamic_cast<OutputNPlus1DImageType*>(this->ProcessObject::GetOutput(3));
		
    outputNPlus1DImage->SetBufferedRegion(outputNPlus1DImage->GetRequestedRegion());
    outputNPlus1DImage->Allocate();
	}

  // Allocate the buffer
  AllocateUpdateBuffer();

  typename InputImageType::ConstPointer input = this->GetInput();

  this->m_HessianFilter->SetInput(input);

	// This is the Lindeberg-style normalization of Gaussian derivatives.
	// Here, we don't need it.
	// For more info, esee the documentation of RecursiveGaussianImageFilter class. 
  this->m_HessianFilter->SetNormalizeAcrossScale(false);

  unsigned int scaleLevel = 0;
  double sigma = m_SigmaMinimum;

  while( (sigma - m_SigmaMaximum) < NumericTraits<double>::epsilon() )
	{
    if ( m_NumberOfSigmaSteps == 0 )
		{
      break;
		}

    itkDebugMacro ( << "Computing measure for scale with sigma = " << sigma );

		if( !m_UseAFixedSigmaForComputingHessianImage )
		{
			m_HessianFilter->SetSigma( sigma );
		}
		else
		{
			m_HessianFilter->SetSigma( m_FixedSigmaForHessianImage );			
		}

    m_HessianToMeasureFilter->SetInput ( m_HessianFilter->GetOutput() );
	
		CallHessianToMeasureFilterSigmaMethod(m_HessianToMeasureFilter, sigma, Int2Type<
		DoesHessianToMeasureFilterHaveSetSigmaMethod<HessianToMeasureFilterType>::Has>());

    m_HessianToMeasureFilter->Update();

    this->UpdateMaximumResponse(sigma, scaleLevel);
		
		scaleLevel++;
    sigma  = this->ComputeSigmaValue( scaleLevel );

    if ( m_NumberOfSigmaSteps == 1 )
		{
      break;
		}
	}

  // Write out the best response to the output image
  // we can assume that the meta-data should match between these two
  // images, therefore we iterate over the desired output region
  OutputNDRegionType outputRegion = this->GetOutput()->GetBufferedRegion();
  ImageRegionIterator<UpdateBufferType> it( m_UpdateBuffer, outputRegion );
  it.GoToBegin();

  ImageRegionIterator<TOutputNDImage> oit( this->GetOutput(), outputRegion );
  oit.GoToBegin();

  while(!oit.IsAtEnd())
	{
    oit.Value() = static_cast< OutputNDPixelType >( it.Get() );
    ++oit;
    ++it;
	}

  // Release data from the update buffer.
  m_UpdateBuffer->ReleaseData();
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::UpdateMaximumResponse(double sigma, unsigned int scaleLevel)
{
  // the meta-data should match between these images, therefore we
  // iterate over the desired output region 
  OutputNDRegionType outputRegion = this->GetOutput()->GetBufferedRegion();
  ImageRegionIterator<UpdateBufferType> oit( m_UpdateBuffer, outputRegion );

  typename ScaleImageType::Pointer scaleImage = static_cast<ScaleImageType*>(this->ProcessObject::GetOutput(1));
  ImageRegionIterator<ScaleImageType> osit;

  typename HessianImageType::Pointer hessianImage = static_cast<HessianImageType*>(this->ProcessObject::GetOutput(2));
  ImageRegionIterator<HessianImageType> ohit;
		
	typename OutputNPlus1DImageType::Pointer outputNPlus1DImage = 
	static_cast<OutputNPlus1DImageType*>(this->ProcessObject::GetOutput(3));
	ImageRegionIterator<OutputNPlus1DImageType> o2it;
	
	typename NPlus1DHessianImageType::Pointer nPlus1DHessianImage = 
	static_cast<NPlus1DHessianImageType*>(this->ProcessObject::GetOutput(4));
	ImageRegionIterator<NPlus1DHessianImageType> o3it;	

  oit.GoToBegin();
  if( m_GenerateScaleOutput )
	{
    osit = ImageRegionIterator<ScaleImageType> ( scaleImage, outputRegion );
    osit.GoToBegin();
	}
  if( m_GenerateHessianOutput )
	{
    ohit = ImageRegionIterator<HessianImageType> ( hessianImage, outputRegion );
    ohit.GoToBegin();
	}
	
	// Create an iterator for the given sigma layer of the (N+1)-D image.
	OutputNPlus1DRegionType outputNPlus1DRegion;
	if( m_GenerateNPlus1DHessianMeasureOutput || m_GenerateNPlus1DHessianOutput )
	{
	  this->CallCopyInputRegionToOutputRegion(outputNPlus1DRegion, outputRegion);
		typename OutputNPlus1DImageType::IndexType outputNPlus1DRegionIndex = outputNPlus1DRegion.GetIndex();
		typename OutputNPlus1DImageType::SizeType outputNPlus1DRegionSize = outputNPlus1DRegion.GetSize();
		outputNPlus1DRegionIndex[OutputNPlus1DImageType::ImageDimension-1] = scaleLevel;
		outputNPlus1DRegionSize[OutputNPlus1DImageType::ImageDimension-1] = 1;
		outputNPlus1DRegion.SetIndex( outputNPlus1DRegionIndex );
		outputNPlus1DRegion.SetSize( outputNPlus1DRegionSize );
	}
	
	if( m_GenerateNPlus1DHessianMeasureOutput )
	{
    o2it = ImageRegionIterator<OutputNPlus1DImageType> ( outputNPlus1DImage, outputNPlus1DRegion );
    o2it.GoToBegin();
	}
	
	if( m_GenerateNPlus1DHessianOutput )
	{
		o3it = ImageRegionIterator<NPlus1DHessianImageType> ( nPlus1DHessianImage, outputNPlus1DRegion );
		o3it.GoToBegin();
	}
	

  typedef typename HessianToMeasureFilterType::OutputImageType HessianToMeasureOutputImageType;

  ImageRegionIterator<HessianToMeasureOutputImageType> it( m_HessianToMeasureFilter->GetOutput(), outputRegion );
  ImageRegionIterator<HessianImageType> hit( m_HessianFilter->GetOutput(), outputRegion );

  it.GoToBegin();
  hit.GoToBegin();

  while(!oit.IsAtEnd())
	{
    if( oit.Value() < it.Value() )
		{
      oit.Value() = it.Value();
      if( m_GenerateScaleOutput )
			{
        osit.Value() = static_cast< ScalePixelType >( sigma );
			}
      if( m_GenerateHessianOutput )
			{
        ohit.Value() = hit.Value();
			}
		}
		if( m_GenerateNPlus1DHessianMeasureOutput )
		{
			o2it.Value() = it.Value();
			++o2it;
		}
		if( m_GenerateNPlus1DHessianOutput )
		{
			o3it.Value() = hit.Value();
			++o3it;
		}
    ++oit;
    ++it;
    if( m_GenerateScaleOutput )
		{
      ++osit;
		}
    if( m_GenerateHessianOutput )
		{
      ++ohit;
		}
		if( m_GenerateHessianOutput || m_GenerateNPlus1DHessianOutput )
		{
      ++hit;
		}
	}
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
double
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::ComputeSigmaValue(int scaleLevel)
{
  if (m_NumberOfSigmaSteps < 2)
	{
    return m_SigmaMinimum;
	}

	const double stepSize = vnl_math_max(1e-10, ( m_SigmaMaximum - m_SigmaMinimum ) / (m_NumberOfSigmaSteps - 1));
  return m_SigmaMinimum + stepSize * scaleLevel;
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
typename MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>::OutputNPlus1DImageType * 
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::GetNPlus1DImageOutput()
{
  return static_cast<const OutputNPlus1DImageType*>(this->ProcessObject::GetOutput(3));
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
typename MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>::HessianImageType * 
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::GetHessianOutput()
{
  return static_cast<HessianImageType*>(this->ProcessObject::GetOutput(2));
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
typename MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>::NPlus1DHessianImageType * 
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::GetNPlus1DHessianOutput()
{
  return static_cast<NPlus1DHessianImageType*>(this->ProcessObject::GetOutput(4));
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
typename MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>::ScaleImageType * 
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::GetScaleOutput()
{
  return static_cast<const ScaleImageType*>(this->ProcessObject::GetOutput(1));
}

template <typename TInputImage,
          typename THessianImage,
          typename TScaleImage,
          typename THessianToMeasureFilter,
          typename TOutputNDImage>
void
MultiScaleHessianBasedMeasureImageFilter2
<TInputImage,THessianImage,TScaleImage,THessianToMeasureFilter,TOutputNDImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "SigmaMinimum:  " << m_SigmaMinimum << std::endl;
  os << indent << "SigmaMaximum:  " << m_SigmaMaximum  << std::endl;
  os << indent << "NumberOfSigmaSteps:  " << m_NumberOfSigmaSteps  << std::endl;
  os << indent << "HessianToMeasureFilter: " << m_HessianToMeasureFilter << std::endl;
  os << indent << "GenerateScaleOutput: " << m_GenerateScaleOutput << std::endl;
  os << indent << "GenerateHessianOutput: " << m_GenerateHessianOutput << std::endl;
	os << indent << "GenerateNPlus1DHessianMeasureOutput: " << m_GenerateNPlus1DHessianMeasureOutput << std::endl;
	os << indent << "GenerateNPlus1DHessianOutput: " << m_GenerateNPlus1DHessianOutput << std::endl;
}


} // end namespace itk

#endif
