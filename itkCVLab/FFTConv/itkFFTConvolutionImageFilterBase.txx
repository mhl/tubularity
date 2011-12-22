/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilterBase.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilterBase_txx
#define __itkFFTConvolutionImageFilterBase_txx

#include "itkFFTConvolutionImageFilterBase.h"
#include "itkFlipImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"
#include "itkClipImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkChangeInformationImageFilter.h"

namespace itk {

template <class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::FFTConvolutionImageFilterBase()
{
  m_Normalize = true;
  m_CropOutput = true;
  m_PaddedRegion = RegionType();
  this->SetNumberOfRequiredInputs(2);
  // don't use the second output of the superclass
  this->SetNumberOfRequiredOutputs(1);
  this->SetNthOutput( 1, NULL );
}

template <class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void 
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput());
  InputImageType * input1 = const_cast<KernelImageType *>(this->GetKernelImage());
  if ( !input0 || !input1 )
    { 
    return;
    }
  input0->SetRequestedRegion( input0->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( input1->GetLargestPossibleRegion() );
}

template <class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void 
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::GenerateOutputInformation()
{
  // call the implementation from FFTPad if the output shouldn't be cropped, and the one
  // from ImageToImageFilter if the output should be cropped
  if( m_CropOutput )
    {
    Superclass::Superclass::GenerateOutputInformation();
    }
  else
    {
    Superclass::GenerateOutputInformation();
    }
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrepareInputs( InternalImagePointerType & paddedInput, InternalImagePointerType & paddedKernel, float progressWeight )
{
  const InputImageType * input = this->GetInput();
  const KernelImageType * kernel = this->GetKernelImage();

  // Create a process accumulator for tracking the progress of this minipipeline
  m_ProgressAccumulator = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  typedef itk::FlipImageFilter< KernelImageType > FlipType;
  typename FlipType::Pointer flip = FlipType::New();
  flip->SetInput( kernel );
  typename FlipType::FlipAxesArrayType axes;
  axes.Fill( true ); // we must flip all the axes
  flip->SetFlipAxes( axes );
  flip->SetNumberOfThreads( this->GetNumberOfThreads() );
  flip->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( flip, 0.33 * progressWeight );
    }

  typedef itk::ImageToImageFilter< KernelImageType, InternalImageType > NormType;
  typedef itk::NormalizeToConstantImageFilter< KernelImageType, InternalImageType > NormConstType;
  typedef itk::CastImageFilter< KernelImageType, InternalImageType > CastType;
  typename NormType::Pointer norm;
  if( m_Normalize )
    {
    norm = NormConstType::New();
    }
  else
    {
    norm = CastType::New();
    }
  norm->SetInput( flip->GetOutput() );
  norm->SetNumberOfThreads( this->GetNumberOfThreads() );
  norm->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( norm,  0.33 * progressWeight );
    }

  typedef itk::FFTPadImageFilter< InputImageType, InternalImageType, InternalImageType, InternalImageType > PadType;
  typename PadType::Pointer pad = PadType::New();
  pad->SetInput( input );
  pad->SetInputKernel( norm->GetOutput() );
  pad->SetNumberOfThreads( this->GetNumberOfThreads() );
  pad->SetReleaseDataFlag( true );
  pad->SetPadMethod( this->GetPadMethod() );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( pad,  0.34 * progressWeight );
    }
  // instantiate a fft filter to know if it is a vnl implementation
  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  // vnl filters need a size which is a power of 2
  if( std::string(fft->GetNameOfClass()).find("Vnl") == 0 )
    {
    pad->SetPadToPowerOfTwo( true );
    // fake the cyclic behavior
    if( this->GetPadMethod() == this->NO_PADDING )
      {
      pad->SetPadMethod( this->WRAP );
      }
    }
  else
    {
    pad->SetGreatestPrimeFactor( this->GetGreatestPrimeFactor() );
    }
  
  pad->Update();
  
  paddedInput = pad->GetOutput();
  paddedInput->DisconnectPipeline();
  paddedKernel = pad->GetOutputKernel();
  paddedKernel->DisconnectPipeline();
  m_PaddedRegion = paddedInput->GetLargestPossibleRegion();
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrepareInputs( InternalImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, float progressWeight )
{
  InternalImagePointerType pk;
  
  this->PrepareInputs( paddedInput, pk, 0.15 * progressWeight );

  typedef itk::FFTShiftImageFilter< InternalImageType, InternalImageType > ShiftType;
  typename ShiftType::Pointer shift = ShiftType::New();
  shift->SetInput( pk );
  shift->SetInverse( true );
  shift->SetNumberOfThreads( this->GetNumberOfThreads() );
  shift->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( shift, 0.05 * progressWeight );
    }
  
  typename FFTFilterType::Pointer fftk = FFTFilterType::New();
  fftk->SetInput( shift->GetOutput() );
  fftk->SetNumberOfThreads( this->GetNumberOfThreads() );
  fftk->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( fftk, 0.8 * progressWeight );
    }
  
  fftk->Update();
  paddedKernel = fftk->GetOutput();
  paddedKernel->DisconnectPipeline();
  
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrepareInputs( ComplexImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, float progressWeight )
{
  InternalImagePointerType pi;
  
  this->PrepareInputs( pi, paddedKernel, progressWeight * 0.6 );
  
  RegionType region = pi->GetLargestPossibleRegion();

  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  fft->SetInput( pi );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( fft, progressWeight * 0.4 );
    }

  fft->Update();
  paddedInput = fft->GetOutput();
  paddedInput->DisconnectPipeline();
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::RegisterInternalFilter( ProcessObject * filter, float progressWeight )
{
  m_ProgressAccumulator->RegisterInternalFilter( filter, progressWeight );
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
bool
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::GetXIsOdd() const
{
  return m_PaddedRegion.GetSize()[0] % 2;
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::ProduceOutput( ComplexImageType * paddedOutput, float progressWeight )
{
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  typename IFFTFilterType::Pointer ifft = IFFTFilterType::New();
  ifft->SetInput( paddedOutput );
  ifft->SetActualXDimensionIsOdd( GetXIsOdd() );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( ifft, progressWeight * 0.8 );
    }
  
  this->ProduceOutput( ifft->GetOutput(), progressWeight * 0.2 );
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::ProduceOutput( InternalImageType * paddedOutput, float progressWeight )
{
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  typedef itk::ClipImageFilter< InternalImageType, OutputImageType > ClipType;
  typename ClipType::Pointer clip = ClipType::New();
  clip->SetInput( paddedOutput );
  clip->SetNumberOfThreads( this->GetNumberOfThreads() );
  clip->SetReleaseDataFlag( true );
  clip->SetInPlace( true );
  
  if( !m_CropOutput )
    {
    if( progressWeight != 0 )
      {
      m_ProgressAccumulator->RegisterInternalFilter( clip, 0.2 * progressWeight );
      }
    this->AllocateOutputs();
    clip->GraftOutput( this->GetOutput() );
    clip->Update();
    this->GraftOutput( clip->GetOutput() );  
    }
  else
    {
    if( progressWeight != 0 )
      {
      m_ProgressAccumulator->RegisterInternalFilter( clip, 0.1 * progressWeight );
      }
  
    typedef itk::RegionFromReferenceImageFilter< OutputImageType, OutputImageType > CropType;
    typename CropType::Pointer crop = CropType::New();
    crop->SetInput( clip->GetOutput() );
    crop->SetReferenceImage( this->GetInput() );
    crop->SetNumberOfThreads( this->GetNumberOfThreads() );
    if( progressWeight != 0 )
      {
      m_ProgressAccumulator->RegisterInternalFilter( crop, 0.1 * progressWeight );
      }
    this->AllocateOutputs();
    crop->GraftOutput( this->GetOutput() );
    crop->Update();
    this->GraftOutput( crop->GetOutput() );
    }
  m_ProgressAccumulator = NULL;
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
typename FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>::InternalImageConstPointerType
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::InternalInput( int pos ) const
{
  // we have 3 kind of possible input image type - lets try them all
  const InputImageType * inputImage = dynamic_cast<const InputImageType *>( this->ProcessObject::GetInput( pos ) );
  const KernelImageType * kernelImage = dynamic_cast<const KernelImageType *>( this->ProcessObject::GetInput( pos ) );
  const InternalImageType * internalImage = dynamic_cast<const InternalImageType *>( this->ProcessObject::GetInput( pos ) );

  typedef CastImageFilter< InputImageType, InternalImageType > InputCastType;
  typename InputCastType::Pointer inputCast = InputCastType::New();

  typedef CastImageFilter< KernelImageType, InternalImageType > KernelCastType;
  typename KernelCastType::Pointer kernelCast = KernelCastType::New();

  if( inputImage != NULL )
    {
    inputCast->SetInput( inputImage );
    inputCast->SetNumberOfThreads( this->GetNumberOfThreads() );
    inputCast->SetReleaseDataFlag( true );
    inputCast->SetInPlace( false );
    inputCast->Update();
    internalImage = inputCast->GetOutput();
    }
  else if( kernelImage != NULL )
    {
    kernelCast->SetInput( kernelImage );
    kernelCast->SetNumberOfThreads( this->GetNumberOfThreads() );
    kernelCast->SetReleaseDataFlag( true );
    kernelCast->SetInPlace( false );
    kernelCast->Update();
    internalImage = kernelCast->GetOutput();
    }
  else if( internalImage != NULL )
    {
    // just do nothing - it's already of the right type
    }
  else
    {
    // or should it send an exception?
    return NULL;
    }
  return internalImage;
  }

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrepareImage( InternalImagePointerType & paddedInput, const InternalImageType * img, bool flipImage, bool normalizeImage, float progressWeight )
{
  if( img == NULL )
    {
    paddedInput = NULL;
    return;
    }
  
  int nbOfFilters = 1; if( flipImage ) nbOfFilters++; if( normalizeImage ) nbOfFilters++;
  const InternalImageType * internalImage = img;
    
  typedef itk::FlipImageFilter< InternalImageType > FlipType;
  typename FlipType::Pointer flip = FlipType::New();
  if( flipImage )
    {
    flip->SetInput( internalImage );
    typename FlipType::FlipAxesArrayType axes;
    axes.Fill( true ); // we must flip all the axes
    flip->SetFlipAxes( axes );
    flip->SetNumberOfThreads( this->GetNumberOfThreads() );
    flip->SetReleaseDataFlag( true );
    if( progressWeight != 0 )
      {
      m_ProgressAccumulator->RegisterInternalFilter( flip, 1.0 / nbOfFilters * progressWeight );
      }
    internalImage = flip->GetOutput();
    }

  typedef itk::NormalizeToConstantImageFilter< InternalImageType, InternalImageType > NormConstType;
  typename NormConstType::Pointer norm;
  if( normalizeImage )
    {
    norm = NormConstType::New();
    norm->SetInput( internalImage );
    norm->SetNumberOfThreads( this->GetNumberOfThreads() );
    norm->SetReleaseDataFlag( true );
    if( progressWeight != 0 )
      {
      m_ProgressAccumulator->RegisterInternalFilter( norm, 1.0 / nbOfFilters * progressWeight );
      }
    internalImage = norm->GetOutput();
    }

  RegionType referenceRegion = this->GetPaddedRegion();
  RegionType inputRegion = img->GetLargestPossibleRegion();
  SizeType s;

  typedef typename itk::ConstantPadImageFilter< InternalImageType, InternalImageType > KernelPadType;
  typename KernelPadType::Pointer pad = KernelPadType::New();
  pad->SetInput( internalImage );
  pad->SetNumberOfThreads( this->GetNumberOfThreads() );
  for( int i=0; i<ImageDimension; i++ )
    {
    s[i] = ( referenceRegion.GetSize()[i] - inputRegion.GetSize()[i] ) / 2;
    }
  pad->SetPadUpperBound( s );
  for( int i=0; i<ImageDimension; i++ )
    {
    // s[i] = itk::Math::Ceil(( referenceRegion.GetSize()[i] - inputRegion.GetSize()[i] ) / 2.0 );
    // this line should do the same, but without requirement on ITK cvs
    s[i] = ( referenceRegion.GetSize()[i] - inputRegion.GetSize()[i] ) / 2 +  ( referenceRegion.GetSize()[i] - inputRegion.GetSize()[i] ) % 2;
    }
  pad->SetPadLowerBound( s );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( pad, 1.0 / nbOfFilters * progressWeight );
    }
  
  typedef typename itk::ChangeInformationImageFilter< InternalImageType > ChangeType;
  typename ChangeType::Pointer change = ChangeType::New();
  change->SetInput( pad->GetOutput() );
  // can't be used because the reference is not of class ImageBase
  // change->SetUseReferenceImage( true );
  // change->SetReferenceImage( reference );
  change->SetOutputOffset( const_cast<long *>((referenceRegion.GetIndex() - inputRegion.GetIndex() + s).m_Offset) );
  change->SetChangeRegion( true );
  // no progress for change - it does almost nothing
  
  paddedInput = change->GetOutput();
  paddedInput->Update();
  paddedInput->DisconnectPipeline();
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrepareImage( ComplexImagePointerType & paddedInput, const InternalImageType * img, bool flipImage, bool normalizeImage, float progressWeight )
{
  if( img == NULL )
    {
    paddedInput = NULL;
    return;
    }

  InternalImagePointerType pi;
  
  this->PrepareImage( pi, img, flipImage, normalizeImage, progressWeight * 0.1 );

  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  fft->SetInput( pi );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( fft, progressWeight * 0.9 );
    }

  fft->Update();
  paddedInput = fft->GetOutput();
  paddedInput->DisconnectPipeline();

}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Normalize: "  << m_Normalize << std::endl;
  os << indent << "CropOutput: "  << m_CropOutput << std::endl;
}

}// end namespace itk
#endif
