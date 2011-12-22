/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilter_txx
#define __itkFFTConvolutionImageFilter_txx

#include "itkFFTConvolutionImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkMultiplyImageFilter.h"

namespace itk {

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::GenerateData()
{
  ComplexImagePointerType input;
  ComplexImagePointerType kernel;
  
  this->PrepareInputs( input, kernel, 0.6f );

  typedef itk::MultiplyImageFilter< ComplexImageType,
                                    ComplexImageType,
                                    ComplexImageType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, input );
  mult->SetInput( 1, kernel );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  this->RegisterInternalFilter( mult, 0.1f );

  this->ProduceOutput( mult->GetOutput(), 0.3f );
}

}// end namespace itk
#endif
