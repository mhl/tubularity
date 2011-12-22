/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNormalizeToConstantImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNormalizeToConstantImageFilter_txx
#define __itkNormalizeToConstantImageFilter_txx

#include "itkNormalizeToConstantImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkNumericTraits.h"
#include "itkStatisticsImageFilter.h"
#include "itkDivideByConstantImageFilter.h"

namespace itk {

template <class TInputImage, class TOutputImage>
NormalizeToConstantImageFilter<TInputImage, TOutputImage>
::NormalizeToConstantImageFilter()
{
  m_Constant = NumericTraits<RealType>::One;
}

template <class TInputImage, class TOutputImage>
void 
NormalizeToConstantImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput(0));
  if ( !input0 )
    { 
    return;
    }
  
  input0->SetRequestedRegion( input0->GetLargestPossibleRegion() );
}


template<class TInputImage, class TOutputImage>
void
NormalizeToConstantImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();
  const InputImageType * input0 = this->GetInput(0);
  OutputImageType * output0 = this->GetOutput(0);

  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef typename itk::StatisticsImageFilter< InputImageType > StatType;  
  typename StatType::Pointer stat = StatType::New();
  stat->SetInput( input0 );
  progress->RegisterInternalFilter( stat, 0.5f );
  stat->SetNumberOfThreads( this->GetNumberOfThreads() );
  stat->Update();
  
  typedef typename itk::DivideByConstantImageFilter< InputImageType, double, OutputImageType > DivType;
  typename DivType::Pointer div = DivType::New();
  div->SetInput( input0 );
  div->SetConstant( stat->GetSum() / m_Constant );
  progress->RegisterInternalFilter( div, 0.5f );
  div->SetNumberOfThreads( this->GetNumberOfThreads() );

  div->GraftOutput( output0 );
  div->Update();
  this->GraftOutput( div->GetOutput() );
}


template<class TInputImage, class TOutputImage>
void
NormalizeToConstantImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Constant: "  << m_Constant << std::endl;
}
  
}// end namespace itk
#endif
