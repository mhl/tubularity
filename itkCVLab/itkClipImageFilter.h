/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkClipImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2010-04-18 16:14:24 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkClipImageFilter_h
#define __itkClipImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
  
/** \class ClipImageFilter
 *
 * \brief Clips input pixels to output pixel type.
 *
 * This filter is templated over the input image type
 * and the output image type.
 * 
 * \author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 * \sa CastImageFilter
 */
namespace Functor {  
  
template< class TInput, class TOutput>
class Clip
{
public:
  Clip() {};
  virtual ~Clip() {};
  bool operator!=( const Clip & ) const
    {
    return false;
    }
  bool operator==( const Clip & other ) const
    {
    return !(*this != other);
    }
  inline TOutput operator()( const TInput & A ) const
    {
    if( (double)A > (double)NumericTraits<TOutput>::max() )
      {
      return NumericTraits<TOutput>::max();
      }
    if( (double)A < (double)NumericTraits<TOutput>::NonpositiveMin() )
      {
      return NumericTraits<TOutput>::NonpositiveMin();
      }
    return static_cast<TOutput>( A );
    }
};
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ClipImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::Clip< 
  typename TInputImage::PixelType, 
  typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef ClipImageFilter               Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::Clip< 
                          typename TInputImage::PixelType, 
                          typename TOutputImage::PixelType>   
                                     >  Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ClipImageFilter, UnaryFunctorImageFilter);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<typename TInputImage::PixelType,
                          typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  ClipImageFilter() {}
  virtual ~ClipImageFilter() {}

  void GenerateData()
    {
    if( this->GetInPlace() && this->CanRunInPlace() )
      {
      // nothing to do, so avoid iterating over all the pixels
      // for nothing! Allocate the output, generate a fake progress and exit
      this->AllocateOutputs();
      ProgressReporter progress(this, 0, 1);
      return;
      }
    Superclass::GenerateData();
    }
  

  
private:
  ClipImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


} // end namespace itk

#endif
