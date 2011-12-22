/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilterBase.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilterBase_h
#define __itkFFTConvolutionImageFilterBase_h

#include "itkFFTPadImageFilter.h"
#include "itkConceptChecking.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkProgressAccumulator.h"

namespace itk {

/** \class FFTConvolutionImageFilterBase
 * \brief Pad two images with zeros to make them suitable for a convolution in the
 * frequency domain.
 *
 * The requires two input images. The first one is supposed to be the image to convolve
 * and the second one the kernel.
 *
 * The filter pad with zeros the first image by the size of the second image. The second image
 * is centered and padded with zeros to have the same size than the first output image.
 * After the transform, both images have the same LargestPossibleRegion.
 *
 * The option PadToPowerOfTwo can be set to true to force the size of the
 * images be a power of two - if the size of the padded image is 274 without this option, it
 * would be increased to 512 when PadToPowerOfTwo is true.
 * This option is makes the images usable with vnl's implementation of FFT.
 * PadToPowerOfTwo is false by default.
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TKernelImage=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT FFTConvolutionImageFilterBase : 
    public FFTPadImageFilter<TInputImage, TKernelImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FFTConvolutionImageFilterBase Self;

  typedef FFTPadImageFilter<TInputImage, TKernelImage, TOutputImage> Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TKernelImage                             KernelImageType;
  typedef KernelImageType                          PointSpreadFunctionType;
  typedef TOutputImage                             OutputImageType;
  typedef TInternalPrecision                       InternalPrecisionType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename KernelImageType::Pointer        KernelImagePointer;
  typedef typename KernelImageType::ConstPointer   KernelImageConstPointer;
  typedef typename KernelImageType::PixelType      KernelImagePixelType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename InputImageType::RegionType      RegionType;
  typedef typename InputImageType::IndexType       IndexType;
  typedef typename InputImageType::SizeType        SizeType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef ImageBase< ImageDimension >              BaseImageType;
  
  typedef itk::FFTRealToComplexConjugateImageFilter< InternalPrecisionType, ImageDimension > FFTFilterType;
  typedef itk::FFTComplexConjugateToRealImageFilter< InternalPrecisionType, ImageDimension > IFFTFilterType;
  typedef typename FFTFilterType::OutputImageType                                            ComplexImageType;
  typedef typename ComplexImageType::PixelType                                               ComplexType;
  typedef typename ComplexImageType::Pointer                                                 ComplexImagePointerType;

  typedef typename itk::Image< InternalPrecisionType, ImageDimension > InternalImageType;
  typedef typename InternalImageType::Pointer                          InternalImagePointerType;
  typedef typename InternalImageType::ConstPointer                     InternalImageConstPointerType;

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(FFTConvolutionImageFilterBase, FFTPadImageFilter);

  /** Set the kernel image */
  void SetKernelImage(const TKernelImage *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<TKernelImage *>(input) );
    }

  /** Get the kernel image */
  const KernelImageType * GetKernelImage() const
    {
    return static_cast<KernelImageType*>(
      const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

  /** Set the kernel image */
  void SetPointSpreadFunction(const PointSpreadFunctionType *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<PointSpreadFunctionType *>(input) );
    }

  /** Get the kernel image */
  const PointSpreadFunctionType * GetPointSpreadFunction() const
    {
    return static_cast<PointSpreadFunctionType*>(
      const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

  /**
   * Set/Get whether the kernel should be normalized to one or not.
   * Default is true.
   */
  itkGetConstMacro(Normalize, bool);
  itkSetMacro(Normalize, bool);
  itkBooleanMacro(Normalize);
  
  /**
   * Set/Get whether the output image should be cropped to the original input size or not.
   * Default is true.
   */
  itkGetConstMacro(CropOutput, bool);
  itkSetMacro(CropOutput, bool);
  itkBooleanMacro(CropOutput);
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  FFTConvolutionImageFilterBase();
  ~FFTConvolutionImageFilterBase() {};

  void GenerateInputRequestedRegion();
  void GenerateOutputInformation();
  
  void PrepareInputs( InternalImagePointerType & paddedInput, InternalImagePointerType & paddedKernel, float progressWeight=0.66 );
  void PrepareInputs( InternalImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, float progressWeight=0.66 );
  void PrepareInputs( ComplexImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, float progressWeight=0.66 );
  
  InternalImageConstPointerType InternalInput( int pos ) const;
  void PrepareImage( InternalImagePointerType & paddedInput, const InternalImageType * img,
                     bool flip=false, bool normalize=false, float progressWeight=0.66 );
  void PrepareImage( ComplexImagePointerType & paddedInput, const InternalImageType * img,
                     bool flip=false, bool normalize=false, float progressWeight=0.66 );
  
  void RegisterInternalFilter( ProcessObject * filter, float progressWeight );
  
  itkGetConstMacro(PaddedRegion, RegionType);
  itkSetMacro(PaddedRegion, RegionType);
  bool GetXIsOdd() const;
  
  void ProduceOutput( ComplexImageType * paddedOutput, float progressWeight=0.34 );
  void ProduceOutput( InternalImageType * paddedOutput, float progressWeight=0.34 );

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  FFTConvolutionImageFilterBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool                         m_Normalize;
  bool                         m_CropOutput;
  ProgressAccumulator::Pointer m_ProgressAccumulator;
  RegionType                   m_PaddedRegion;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTConvolutionImageFilterBase.txx"
#endif

#endif
