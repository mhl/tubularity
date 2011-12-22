/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTWLock.h,v $
  Language:  C++
  Date:      $Date: 2008-12-21 19:13:11 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTWLock_h
#define __itkFFTWLock_h

#include "itkSimpleFastMutexLock.h"
#if defined(USE_FFTWF) || defined(USE_FFTWD)
#include "fftw3.h"
#endif

namespace itk
{
/**
 * A simple global lock which must be called before calling FFTW unsafe functions.
 * It also handle cleanly the initialization and cleanup of FFTW.
 */
class  FFTWLock
{

  public:
  
  static void Lock();
  static void Unlock();
  static void NewWisdomAvailable();
  
  enum { ESTIMATE=FFTW_ESTIMATE,
         PATIENT=FFTW_PATIENT,
         EXHAUSTIVE=FFTW_EXHAUSTIVE } Optimization;
         
  static int GetGlobalOptimizationLevel();
  static void SetGlobalOptimizationLevel(int opt);
  
  static std::string GetWisdomFileDefaultBaseName();
  // static void SetWisdomFileDefaultBaseName( std::string fname );
  
  bool ImportWisdomFileDouble( std::string fname );
  bool ExportWisdomFileDouble( std::string fname );
  bool ImportWisdomFileDouble();
  bool ExportWisdomFileDouble();

  bool ImportWisdomFileFloat( std::string fname );
  bool ExportWisdomFileFloat( std::string fname );
  bool ImportWisdomFileFloat();
  bool ExportWisdomFileFloat();
  
  protected:
  
  FFTWLock();
  ~FFTWLock();
  
  static FFTWLock            m_Singleton;
  static SimpleFastMutexLock m_Lock;
  static bool                m_NewWisdomAvailable;
  static int                 m_GlobalOptimizationLevel;
  // static std::string         m_WisdomFileDefaultBaseName;

};
}
#endif
