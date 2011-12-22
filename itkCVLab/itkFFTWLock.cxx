/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFastMutexLock.cxx,v $
  Language:  C++
  Date:      $Date: 2003-09-10 14:29:07 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkFFTWLock.h"
#include "itksys/SystemTools.hxx"
#include <sys/file.h>

namespace itk
{
SimpleFastMutexLock    FFTWLock::m_Lock;
FFTWLock               FFTWLock::m_Singleton;
bool                   FFTWLock::m_NewWisdomAvailable;
int                    FFTWLock::m_GlobalOptimizationLevel;
// std::string            FFTWLock::m_WisdomFileDefaultBaseName;

FFTWLock
::FFTWLock()
{ 
  // std::cout << "======== init fftw stuff =========" << std::endl;
  // initialize static stuff
  m_NewWisdomAvailable = false;
  
  SetGlobalOptimizationLevel(ESTIMATE); // the default value
  std::string opt;
  if( itksys::SystemTools::GetEnv("ITK_FFTW_OPTIMIZATION_LEVEL", opt) )
    {
    if( opt == "PATIENT" )
      {
      SetGlobalOptimizationLevel(PATIENT);
      }
    else if( opt == "EXHAUSTIVE" )
      {
      SetGlobalOptimizationLevel(EXHAUSTIVE);
      }
    }
  // SetWisdomFileDefaultBaseName(std::string(itksys::SystemTools::GetEnv("HOME")) + "/.itkwisdom");

#if defined(USE_FFTWF)
  fftwf_init_threads();
  fftwf_init_threads();
#endif
#if defined(USE_FFTWD)
  fftw_init_threads();
  fftw_init_threads();
#endif

  std::string auto_import_env;
  if( !itksys::SystemTools::GetEnv("ITK_FFTW_WISDOM_AUTO_IMPORT", auto_import_env) ||
      ( auto_import_env != "no" &&
      auto_import_env != "NO" && 
      auto_import_env != "off" && 
      auto_import_env != "OFF" &&
      auto_import_env != "0" ) )
    {
#if defined(USE_FFTWF)
      fftwf_import_system_wisdom();
      ImportWisdomFileFloat();
#endif
#if defined(USE_FFTWD)
      fftw_import_system_wisdom();
      ImportWisdomFileDouble();
#endif
    }
}

FFTWLock
::~FFTWLock()
{
  // std::cout << "======== cleanup fftw stuff =========" << std::endl;
  std::string auto_export_env;
  if( m_NewWisdomAvailable && ( !itksys::SystemTools::GetEnv("ITK_FFTW_WISDOM_AUTO_IMPORT", auto_export_env) ||
      ( auto_export_env != "no" &&
      auto_export_env != "NO" && 
      auto_export_env != "off" && 
      auto_export_env != "OFF" &&
      auto_export_env != "0" ) ) )
    {
    // import the wisdom files again to be sure to not erase the wisdom saved in another process
    ImportWisdomFileFloat();
    ExportWisdomFileFloat();
    
    ImportWisdomFileDouble();
    ExportWisdomFileDouble();
    }
#if defined(USE_FFTWF)
  fftwf_cleanup_threads();
  fftwf_cleanup();
#endif
#if defined(USE_FFTWD)
  fftw_cleanup_threads();
  fftw_cleanup();
#endif
}

// void
// FFTWLock
// ::SetWisdomFileDefaultBaseName(std::string s)
// {
//   m_WisdomFileDefaultBaseName = s;
// }

std::string
FFTWLock
::GetWisdomFileDefaultBaseName()
{
  // return m_WisdomFileDefaultBaseName;
  std::string name;
  if( !itksys::SystemTools::GetEnv("ITK_FFTW_WISDOM_FILE_BASE_NAME", name) )
    {
    name = std::string(itksys::SystemTools::GetEnv("HOME")) + "/.itkwisdom";
    }
  return name;
}

bool
FFTWLock
::ImportWisdomFileFloat()
{
  return ImportWisdomFileFloat( GetWisdomFileDefaultBaseName() + "f" );
}

bool
FFTWLock
::ImportWisdomFileDouble()
{
  return ImportWisdomFileDouble( GetWisdomFileDefaultBaseName() );
}

bool
FFTWLock
::ExportWisdomFileFloat()
{
  return ExportWisdomFileFloat( GetWisdomFileDefaultBaseName() + "f" );
}

bool
FFTWLock
::ExportWisdomFileDouble()
{
  return ExportWisdomFileDouble( GetWisdomFileDefaultBaseName() );
}

bool
FFTWLock
::ImportWisdomFileFloat( std::string path )
{
  bool ret = false;
#if defined(USE_FFTWF)
  FILE * f = fopen( path.c_str(), "r" );
  if( f )
    {
    flock( fileno(f), LOCK_SH );
    ret = fftwf_import_wisdom_from_file( f );
    flock( fileno(f), LOCK_UN );
    fclose( f );
    }
#endif
  return ret;
}

bool
FFTWLock
::ImportWisdomFileDouble( std::string path )
{
  bool ret = false;
#if defined(USE_FFTWD)
  FILE * f = fopen( path.c_str(), "r" );
  if( f )
    {
    flock( fileno(f), LOCK_SH );
    ret = fftw_import_wisdom_from_file( f );
    flock( fileno(f), LOCK_UN );
    fclose( f );
    }
#endif
  return ret;
}
  
bool
FFTWLock
::ExportWisdomFileFloat( std::string path )
{
  bool ret = false;
#if defined(USE_FFTWF)
  FILE * f = fopen( path.c_str(), "w" );
  if( f )
    {
    flock( fileno(f), LOCK_EX );
    fftwf_export_wisdom_to_file( f );
    flock( fileno(f), LOCK_UN );
    ret = fclose( f );
    }
#endif
  return ret;
}

bool
FFTWLock
::ExportWisdomFileDouble( std::string path )
{
  bool ret = false;
#if defined(USE_FFTWD)
  FILE * f = fopen( path.c_str(), "w" );
  if( f )
    {
    flock( fileno(f), LOCK_EX );
    fftw_export_wisdom_to_file( f );
    flock( fileno(f), LOCK_UN );
    ret = fclose( f );
    }
#endif
  return ret;
}

  
void
FFTWLock
::Lock()
{
  FFTWLock::m_Lock.Lock();
}
  
void
FFTWLock
::Unlock()
{
  FFTWLock::m_Lock.Unlock();
}
  
void
FFTWLock
::NewWisdomAvailable()
{
  FFTWLock::m_NewWisdomAvailable = true;
}
  
int
FFTWLock
::GetGlobalOptimizationLevel()
{
  return m_GlobalOptimizationLevel;
}

void
FFTWLock
::SetGlobalOptimizationLevel(int opt)
{
  m_GlobalOptimizationLevel = opt;
}


}//end namespace itk
