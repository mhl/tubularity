/* A trivial C++ program for testing JNI calls */

#include <iostream>

#include "FijiITKInterface_OOFTubularityMeasure.h"
#include "itkImageFileWriter.h"
#include "itkMultiScaleOrientedFluxBasedMeasureFFTImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkOrientedFluxTraceMeasure.h"
#include "itkOrientedFluxCrossSectionTraceMeasure.h"


#define SwitchCase(CaseValue, DerivedFilterType, BaseFilterObjectPtr, Call ) \
case CaseValue: \
{ \
typedef DerivedFilterType FilterObjectType; \
typename FilterObjectType::Pointer FilterObjectPtr = static_cast<FilterObjectType*>(BaseFilterObjectPtr.GetPointer()); \
Call; \
break; \
} \

#define MultiScaleEnhancementFilterSwitch3D(HessianFilterTypeValue, BaseFilterObjectPtr, Call)  \
switch( HessianFilterTypeValue ) \
{ \
SwitchCase(OrientedFluxTrace, OrientedFluxTraceMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
SwitchCase(OrientedFluxCrossSectionCurvature, OrientedFluxMainCurvatureMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
}\

#define MultiScaleEnhancementFilterSwitchND(HessianFilterTypeValue, BaseFilterObjectPtr, Call)  \
switch( HessianFilterTypeValue ) \
{ \
SwitchCase(OrientedFluxTrace, OrientedFluxTraceMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
SwitchCase(OrientedFluxCrossSectionCurvature, OrientedFluxMainCurvatureMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
}\
\

#define GRAY8 0
#define GRAY16 2
#define GRAY32 4

using std::cout;
using std::endl;

const unsigned int maxDimension = 3;




// Main code goes here! 
template<class TInputPixel, unsigned int VDimension> 
typename itk::Image<float,VDimension+1>::Pointer
Execute(typename itk::Image<TInputPixel,VDimension>::Pointer Input_Image, double sigmaMin, double sigmaMax, unsigned int numberOfScales)
{	
	// Define the dimension of the images
	const unsigned int Dimension = VDimension;
	
	// Default hessian filter type
	typedef enum
	{
		OrientedFluxTrace = 2,
		OrientedFluxCrossSectionCurvature = 3
	}HessianFilterTypeEnum;
	HessianFilterTypeEnum HessianFilterTypeValue;
	
	
	// Typedefs
	typedef TInputPixel											InputPixelType;
	typedef itk::Image<InputPixelType,Dimension>								InputImageType;
	typedef typename InputImageType::SpacingType								SpacingType;
	
	
	typedef float												OutputPixelType;
	typedef itk::Image<OutputPixelType,Dimension>								OutputImageType;
	typedef itk::Image<OutputPixelType,Dimension+1>								OutputScaleSpaceImageType;
	
	
	typedef float												HessianPixelScalarType;
	typedef itk::SymmetricSecondRankTensor< HessianPixelScalarType, Dimension > 				HessianPixelType;
	typedef itk::Image< HessianPixelType, Dimension >							HessianImageType;
	
	typedef float												ScalesPixelType;
	typedef itk::Image<ScalesPixelType, Dimension>								ScalesImageType;

	
	typedef itk::ShiftScaleImageFilter<OutputImageType, OutputImageType>					ShiftScaleFilterType;
	typedef itk::MinimumMaximumImageCalculator<OutputImageType>						MinMaxCalculatorType;
	typedef itk::ShiftScaleImageFilter<OutputScaleSpaceImageType, OutputScaleSpaceImageType>		ShiftScaleFilterForScaleSpaceImageType;
	typedef itk::MinimumMaximumImageCalculator<OutputScaleSpaceImageType>					MinMaxCalculatorForScaleSpaceImageType;
	
	
	// Declare the type of enhancement filter
	typedef itk::ProcessObject ObjectnessBaseFilterType;
	typedef itk::OrientedFluxTraceMeasureFilter<HessianImageType,OutputImageType> 		HessianToOrientedFluxTraceObjectnessFilterType;	
	typedef itk::OrientedFluxCrossSectionTraceMeasureFilter<HessianImageType,OutputImageType> 	HessianToOrientedFluxMainCurvatureObjectnessFilterType;	
	
	// Declare the type of multiscale enhancement filter
	typedef itk::ProcessObject 										MultiScaleEnhancementBaseFilterType;

	typedef itk::MultiScaleOrientedFluxBasedMeasureFFTImageFilter< InputImageType, 
								HessianImageType, 
								ScalesImageType,
								HessianToOrientedFluxTraceObjectnessFilterType, 
								OutputImageType > 				OrientedFluxTraceMultiScaleEnhancementFilterType;
	typedef itk::MultiScaleOrientedFluxBasedMeasureFFTImageFilter< InputImageType, 
								HessianImageType, 
								ScalesImageType,
								HessianToOrientedFluxMainCurvatureObjectnessFilterType, 
								OutputImageType > 				OrientedFluxMainCurvatureMultiScaleEnhancementFilterType;	


	SpacingType spacing = Input_Image->GetSpacing();
	double maxSpacing = spacing[0];
	double minSpacing = spacing[0];
	for(unsigned int i = 1; i < Dimension; i++)
	{
		maxSpacing = vnl_math_max(maxSpacing, spacing[i]);
		minSpacing = vnl_math_min(minSpacing, spacing[i]);
	}
	// Parse the input arguments.

	double fixedSigmaForHessianComputation = minSpacing;//TODO : use the minimal ImageSpacing
	
	bool brightObject = true;
	//true	
	//bool normalizeTubularityImageBtw0And1 = true;
	// for first try, fix to 3 (cross section eigen trace)	
	HessianFilterTypeValue =  OrientedFluxCrossSectionCurvature;//OrientedFluxCrossSectionCurvature;
	bool useAFixedSigmaForComputingHessianImage = true;
	//bool scaleObjectnessMeasure = false;
	//unsigned int objectDimension = 1; // default values just to avoid compiler 


	//bool generateScaleImage = false;
	bool generateScaleSpaceTubularityScoreImage = true;

       	// typenames needed
	ObjectnessBaseFilterType::Pointer objectnessFilter;
	MultiScaleEnhancementBaseFilterType::Pointer multiScaleEnhancementFilter;
	typename HessianToOrientedFluxMainCurvatureObjectnessFilterType::Pointer orientedFluxMainCurvatureObjectnessFilter = HessianToOrientedFluxMainCurvatureObjectnessFilterType::New();
	orientedFluxMainCurvatureObjectnessFilter->SetBrightObject( brightObject );
	objectnessFilter = orientedFluxMainCurvatureObjectnessFilter;
	
	typename OrientedFluxMainCurvatureMultiScaleEnhancementFilterType::Pointer orientedFluxMainCurvatureMultiScaleEnhancementFilter = OrientedFluxMainCurvatureMultiScaleEnhancementFilterType::New();
	orientedFluxMainCurvatureMultiScaleEnhancementFilter->SetOrientedFluxToMeasureFilter( orientedFluxMainCurvatureObjectnessFilter );
	multiScaleEnhancementFilter = orientedFluxMainCurvatureMultiScaleEnhancementFilter;		
	// main function
	MultiScaleEnhancementFilterSwitchND(
		HessianFilterTypeValue, multiScaleEnhancementFilter,FilterObjectPtr->SetInput(Input_Image);
		FilterObjectPtr->SetSigmaMinimum( sigmaMin ); 
		FilterObjectPtr->SetSigmaMaximum( sigmaMax );  
		FilterObjectPtr->SetNumberOfSigmaSteps( numberOfScales );
		FilterObjectPtr->SetGenerateNPlus1DHessianMeasureOutput(generateScaleSpaceTubularityScoreImage);
  
		if( useAFixedSigmaForComputingHessianImage )
		{
			FilterObjectPtr->SetFixedSigmaForHessianImage( fixedSigmaForHessianComputation );
		}
		FilterObjectPtr->SetGenerateHessianOutput( false );//false
		
		try
		{
			FilterObjectPtr->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cerr << e << std::endl;
		}
										
		// Writing the output image.
		typename OutputScaleSpaceImageType::Pointer tubularityScoreImage;
		typename MinMaxCalculatorForScaleSpaceImageType::Pointer minMaxCalc = MinMaxCalculatorForScaleSpaceImageType::New();
		minMaxCalc->SetImage( FilterObjectPtr->GetNPlus1DImageOutput() );
		minMaxCalc->Compute();
		typename ShiftScaleFilterForScaleSpaceImageType::Pointer shiftScaleFilter = ShiftScaleFilterForScaleSpaceImageType::New();
		shiftScaleFilter->SetInput( FilterObjectPtr->GetNPlus1DImageOutput() );
		shiftScaleFilter->SetShift( -minMaxCalc->GetMinimum() );
		shiftScaleFilter->SetScale( 1 / (minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum()) );
		shiftScaleFilter->Update();
		tubularityScoreImage = shiftScaleFilter->GetOutput();
		return tubularityScoreImage;															
	)		// end MultiScaleEnhancementFilterSwitchND
	

	return EXIT_SUCCESS;
}


JNIEXPORT jint JNICALL Java_FijiITKInterface_OOFTubularityMeasure_OrientedFlux(JNIEnv *env, jobject ignored, jbyteArray jba, jfloatArray jbOut, jint type, jint width, jint height, jint NSlice, jdouble widthpix, jdouble heightpix, jdouble depthpix, jdouble sigmaMin, jdouble sigmaMax, jint numberOfScales, jstring outputFileName)
{
    jboolean isCopy;
    jbyte * jbs   = env->GetByteArrayElements(jba,&isCopy);
    jfloat * jbOutS = env->GetFloatArrayElements(jbOut,&isCopy);

    if( ! jbs )
        return -1;

	double spacing[3], origin[3]; 

	unsigned char * InputImageData = (unsigned char *) jbs;
	//Allocates an itk image with the 2D buffered data
   	typedef itk::Image<unsigned char, 3>  ImageType;
	typedef itk::Image<float, 4> 	      OutputImageType;

	ImageType::Pointer itkImageP = ImageType::New();

	ImageType::SizeType size;
	size[0] = width;size[1] = height;size[2] = NSlice;

	ImageType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	
	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	itkImageP->SetRegions( region);
	itkImageP->Allocate();

	spacing[0] = widthpix;spacing[1] = heightpix;spacing[2] = depthpix;
	itkImageP->SetSpacing( spacing );

	origin[0] = 0;origin[1] = 0;origin[2] = 0;
	itkImageP->SetOrigin( origin );

	//Copies the data from the image buffer to the itkImage
	typedef itk::ImageRegionIterator< ImageType> IteratorType;
	typedef itk::ImageRegionIterator< OutputImageType> OutputIteratorType;
	IteratorType it( itkImageP, region);

	it.GoToBegin();
	unsigned char * dataPointer = InputImageData;
	while( ! it.IsAtEnd() )
	{
		it.Set( *dataPointer);
		++it;
		++dataPointer;
	 }

	OutputImageType::Pointer outputImage = Execute<unsigned char, 3>(itkImageP, sigmaMin, sigmaMax, numberOfScales);

	OutputImageType::RegionType Outputregion;
	OutputImageType::SizeType Outputsize;
	Outputsize[0] = width;Outputsize[1] = height; Outputsize[2] = NSlice; Outputsize[3] = numberOfScales;
	OutputImageType::IndexType Outputstart;
	Outputstart[0] = 0;Outputstart[1] = 0;Outputstart[2] = 0; Outputstart[3] = 0;
	Outputregion.SetSize( Outputsize );
	Outputregion.SetIndex( Outputstart );
	
	float* outputImageData = (float*) jbOutS;
	std::cout << Outputregion << std::endl;	
	OutputIteratorType outit( outputImage, Outputregion);
	std::cout << outputImage << std::endl;
	outit.GoToBegin();
	int length = width * height * NSlice;
	for( int i = 0; i < length; ++i ) {
	  	outputImageData[i] = outit.Get();
	  	++outit;
	}
	

	 const char *s = env->GetStringUTFChars(outputFileName,NULL);
	 if( ! s )
		std::cerr << "Converting and allocating the filename string failed" << std::endl;

	typedef itk::ImageFileWriter<OutputImageType>  ImageWriter;
	ImageWriter::Pointer writer = ImageWriter::New();
	writer->SetInput(outputImage);
	writer->SetFileName(s);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}

	

  env->ReleaseByteArrayElements(jba,jbs,0);
  env->ReleaseFloatArrayElements(jbOut, jbOutS,0);
	env->ReleaseStringUTFChars( outputFileName, s );
    return 0;
}

