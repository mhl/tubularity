/* A trivial C++ program for testing JNI calls */

#include <iostream>

#include "FijiITKInterface_TubularityMeasure.h"

//#include "itkHessianMainPrincipleCurvatureObjectnessImageFilter.h"
//#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter2.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"

#include "itkHessianToOrientedFluxTraceMeasureImageFilter.h"
#include "itkHessianToOrientedFluxMainCurvatureMeasureImageFilter.h"

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
SwitchCase(Lindeberg, LindebergMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
SwitchCase(OrientedFluxTrace, OrientedFluxTraceMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
SwitchCase(OrientedFluxCrossSectionCurvature, OrientedFluxMainCurvatureMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
}\

#define MultiScaleEnhancementFilterSwitchND(HessianFilterTypeValue, BaseFilterObjectPtr, Call)  \
switch( HessianFilterTypeValue ) \
{ \
SwitchCase(Lindeberg, LindebergMultiScaleEnhancementFilterType, BaseFilterObjectPtr, Call ) \
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
typename itk::Image<float,VDimension>::Pointer
Execute(typename itk::Image<TInputPixel,VDimension>::Pointer Input_Image, double sigmaMin, double sigmaMax, unsigned int numberOfScales)
{	
	// Define the dimension of the images
	const unsigned int Dimension = VDimension;
	
	// Default hessian filter type
	typedef enum
	{
		Lindeberg = 0,
		Frangi = 1,
		OrientedFluxTrace = 2,
		OrientedFluxCrossSectionCurvature = 3
	}HessianFilterTypeEnum;
	HessianFilterTypeEnum HessianFilterTypeValue = Lindeberg;
	
	
	// Typedefs
	typedef TInputPixel											InputPixelType;
	typedef itk::Image<InputPixelType,Dimension>								InputImageType;
	typedef typename InputImageType::SpacingType								SpacingType;
	
	
	typedef float												OutputPixelType;
	typedef itk::Image<OutputPixelType,Dimension>				OutputImageType;
	typedef itk::Image<OutputPixelType,Dimension+1>				OutputScaleSpaceImageType;
	
	
	typedef float																						HessianPixelScalarType;
	typedef itk::SymmetricSecondRankTensor< HessianPixelScalarType, Dimension > HessianPixelType;
	typedef itk::Image< HessianPixelType, Dimension >		HessianImageType;
	
	typedef float																				ScalesPixelType;
	typedef itk::Image<ScalesPixelType, Dimension>			ScalesImageType;

	
	typedef itk::ShiftScaleImageFilter<OutputImageType, OutputImageType>	ShiftScaleFilterType;
	typedef itk::MinimumMaximumImageCalculator<OutputImageType>						MinMaxCalculatorType;
	typedef itk::ShiftScaleImageFilter<OutputScaleSpaceImageType, OutputScaleSpaceImageType>	ShiftScaleFilterForScaleSpaceImageType;
	typedef itk::MinimumMaximumImageCalculator<OutputScaleSpaceImageType>											MinMaxCalculatorForScaleSpaceImageType;
	
	
	// Declare the type of enhancement filter
	typedef itk::ProcessObject ObjectnessBaseFilterType;
	typedef itk::HessianMainPrincipleCurvatureObjectnessImageFilter <HessianImageType,OutputImageType> LindebergObjectnessFilterType;
	//typedef itk::HessianToObjectnessMeasureImageFilter <HessianImageType,OutputImageType> FrangiObjectnessFilterType;
	typedef itk::HessianToOrientedFluxTraceMeasureFilter <HessianImageType,OutputImageType> HessianToOrientedFluxTraceObjectnessFilterType;	
	typedef itk::HessianToOrientedFluxMainCurvatureMeasureFilter <HessianImageType,OutputImageType> HessianToOrientedFluxMainCurvatureObjectnessFilterType;	
	
	// Declare the type of multiscale enhancement filter
	typedef itk::ProcessObject MultiScaleEnhancementBaseFilterType;
	typedef itk::MultiScaleHessianBasedMeasureImageFilter2< InputImageType, 
															HessianImageType, 
															ScalesImageType,
															LindebergObjectnessFilterType, 
															OutputImageType > LindebergMultiScaleEnhancementFilterType;
	typedef itk::MultiScaleHessianBasedMeasureImageFilter2< InputImageType, 
															HessianImageType, 
															ScalesImageType,
															HessianToOrientedFluxTraceObjectnessFilterType, 
															OutputImageType > OrientedFluxTraceMultiScaleEnhancementFilterType;
	typedef itk::MultiScaleHessianBasedMeasureImageFilter2< InputImageType, 
															HessianImageType, 
															ScalesImageType,
															HessianToOrientedFluxMainCurvatureObjectnessFilterType, 
															OutputImageType > OrientedFluxMainCurvatureMultiScaleEnhancementFilterType;	


	SpacingType spacing = Input_Image->GetSpacing();
	double maxSpacing = spacing[0];
	double minSpacing = spacing[0];
	for(unsigned int i = 1; i < Dimension; i++)
	{
		maxSpacing = vnl_math_max(maxSpacing, spacing[i]);
		minSpacing = vnl_math_min(minSpacing, spacing[i]);
	}
	// Parse the input arguments.
	unsigned int argumentOffset = 1;

	double fixedSigmaForHessianComputation = minSpacing;//TODO : use the minimal ImageSpacing
	
	bool brightObject = true;
	//true	
	bool normalizeTubularityImageBtw0And1 = true;
	// for first try, fix to 3 (cross section eigen trace)	
	HessianFilterTypeValue =  OrientedFluxCrossSectionCurvature;//OrientedFluxCrossSectionCurvature;
	bool useAFixedSigmaForComputingHessianImage = true;
	bool scaleObjectnessMeasure = false;
	unsigned int objectDimension = 1; // default values just to avoid compiler 


	bool generateScaleImage = false;
	bool generateHessianMatrixImage = false;
	bool generateScaleSpaceTubularityScoreImage = false;

	
       	// typenames needed
	ObjectnessBaseFilterType::Pointer objectnessFilter;
	MultiScaleEnhancementBaseFilterType::Pointer multiScaleEnhancementFilter;
	typename HessianToOrientedFluxMainCurvatureObjectnessFilterType::Pointer orientedFluxMainCurvatureObjectnessFilter = HessianToOrientedFluxMainCurvatureObjectnessFilterType::New();
	orientedFluxMainCurvatureObjectnessFilter->SetBrightObject( brightObject );
	objectnessFilter = orientedFluxMainCurvatureObjectnessFilter;
	
	typename OrientedFluxMainCurvatureMultiScaleEnhancementFilterType::Pointer orientedFluxMainCurvatureMultiScaleEnhancementFilter = OrientedFluxMainCurvatureMultiScaleEnhancementFilterType::New();
	orientedFluxMainCurvatureMultiScaleEnhancementFilter->SetHessianToMeasureFilter( orientedFluxMainCurvatureObjectnessFilter );
	multiScaleEnhancementFilter = orientedFluxMainCurvatureMultiScaleEnhancementFilter;		

	MultiScaleEnhancementFilterSwitchND(
		HessianFilterTypeValue, multiScaleEnhancementFilter,FilterObjectPtr->SetInput(Input_Image);
		FilterObjectPtr->SetSigmaMinimum( sigmaMin ); 
		FilterObjectPtr->SetSigmaMaximum( sigmaMax );  
		FilterObjectPtr->SetNumberOfSigmaSteps( numberOfScales );
		FilterObjectPtr->SetUseAFixedSigmaForComputingHessianImage( useAFixedSigmaForComputingHessianImage );

  
	if( useAFixedSigmaForComputingHessianImage )
	{
		FilterObjectPtr->SetFixedSigmaForHessianImage( fixedSigmaForHessianComputation );
	}
		FilterObjectPtr->SetGenerateHessianOutput( generateHessianMatrixImage );//false
		
		std::cout << Input_Image.GetPointer() << std::endl;

		std::cout << Input_Image->GetSpacing() << std::endl;
		
		try
		{
			
			FilterObjectPtr->Update();
			std::cout << "FilterObjectPtr->Update();" << std::endl;
		}
		catch (itk::ExceptionObject &e)
		{
			std::cerr << e << std::endl;
		}
										
		// Writing the output image.
		typename OutputImageType::Pointer tubularityScoreImage;
		if( normalizeTubularityImageBtw0And1 ){
			typename MinMaxCalculatorType::Pointer minMaxCalc = MinMaxCalculatorType::New();
			minMaxCalc->SetImage( FilterObjectPtr->GetOutput() );
			minMaxCalc->Compute();
																				
			typename ShiftScaleFilterType::Pointer shiftScaleFilter = ShiftScaleFilterType::New();
			shiftScaleFilter->SetInput( FilterObjectPtr->GetOutput() );
			shiftScaleFilter->SetShift( -minMaxCalc->GetMinimum() );
			shiftScaleFilter->SetScale( 1 / (minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum()) );
			shiftScaleFilter->Update();
			tubularityScoreImage = shiftScaleFilter->GetOutput();
		}
		else{
			tubularityScoreImage = FilterObjectPtr->GetOutput();
		}

		return tubularityScoreImage;
																			
	)		// end MultiScaleEnhancementFilterSwitchND
	std::cout << "Exiting with success." << std::endl;
	
	
}


JNIEXPORT jint JNICALL Java_FijiITKInterface_TubularityMeasure_nonsense(JNIEnv * env, jobject ignored, jstring js)
{
    const char *s = env->GetStringUTFChars(js,NULL);
    if( ! s )
        return -1;
    cout << "We were passed the nonsense string: " << s << endl;
    env->ReleaseStringUTFChars( js, s );
    return 42;
}

JNIEXPORT jint JNICALL Java_FijiITKInterface_TubularityMeasure_OrientedFlux(JNIEnv *env, jobject ignored, jbyteArray jba, jfloatArray jbOut, jint type, jint width, jint height, jint NSlice, jdouble widthpix, jdouble heightpix, jdouble depthpix, jdouble sigmaMin, jdouble sigmaMax, jint numberOfScales)
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
	typedef itk::Image<float, 3> 	      OutputImageType;

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

	// Now it comes the itk filtering structure
	/*
	typedef itk::SimpleOperationFilter<ImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( itkImageP );
	filter->Update();*/

	OutputImageType::Pointer outputImage = Execute<unsigned char, 3>(itkImageP, sigmaMin, sigmaMax, numberOfScales);
	float* outputImageData = (float*) jbOutS;
	OutputIteratorType outit( outputImage, region);

	outit.GoToBegin();
	int length = width * height * NSlice;
	for( int i = 0; i < length; ++i ) {
	  	outputImageData[i] = outit.Get();
	  	++outit;
	}

    env->ReleaseByteArrayElements(jba,jbs,0);
    env->ReleaseFloatArrayElements(jbOut, jbOutS,0);

    return 0;
}


