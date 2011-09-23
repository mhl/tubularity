/* A trivial C++ program for testing JNI calls */

#include <iostream>

#include "FijiITKInterface_MultiScaleTubularityMeasure.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter2.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkTubularMetricToPathFilter3.h"
#include <itkImageFileWriter.h>
#include "itkLogImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"

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
using std::cerr;
using std::endl;

const unsigned int maxDimension = 3;
std::vector< float > Outputpath;


// Main code goes here! 
template<class TInputPixel, unsigned int VDimension> 
typename itk::Image<float,VDimension+1>::Pointer
Execute(typename itk::Image<TInputPixel,VDimension>::Pointer Input_Image, double sigmaMin, double sigmaMax, unsigned int numberOfScales, float* pt1, float* pt2)
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
	typedef TInputPixel										InputPixelType;
	typedef itk::Image<InputPixelType,Dimension>							InputImageType;
	typedef typename InputImageType::SpacingType							SpacingType;
	
	
	typedef float											OutputPixelType;
	typedef itk::Image<OutputPixelType,Dimension>							OutputImageType;
	typedef itk::Image<OutputPixelType,Dimension+1>							OutputScaleSpaceImageType;
	
	
	typedef float											HessianPixelScalarType;
	typedef itk::SymmetricSecondRankTensor< HessianPixelScalarType, Dimension > 			HessianPixelType;
	typedef itk::Image< HessianPixelType, Dimension >						HessianImageType;
	
	typedef float											ScalesPixelType;
	typedef itk::Image<ScalesPixelType, Dimension>							ScalesImageType;


	typedef itk::ShiftScaleImageFilter<OutputImageType, OutputImageType>				ShiftScaleFilterType;
	typedef itk::MinimumMaximumImageCalculator<OutputImageType>					MinMaxCalculatorType;
	typedef itk::ShiftScaleImageFilter<OutputScaleSpaceImageType, OutputScaleSpaceImageType>	ShiftScaleFilterForScaleSpaceImageType;
	typedef itk::MinimumMaximumImageCalculator<OutputScaleSpaceImageType>				MinMaxCalculatorForScaleSpaceImageType;
	
	
	// Declare the type of enhancement filter
	typedef itk::ProcessObject ObjectnessBaseFilterType;
	typedef itk::HessianToOrientedFluxTraceMeasureFilter <HessianImageType,OutputImageType> HessianToOrientedFluxTraceObjectnessFilterType;	
	typedef itk::HessianToOrientedFluxMainCurvatureMeasureFilter <HessianImageType,OutputImageType> HessianToOrientedFluxMainCurvatureObjectnessFilterType;	
	
	// Declare the type of multiscale enhancement filter
	typedef itk::ProcessObject MultiScaleEnhancementBaseFilterType;
	
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


	typedef typename OutputScaleSpaceImageType::IndexType IndexType;
	typedef itk::TubularMetricToPathFilter3< OutputScaleSpaceImageType >  PathFilterType;
	typedef typename PathFilterType::VertexType										VertexType;
	

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
	// for first try, fix to 3 (cross section eigen trace)	
	HessianFilterTypeValue =  OrientedFluxCrossSectionCurvature;//OrientedFluxCrossSectionCurvature;
	bool useAFixedSigmaForComputingHessianImage = true;
 	bool generateHessianMatrixImage = false;
	bool generateScaleSpaceTubularityScoreImage = true;

	
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
		FilterObjectPtr->SetGenerateNPlus1DHessianMeasureOutput(generateScaleSpaceTubularityScoreImage);

  
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
		typename OutputScaleSpaceImageType::Pointer tubularityScoreImage;
		
			
		  typename MinMaxCalculatorForScaleSpaceImageType::Pointer minMaxCalc = MinMaxCalculatorForScaleSpaceImageType::New();
			minMaxCalc->SetImage( (FilterObjectPtr->GetNPlus1DImageOutput()) );
			minMaxCalc->Compute();
																				
			typename ShiftScaleFilterForScaleSpaceImageType::Pointer shiftScaleFilter = ShiftScaleFilterForScaleSpaceImageType::New();
			shiftScaleFilter->SetInput( FilterObjectPtr->GetNPlus1DImageOutput() );
			shiftScaleFilter->SetShift( -minMaxCalc->GetMinimum() );
			shiftScaleFilter->SetScale( 1 / (minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum()) );
			shiftScaleFilter->Update();
			tubularityScoreImage = shiftScaleFilter->GetOutput();
		
		


	////code test///
	typename PathFilterType::Pointer pathFilter = PathFilterType::New();
	typename PathFilterType::Pointer path = PathFilterType::New();

	pathFilter->SetInput( tubularityScoreImage );

	IndexType startPoint;
	startPoint[0] = pt1[0];
	startPoint[1] = pt1[1];
	startPoint[2] = pt1[2];
	startPoint[3] = 0;

	std::cout << "startpoint " << startPoint[0] << ", " << startPoint[1] << ", " << startPoint[2] << std::endl;

	IndexType endPoint;
	endPoint[0] = pt2[0];
	endPoint[1] = pt2[1];
	endPoint[2] = pt2[2];
	endPoint[3] = 0;

	std::cout << "endpoint " << endPoint[0] << ", " << endPoint[1] << ", " << endPoint[2] << std::endl;

	pathFilter->SetStartPoint(startPoint);
	pathFilter->AddPathEndPoint(endPoint);
	
	std::cout<< "path update" << std::endl;

	pathFilter->Update();

	typedef typename OutputScaleSpaceImageType::SpacingType ScaleSpaceSpacingType;
	typedef typename OutputScaleSpaceImageType::PointType  ScaleSpaceOriginType;

	ScaleSpaceSpacingType ss_spacing = tubularityScoreImage->GetSpacing();
	ScaleSpaceOriginType  ss_origin = tubularityScoreImage->GetOrigin();

	for(unsigned int k = 0; k < pathFilter->GetPath(0)->GetVertexList()->Size(); k++)
	{      
	  VertexType vertex = pathFilter->GetPath(0)->GetVertexList()->GetElement(k);	     
	  for (unsigned int i = 0; i < Dimension+1; i++)
	    {	     
	      Outputpath.push_back(vertex[i]*ss_spacing[i]+ss_origin[i]);
	    }
	}

	pathFilter->WritePathsToFile("testPathnew");
	////////////////

	return tubularityScoreImage;
																			
	)		// end MultiScaleEnhancementFilterSwitchND
	std::cout << "Exiting with success." << std::endl;
	return EXIT_SUCCESS;
	
}

JNIEXPORT jint JNICALL Java_FijiITKInterface_MultiScaleTubularityMeasure_GetPath(JNIEnv * env, jobject ignored, jfloatArray jba)
{
	jboolean isCopy;
	jfloat * jbOutS = env->GetFloatArrayElements(jba,&isCopy);
	float* outputpath = (float*) jbOutS;
	

	//vertex[SetDimension-1] * spacing[SetDimension-1] + origin[SetDimension-1]
	for(unsigned int i=0;i<Outputpath.size();i++)
		outputpath[i] = Outputpath[i];

	 env->ReleaseFloatArrayElements(jba, jbOutS,0);
	
	return Outputpath.size();
}

JNIEXPORT jint JNICALL Java_FijiITKInterface_MultiScaleTubularityMeasure_FindPath(JNIEnv * env, jobject ignored, jfloatArray jba, jintArray point1, jintArray point2, jfloatArray Path, jint width, jint height, jint NSlice, jint NScale, jdouble widthpix, jdouble heightpix, jdouble depthpix){
	
	double spacing[3], origin[3];
	// compute the min path	
	typedef float	      TubularityScorePixelType;
	typedef itk::Image<TubularityScorePixelType,4>							TubularityScoreImageType;

	typedef itk::TubularMetricToPathFilter3< TubularityScoreImageType >  PathFilterType;
	PathFilterType::Pointer pathFilter = PathFilterType::New();


	//Copies the data from the image buffer to the itkImage
	jboolean isCopy;
	jfloat * jbs   = env->GetFloatArrayElements(jba,&isCopy);
	TubularityScorePixelType * InputImageData = (TubularityScorePixelType *) jbs;

	jfloat * output   = env->GetFloatArrayElements(Path,&isCopy);
	//float * OutputPath = (float *) output;

	TubularityScoreImageType::Pointer itkImageP = TubularityScoreImageType::New();
	TubularityScoreImageType::SizeType size;
	size[0] = width;size[1] = height;size[2] = NSlice; size[3] = NScale;

	TubularityScoreImageType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0; start[3] = 0;
	
	TubularityScoreImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	itkImageP->SetRegions( region);
	itkImageP->Allocate();

	spacing[0] = widthpix;spacing[1] = heightpix;spacing[2] = depthpix; //spacing[3] = 
	itkImageP->SetSpacing( spacing );

	origin[0] = 0;origin[1] = 0;origin[2] = 0; origin[3] = 1;
	itkImageP->SetOrigin( origin );
	

	typedef itk::ImageRegionIterator< TubularityScoreImageType> IteratorType;
	IteratorType it( itkImageP, region);

	it.GoToBegin();
	float * dataPointer = InputImageData;
	while( ! it.IsAtEnd() )
	{
		it.Set( *dataPointer);
		++it;
		++dataPointer;
	 }
	
	typedef itk::ImageFileWriter<TubularityScoreImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(itkImageP);	
	writer->SetFileName("ZopAvantLog.nrrd");	
	writer->Update();

	// Take the negative log-likelihood of the probability image.

	typedef itk::LogImageFilter
		<TubularityScoreImageType, TubularityScoreImageType> 		LogImageFilterType;
	typedef itk::MultiplyByConstantImageFilter<TubularityScoreImageType,
		double, TubularityScoreImageType>				MultiplyImageFilterType;

	LogImageFilterType::Pointer logFilter = LogImageFilterType::New();
	logFilter->SetInput( itkImageP );
	MultiplyImageFilterType::Pointer multipFilter = MultiplyImageFilterType::New();
	multipFilter->SetInput( logFilter->GetOutput() );
	multipFilter->SetConstant( -1.0 );
	multipFilter->Update();
	itkImageP = multipFilter->GetOutput();
	itkImageP->DisconnectPipeline();
	//logFilter = 0;
	//multipFilter = 0;
	writer->SetInput(itkImageP);
	writer->SetFileName("ZopApresLog.nrrd");	
	writer->Update();


	it.GoToBegin();
	while( ! it.IsAtEnd() )
	{
		it.Set( 1 / it.Get() );
		++it;
	 }

	writer->SetFileName("ZopApresInv.nrrd");	
	writer->Update();


	std::cout<< "path filter set input" << std::endl;

	pathFilter->SetInput( itkImageP );
	
	typedef TubularityScoreImageType::IndexType IndexType;
	typedef std::vector<IndexType> ListIndices;

	jint* pt1   = env->GetIntArrayElements(point1,&isCopy);
	jint* pt2   = env->GetIntArrayElements(point2,&isCopy);

	IndexType startPoint;
	startPoint[0] = pt1[0];
	startPoint[1] = pt1[1];
	startPoint[2] = pt1[2];
	startPoint[3] = 1;

	std::cout << "startpoint " << startPoint[0] << ", " << startPoint[1] << ", " << startPoint[2] << std::endl;

	IndexType endPoint;
	endPoint[0] = pt2[0];
	endPoint[1] = pt2[1];
	endPoint[2] = pt2[2];
	endPoint[3] = 1;

	std::cout << "endpoint " << endPoint[0] << ", " << endPoint[1] << ", " << endPoint[2] << std::endl;

	pathFilter->SetStartPoint(startPoint);
	pathFilter->AddPathEndPoint(endPoint);
	//pathFilter->SetRegionToProcess(region);

	std::cout<< "path update" << std::endl;

	pathFilter->Update();
	pathFilter->WritePathsToFile("testPath");

	env->ReleaseFloatArrayElements(jba, jbs,0);
	env->ReleaseFloatArrayElements(Path, output,0);

	return 0;
}


JNIEXPORT jint JNICALL Java_FijiITKInterface_MultiScaleTubularityMeasure_OrientedFlux(JNIEnv *env, jobject ignored, jbyteArray jba, jfloatArray jbOut, jint type, jint width, jint height, jint NSlice, jdouble widthpix, jdouble heightpix, jdouble depthpix, jdouble sigmaMin, jdouble sigmaMax, jint numberOfScales, jfloatArray pt1, jfloatArray pt2)
{
    jboolean isCopy;
    jbyte * jbs   = env->GetByteArrayElements(jba,&isCopy);
    jfloat * jbOutS = env->GetFloatArrayElements(jbOut,&isCopy);

    jfloat * jSeed = env->GetFloatArrayElements(pt1,&isCopy);
    jfloat * jEndPoint = env->GetFloatArrayElements(pt2,&isCopy);

    float * Seed = (float*) jSeed;
    float * EndPoint = (float*) jEndPoint;

    if( ! jbs )
        return -1;

	double spacing[4], origin[4]; 

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

	spacing[0] = widthpix;spacing[1] = heightpix;spacing[2] = depthpix; spacing[3] = (sigmaMax-sigmaMin)/(numberOfScales-1) ;
	itkImageP->SetSpacing( spacing );

	origin[0] = 0;origin[1] = 0;origin[2] = 0; origin[3] = sigmaMin;
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

	OutputImageType::Pointer outputImage = Execute<unsigned char, 3>(itkImageP, sigmaMin, sigmaMax, numberOfScales, Seed, EndPoint);

	OutputImageType::RegionType Outputregion;
	OutputImageType::SizeType Outputsize;
	Outputsize[0] = width;Outputsize[1] = height;Outputsize[2] = NSlice; 
	OutputImageType::IndexType Outputstart;
	Outputstart[0] = 0;Outputstart[1] = 0;Outputstart[2] = 0;
	Outputregion.SetSize( Outputsize );
	Outputregion.SetIndex( Outputstart );

	outputImage->SetRegions( Outputregion);
	outputImage->Allocate();

	spacing[0] = widthpix;spacing[1] = heightpix;spacing[2] = depthpix;
	outputImage->SetSpacing( spacing );

	origin[0] = 0;origin[1] = 0;origin[2] = 0;
	outputImage->SetOrigin( origin );

	float* outputImageData = (float*) jbOutS;
	OutputIteratorType outit( outputImage, Outputregion);

	outit.GoToBegin();
	int length, ind = 0; double w = width; double h = height;

	for(int k=0;k<numberOfScales;k++){
		length =  w * h * NSlice;
		for( int i = 0; i < length; ++i ) {
	  		outputImageData[ind] = outit.Get();
			ind++;
	  		++outit;
		}
	}

    env->ReleaseByteArrayElements(jba,jbs,0);
    env->ReleaseFloatArrayElements(jbOut, jbOutS,0);

    return 0;
}

void setSuccess(JNIEnv * env,
                bool succeeded,
                jobject pathResultObject,
                jclass pathResultClass)
{
    jmethodID mid = env->GetMethodID(pathResultClass,
                                     "setSuccess",
                                     "(Z)V");
    if (!mid) {
        std::cout << "Failed to find the setSuccess method";
        return;
    }

    jboolean jSucceeded = succeeded ? 1 : 0;

    env->CallVoidMethod(pathResultObject,
                        mid,
                        jSucceeded);
}

void setErrorMessage(JNIEnv * env,
                     const char * message,
                     jobject pathResultObject,
                     jclass pathResultClass)
{
    jmethodID mid = env->GetMethodID(pathResultClass,
                                     "setErrorMessage",
                                     "(Ljava/lang/String;)V");
    if (!mid) {
        std::cout << "Failed to find the setErrorMessage method";
        return;
    }

    jstring jMessage = env->NewStringUTF(message);

    env->CallVoidMethod(pathResultObject,
                        mid,
                        jMessage);

    env->DeleteLocalRef(jMessage);

    setSuccess(env,
               false,
               pathResultObject,
               pathResultClass);
}

/*
 * Class:     FijiITKInterface_MultiScaleTubularityMeasure
 * Method:    getPathResult
 * Signature: (Ljava/lang/String;[F[FLFijiITKInterface/PathResult;)V
 */
JNIEXPORT void JNICALL Java_FijiITKInterface_MultiScaleTubularityMeasure_getPathResult
 (JNIEnv * env,
  jobject ignored,
  jstring jTubularityFilename,
  jfloatArray jPoint1,
  jfloatArray jPoint2,
  jobject pathResultObject)
{

    jclass pathResultClass = env->GetObjectClass(pathResultObject);

    // Check that the float arrays are of the right length:

    int pt1Length = env->GetArrayLength(jPoint1);
    int pt2Length = env->GetArrayLength(jPoint2);

    if (pt1Length != 3) {
        setErrorMessage(env,
                        "pt1 was not of length 3",
                        pathResultObject,
                        pathResultClass);
        return;
    }

    if (pt2Length != 3) {
        setErrorMessage(env,
                        "pt2 was not of length 3",
                        pathResultObject,
                        pathResultClass);
        return;
    }

    // And convert them to C types:

    jboolean isPoint1Copy, isPoint2Copy;
    jfloat * pt1 = env->GetFloatArrayElements(jPoint1, &isPoint1Copy);
    jfloat * pt2 = env->GetFloatArrayElements(jPoint2, &isPoint2Copy);

    // Some debugging, just print out the start and end points:

    for(int i = 0; i < 3; ++i ) {
        std::cout << "pt1[" << i << "] is " << pt1[i] << std::endl;
    }
    for(int i = 0; i < 3; ++i ) {
        std::cout << "pt2[" << i << "] is " << pt2[i] << std::endl;
    }

    // Now get the filename:

    const char * filename = env->GetStringUTFChars(jTubularityFilename, NULL);
    if (!filename) {
        setErrorMessage(env,
                        "Failed to convert the filename",
                        pathResultObject,
                        pathResultClass);
        env->ReleaseFloatArrayElements(jPoint1, pt1, 0);
        env->ReleaseFloatArrayElements(jPoint2, pt2, 0);
        return;
    }

    // ------------------------------------------------------------------------
    // FIXME: real code should go here, just example data:

    int points = 40;
    float * pathPoints = new float[points*4];

    for(int i = 0; i < points; ++i ) {
        pathPoints[i] = i;
        pathPoints[i+1] = i;
        pathPoints[i+2] = i;
        pathPoints[i+3] = 6 + ((i % 8) - 4);
    }

    // Now convert that to a Java float array:

    jfloatArray jResultArray = env->NewFloatArray(points * 4);
    if (!jResultArray) {
        std::cout << "Failed to allocate a new Java float array";
        env->ReleaseFloatArrayElements(jPoint1, pt1, 0);
        env->ReleaseFloatArrayElements(jPoint2, pt2, 0);
        env->ReleaseStringUTFChars(jTubularityFilename, filename);
        return;
    }

    env->SetFloatArrayRegion(jResultArray,
                             0,
                             points * 4,
                             (jfloat *)pathPoints);

    {
        jmethodID mid = env->GetMethodID(pathResultClass,
                                         "setPath",
                                         "([F)V");
        if (!mid) {
            std::cout << "Failed to find the setPath method";
            env->ReleaseFloatArrayElements(jPoint1, pt1, 0);
            env->ReleaseFloatArrayElements(jPoint2, pt2, 0);
            env->ReleaseStringUTFChars(jTubularityFilename, filename);
            return;
        }

        env->CallVoidMethod(pathResultObject,
                            mid,
                            jResultArray);
    }

    env->DeleteLocalRef(jResultArray);

    delete [] pathPoints;

    // ------------------------------------------------------------------------

    setSuccess(env,
               true,
               pathResultObject,
               pathResultClass);

    /* Release anything allocated... */

    env->ReleaseFloatArrayElements(jPoint1, pt1, 0);
    env->ReleaseFloatArrayElements(jPoint2, pt2, 0);
    env->ReleaseStringUTFChars(jTubularityFilename, filename);
}
