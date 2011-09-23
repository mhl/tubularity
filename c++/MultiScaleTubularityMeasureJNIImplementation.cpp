/* A trivial C++ program for testing JNI calls */

#include <iostream>

#include "FijiITKInterface_MultiScaleTubularityMeasure.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include "itkImageFileReader.h"
#include "itkTubularMetricToPathFilter3.h"
#include "vnl/vnl_math.h"


// Consts and typedefs
const unsigned int Dimension = 3;
const unsigned int SSDimension = Dimension+1;
typedef float	                                                     TubularityScorePixelType;
typedef itk::Image<TubularityScorePixelType,4>	                     TubularityScoreImageType;
typedef TubularityScoreImageType::RegionType                         RegionType;
typedef TubularityScoreImageType::SizeType                           SizeType;

typedef TubularityScoreImageType::IndexValueType                     IndexValueType;
typedef TubularityScoreImageType::IndexType                          IndexType;
typedef TubularityScoreImageType::PointType                          OriginType;
typedef TubularityScoreImageType::SpacingType                        SpacingType;

typedef itk::TubularMetricToPathFilter3< TubularityScoreImageType >  PathFilterType;
typedef PathFilterType::VertexType			             VertexType;

typedef itk::ImageFileReader< TubularityScoreImageType >             ImageReaderType;

// Global variables
TubularityScoreImageType::Pointer tubularityScore;
bool isTubularityScoreLoaded = false;
std::vector< float > Outputpath;


// For a Given Location, Get the Optimal scale, 
/**
 * ITK related methods: 
 * 
 */
void GetOptimalScale(IndexType *point)
{
  TubularityScorePixelType bestScore    =  itk::NumericTraits< TubularityScorePixelType >::min();
  IndexType scaleSpaceSourceVertexIndex = *point;
  RegionType region                     = tubularityScore->GetBufferedRegion();
  IndexValueType noOfScales             = region.GetSize()[Dimension];	       
  IndexValueType scaleStartIndex        = region.GetIndex()[Dimension];	       
  IndexValueType scaleEndIndex          = scaleStartIndex + noOfScales - 1; 
  IndexValueType bestScaleIndex         = 0;
  for(IndexValueType sourceScaleIndex   = scaleStartIndex;
     sourceScaleIndex  <= scaleEndIndex;       			      
     sourceScaleIndex++)
    {
      scaleSpaceSourceVertexIndex[Dimension] = sourceScaleIndex;
      if( bestScore < tubularityScore->GetPixel( scaleSpaceSourceVertexIndex ) )
	{
	  bestScore = tubularityScore->GetPixel( scaleSpaceSourceVertexIndex );				       
	  bestScaleIndex = sourceScaleIndex;
	}
    }
  
  (*point)[Dimension] = bestScaleIndex;  
}

// Computes the minimal path between 2 provided points
/**
 * ITK related methods: 
 * 
 */
void
Execute( float* pt1, float* pt2)
{

  // Instantiate the path filter
  PathFilterType::Pointer pathFilter = PathFilterType::New();
  
  // Set the tubularity score
  pathFilter->SetInput( tubularityScore );
  pathFilter->SetScaleSpeedFactor(10.0);// TODO: shouldn't be hardcoded

  // Get the start and end points and give them to the path filter
  IndexType startPoint;
  IndexType endPoint;
  for(unsigned int i = 0; i < Dimension; i++)
    {
      startPoint[i] = pt1[i];
      endPoint[i]   = pt2[i];
    }
  // Get and assign to them the optimal scale
  GetOptimalScale( &startPoint );
  GetOptimalScale( &endPoint );
  pathFilter->SetStartPoint( startPoint );
  pathFilter->AddPathEndPoint( endPoint );
  
  // Get the sub region to be processed
  // Warning a padding parameter is hardcoded
  RegionType region = tubularityScore->GetBufferedRegion();
 
  RegionType subRegionToProcess;
  IndexType startSubRegion;
  SizeType  sizeSubRegion;
  
  // No Padding or sub-selcetion on the scale dimension
  startSubRegion[Dimension] = region.GetIndex()[Dimension];
  sizeSubRegion[Dimension]  = region.GetSize()[Dimension];
  // extract sub region and pad it in the spatial domain
  int subRegionPad    = 2;
  for(unsigned int i = 0; i < Dimension; i++)
    {
      startSubRegion[i] = vnl_math_max(vnl_math_min( startPoint[i], endPoint[i] ) - subRegionPad, region.GetIndex()[i] );
      sizeSubRegion[i]  = vnl_math_min(int(vnl_math_max( startPoint[i], endPoint[i] ) - startSubRegion[i] + subRegionPad+1), int(region.GetSize()[i] - startSubRegion[i] + 1) );//TODO
    }
  subRegionToProcess.SetIndex( startSubRegion );
  subRegionToProcess.SetSize( sizeSubRegion );

  pathFilter->SetRegionToProcess(subRegionToProcess);

  pathFilter->Update();
  
  SpacingType spacing = tubularityScore->GetSpacing();
  OriginType  origin = tubularityScore->GetOrigin();

  Outputpath.clear();
  for(unsigned int k = 0; k < pathFilter->GetPath(0)->GetVertexList()->Size(); k++)
    {
       	VertexType vertex = pathFilter->GetPath(0)->GetVertexList()->GetElement(k);
	for (unsigned int i = 0; i < Dimension+1; i++) 
	  {
	    Outputpath.push_back(vertex[i]*spacing[i]+origin[i]);
	  }
    }
}

/**
 * JNI related methods: 
 * 
 */
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

/**
 * JNI related methods: 
 * 
 */
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
    /**
     * Fethallah 
     * First, check if the tubularity score image is loaded.
     * If not, load it
     */
    if( !isTubularityScoreLoaded )
      {
	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetFileName( filename );
	try
	  {
	    reader->Update();
	  }
	catch(itk::ExceptionObject &e)   
	  {      
	    std::cerr << e << std::endl; 
	  }

	tubularityScore = reader->GetOutput();
	tubularityScore->DisconnectPipeline();
	isTubularityScoreLoaded = true;
      }
    /**
     * Fethallah 
     * At this point, the tubularity score is supposed to be loaded and
     * the start and end points are given
     * One just needs to call the Execute method and convert the output
     */

    Execute( pt1, pt2 );

    long nb_points  = Outputpath.size() / 4;
    float * pathPoints = new float[nb_points*4];
    for(int i = 0; i < nb_points; ++i ) 
      {
	pathPoints[i  ] = Outputpath[4*i + 0];
        pathPoints[i+1] = Outputpath[4*i + 1];
        pathPoints[i+2] = Outputpath[4*i + 2];
        pathPoints[i+3] = Outputpath[4*i + 3];
      }
    /*
    // ------------------------------------------------------------------------
    // FIXME: real code should go here, just example data:

    int points = 40;
    float * pathPoints = new float[points*4];

    for(int i = 0; i < points; ++i ) {
        pathPoints[i] = i;
        pathPoints[i+1] = i;
        pathPoints[i+2] = i;
        pathPoints[i+3] = 6 + ((i % 8) - 4);
    }*/

    // Now convert that to a Java float array:

    jfloatArray jResultArray = env->NewFloatArray(nb_points * 4);
    if (!jResultArray) {
        std::cout << "Failed to allocate a new Java float array";
        env->ReleaseFloatArrayElements(jPoint1, pt1, 0);
        env->ReleaseFloatArrayElements(jPoint2, pt2, 0);
        env->ReleaseStringUTFChars(jTubularityFilename, filename);
        return;
    }

    env->SetFloatArrayRegion(jResultArray,
                             0,
                             nb_points * 4,
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
