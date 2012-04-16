/* A trivial C++ program for testing JNI calls */

#include <iostream>

#include "FijiITKInterface_TubularGeodesics.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include "itkImageFileReader.h"
#include "itkTubularMetricToPathFilter.h"
#include <itkMultiThreader.h>
#include <itkSemaphore.h>
#include <itkFastMutexLock.h>
#include "vnl/vnl_math.h"
#include <jni.h>
#include <unistd.h>

using std::cout;
using std::endl;
using std::flush;

// Consts and typedefs
const unsigned int Dimension = 3;
const unsigned int SSDimension = Dimension+1;
typedef float																												 TubularityScorePixelType;
typedef itk::Image<TubularityScorePixelType,4>	                     TubularityScoreImageType;
typedef TubularityScoreImageType::RegionType                         RegionType;
typedef TubularityScoreImageType::SizeType                           SizeType;

typedef TubularityScoreImageType::IndexValueType                     IndexValueType;
typedef TubularityScoreImageType::IndexType                          IndexType;
typedef TubularityScoreImageType::PointType                          OriginType;
typedef TubularityScoreImageType::SpacingType                        SpacingType;

typedef itk::TubularMetricToPathFilter< TubularityScoreImageType >  PathFilterType;
typedef PathFilterType::VertexType			             VertexType;

typedef itk::ImageFileReader< TubularityScoreImageType >             ImageReaderType;

// Global variables
TubularityScoreImageType::Pointer tubularityScore;
bool isTubularityScoreLoaded = false;
std::vector< float > Outputpath;

JavaVM * globalJVM = NULL;
jobject pathResultObject;
jobject javaSearchThread;
itk::FastMutexLock::Pointer globalMutex = itk::FastMutexLock::New();
bool pleaseStop;

JNIEXPORT void JNICALL
reportProgress(JNIEnv *env, jobject obj, jfloat proportionDone)
{
  jclass cls = env->GetObjectClass(obj);
  jmethodID mid = env->GetMethodID(cls, "reportProgress", "(F)V");
  if (! mid) {
      cout << "Failed to find the reportProgress method" << endl;
      return;
  }
  env->CallVoidMethod(obj, mid, proportionDone);
}

JNIEXPORT void JNICALL
reportFinished(JNIEnv *env, jobject obj, jboolean success)
{
  jclass cls = env->GetObjectClass(obj);
  jmethodID mid = env->GetMethodID(cls, "reportFinished", "(Z)V");
  if (! mid) {
      cout << "Failed to find the reportFinished method" << endl;
      return;
  }
  env->CallVoidMethod(obj, mid, success);
}

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

enum ExecuteReturnValues {
    eSuccess = 0,
    eInterrupted = 1,
    eFailed = 2
};

// Computes the minimal path between 2 provided points
/**
 * ITK related methods: 
 * 
 */
int
Execute( float* pt1, float* pt2)
{

  // Instantiate the path filter
  PathFilterType::Pointer pathFilter = PathFilterType::New();
  
  // Set the tubularity score
  pathFilter->SetInput( tubularityScore );

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
  //TODO: This shouldn't be hardcoded
  int subRegionPad    = 50;
  //TODO: This shouldn't be hardcoded
  // Anyways, we'll get rid of it soon ;-)
  for(unsigned int i = 0; i < Dimension; i++)
    {
      IndexValueType minIndex = vnl_math_min( startPoint[i], endPoint[i] );
      IndexValueType maxIndex = vnl_math_max( startPoint[i], endPoint[i] );
      startSubRegion[i] = vnl_math_max( minIndex - subRegionPad, region.GetIndex()[i] );
      IndexValueType maxSubRegionIndex = vnl_math_min( int(maxIndex + subRegionPad), int(region.GetIndex()[i] + region.GetSize()[i] -1) );
      sizeSubRegion[i]  = maxSubRegionIndex - startSubRegion[i] + 1;
    }
  subRegionToProcess.SetIndex( startSubRegion );
  subRegionToProcess.SetSize( sizeSubRegion );

  pathFilter->SetRegionToProcess(subRegionToProcess);
  try {
      pathFilter->Update();
  } catch (itk::ProcessAborted &e) {
      return eInterrupted;
  }
  
	SpacingType spacing = tubularityScore->GetSpacing();
  OriginType  origin = tubularityScore->GetOrigin();
	
	// Get the minimum spacing among all the spatial dimensions.
	double minSpacing = spacing[0];
	for(unsigned int i = 1; i < Dimension-1; i++)
	{
		minSpacing = vnl_math_min(minSpacing, spacing[i]);
	}
	
	// Downsample the path and smooth it slightly.
	pathFilter->GetPath(0)->Resample(0.5 * minSpacing, tubularityScore.GetPointer());
	pathFilter->GetPath(0)->SmoothVertexLocationsAndRadii(minSpacing, tubularityScore.GetPointer());

  Outputpath.clear();
  for(unsigned int k = 0; k < pathFilter->GetPath(0)->GetVertexList()->Size(); k++)
	{
		VertexType vertex = pathFilter->GetPath(0)->GetVertexList()->GetElement(k);
		for (unsigned int i = 0; i < Dimension+1; i++) 
	  {
	    Outputpath.push_back(vertex[i]*spacing[i]+origin[i]);
	  }
	}
  return eSuccess;
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
        cout << "Failed to find the setSuccess method" << endl;
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
        cout << "Failed to find the setErrorMessage method" << endl;
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

class UserData {

public:
    jstring jTubularityFilename;
    jfloatArray jPoint1;
    jfloatArray jPoint2;
    itk::Semaphore * semaphoreUserDataCopied;
};

void releaseJVM(UserData * userData) {
    delete userData;
    globalJVM->DetachCurrentThread();
    globalMutex->Lock();
    globalJVM = NULL;
    globalMutex->Unlock();
}

ITK_THREAD_RETURN_TYPE ThreadedFunction(void* param) {

    JNIEnv * env;
    UserData * userData = NULL;

    // From: http://www.adamish.com/blog/archives/327

    int getEnvStat = globalJVM->GetEnv((void **)&env, JNI_VERSION_1_6);
    if (getEnvStat == JNI_EDETACHED) {
        if (globalJVM->AttachCurrentThread((void **)&env, NULL) != 0) {
            cout << "Failed to attach to JVM" << endl;
            return ITK_THREAD_RETURN_VALUE;
        }
    } else if (getEnvStat == JNI_OK) {
    } else if (getEnvStat == JNI_EVERSION) {
        cout << "Failed to attach to JVM with GetEnv: version not supported" << endl;
    }

    itk::MultiThreader::ThreadInfoStruct* threadInfo = static_cast<itk::MultiThreader::ThreadInfoStruct*>(param);

    if (!threadInfo)  {
        cout << "Failed to get the thread info" << endl;
        releaseJVM(userData);
        return ITK_THREAD_RETURN_VALUE;
    }

    userData = static_cast<UserData*>(threadInfo->UserData);

    // Generate global references for the two objects that were passed in:

    jclass pathResultClass = env->GetObjectClass(pathResultObject);

    // Check that the float arrays are of the right length:

    int pt1Length = env->GetArrayLength(userData->jPoint1);
    int pt2Length = env->GetArrayLength(userData->jPoint2);

    if (pt1Length != 3) {
        cout << "wrong length of pt1" << endl;

        setErrorMessage(env,
                        "pt1 was not of length 3",
                        pathResultObject,
                        pathResultClass);
        releaseJVM(userData);
        return ITK_THREAD_RETURN_VALUE;
    }

    if (pt2Length != 3) {
        cout << "wrong length of pt2" << endl;
        setErrorMessage(env,
                        "pt2 was not of length 3",
                        pathResultObject,
                        pathResultClass);
        releaseJVM(userData);
        return ITK_THREAD_RETURN_VALUE;
    }

    // And convert them to C types:

    jboolean isPoint1Copy, isPoint2Copy;
    jfloat * pt1 = env->GetFloatArrayElements(userData->jPoint1, &isPoint1Copy);
    jfloat * pt2 = env->GetFloatArrayElements(userData->jPoint2, &isPoint2Copy);

    float copiedPt1[3], copiedPt2[3];

    for(int i = 0; i < 3; ++i ) {
        copiedPt1[i] = pt1[i];
        copiedPt2[i] = pt2[i];
    }

    env->ReleaseFloatArrayElements(userData->jPoint1, pt1, 0);
    env->ReleaseFloatArrayElements(userData->jPoint2, pt2, 0);

    // Now get the filename:

    const char * filename = env->GetStringUTFChars(userData->jTubularityFilename, NULL);

    if (!filename) {
        setErrorMessage(env,
                        "Failed to convert the filename",
                        pathResultObject,
                        pathResultClass);

        releaseJVM(userData);
        return ITK_THREAD_RETURN_VALUE;
    }

    // Now we can safely allow the function that spawned this thread
    // to return:

    userData->semaphoreUserDataCopied->Up();

    // ------------------------------------------------------------------------

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
          std::cerr << e << endl;
          env->ReleaseStringUTFChars(userData->jTubularityFilename, filename);
          setErrorMessage(env,
                          e.GetDescription(),
                          pathResultObject,
                          pathResultClass);

          releaseJVM(userData);
          return ITK_THREAD_RETURN_VALUE;
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
    try {
        int executeResult = Execute( copiedPt1, copiedPt2 );
        if (eInterrupted == executeResult) {
            // ... then interrupt

            /*
            setSuccess(env,
                       true,
                       pathResultObject,
                       pathResultClass);
            */

            reportFinished(env,
                           javaSearchThread,
                           false);
            releaseJVM(userData);
            return ITK_THREAD_RETURN_VALUE;
        }
    } catch(itk::ExceptionObject &e) {
        env->ReleaseStringUTFChars(userData->jTubularityFilename, filename);
        setErrorMessage(env,
                        e.GetDescription(),
                        pathResultObject,
                        pathResultClass);
        releaseJVM(userData);
        return ITK_THREAD_RETURN_VALUE;
    }

    long nb_points  = Outputpath.size() / 4;
    float * pathPoints = new float[nb_points*4];
    for(int i = 0; i < 4*nb_points; ++i ) 
      {
	pathPoints[i] = Outputpath[i];
      }

    // Now convert that to a Java float array:

    jfloatArray jResultArray = env->NewFloatArray(nb_points * 4);
    if (!jResultArray) {
        cout << "Failed to allocate a new Java float array" << endl;
        env->ReleaseStringUTFChars(userData->jTubularityFilename, filename);
        releaseJVM(userData);
        return ITK_THREAD_RETURN_VALUE;
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
            cout << "Failed to find the setPath method" << endl;
            env->ReleaseStringUTFChars(userData->jTubularityFilename, filename);
            releaseJVM(userData);
            return ITK_THREAD_RETURN_VALUE;
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

    reportFinished(env,
                   javaSearchThread,
                   true);

    env->ReleaseStringUTFChars(userData->jTubularityFilename, filename);

    /* Now we can delete the global references to the two objects that
       were passed in: */
    env->DeleteGlobalRef(pathResultObject);
    env->DeleteGlobalRef(javaSearchThread);

    releaseJVM(userData);
    return ITK_THREAD_RETURN_VALUE;
}


/*
 * Class:     FijiITKInterface_TubularGeodesics
 * Method:    interruptSearch
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_FijiITKInterface_TubularGeodesics_interruptSearch
  (JNIEnv *, jobject)
{
    globalMutex->Lock();
    pleaseStop = true;
    globalMutex->Unlock();
}

/*
 * Class:     FijiITKInterface_TubularGeodesics
 * Method:    startSearch
 * Signature: (Ljava/lang/String;[F[FLtracing/PathResult;Ltracing/TubularGeodesicsTracer;)V
 */
JNIEXPORT void JNICALL Java_FijiITKInterface_TubularGeodesics_startSearch
 (JNIEnv * env,
  jobject ignored,
  jstring jTubularityFilename,
  jfloatArray jPoint1,
  jfloatArray jPoint2,
  jobject passedPathResultObject,
  jobject passedJavaSearchThread)
{
    globalMutex->Lock();
    if (globalJVM) {
        reportFinished(env, passedJavaSearchThread, false);
        globalMutex->Unlock();
        env->DeleteGlobalRef(pathResultObject);
        env->DeleteGlobalRef(javaSearchThread);
        return;
    } else {
        env->GetJavaVM(&globalJVM);
        globalMutex->Unlock();
    }

    /* Keep these two in global variables, now that we know we're the
       only thread running (due to checking jvm under a mutex): */
    pathResultObject = env->NewGlobalRef(passedPathResultObject);
    javaSearchThread = env->NewGlobalRef(passedJavaSearchThread);

    globalMutex->Lock();
    pleaseStop = false;
    globalMutex->Unlock();

    /* The code here that spawns this as a separate thread is based on
       the example here:
       http://docs.mitk.org/nightly-qt4/mitkITKThreadingTest_8cpp_source.html
    */

    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    itk::Semaphore::Pointer semaphorePointer = itk::Semaphore::New();
    semaphorePointer->Initialize(0);

    UserData * userData = new UserData();
    userData->jTubularityFilename = jTubularityFilename;
    userData->jPoint1 = jPoint1;
    userData->jPoint2 = jPoint2;
    userData->semaphoreUserDataCopied = semaphorePointer.GetPointer();

    itk::ThreadFunctionType functionPointer = &ThreadedFunction;
    threader->SpawnThread( functionPointer, userData );

    // Wait until the userData has been copied by the other thread:
    semaphorePointer->Down();

    /* Now we know that the userData has been copied, and the thread
       is running, so we can leave the function. */
}
