/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class FijiITKInterface_TubularGeodesics */

#ifndef _Included_FijiITKInterface_TubularGeodesics
#define _Included_FijiITKInterface_TubularGeodesics
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     FijiITKInterface_TubularGeodesics
 * Method:    GetPath
 * Signature: ([F)I
 */
JNIEXPORT jint JNICALL Java_FijiITKInterface_TubularGeodesics_GetPath
  (JNIEnv *, jobject, jfloatArray);

/*
 * Class:     FijiITKInterface_TubularGeodesics
 * Method:    FindPath
 * Signature: ([F[I[I[FIIIIDDD)I
 */
JNIEXPORT jint JNICALL Java_FijiITKInterface_TubularGeodesics_FindPath
  (JNIEnv *, jobject, jfloatArray, jintArray, jintArray, jfloatArray, jint, jint, jint, jint, jdouble, jdouble, jdouble);

/*
 * Class:     FijiITKInterface_TubularGeodesics
 * Method:    OrientedFlux
 * Signature: ([B[FIIIIDDDDDI[F[F)I
 */
JNIEXPORT jint JNICALL Java_FijiITKInterface_TubularGeodesics_OrientedFlux
  (JNIEnv *, jobject, jbyteArray, jfloatArray, jint, jint, jint, jint, jdouble, jdouble, jdouble, jdouble, jdouble, jint, jfloatArray, jfloatArray);

/*
 * Class:     FijiITKInterface_TubularGeodesics
 * Method:    getPathResult
 * Signature: (Ljava/lang/String;[F[FLtracing/PathResult;)V
 */
JNIEXPORT void JNICALL Java_FijiITKInterface_TubularGeodesics_getPathResult
  (JNIEnv *, jobject, jstring, jfloatArray, jfloatArray, jobject);

#ifdef __cplusplus
}
#endif
#endif