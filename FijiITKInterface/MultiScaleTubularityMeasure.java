package FijiITKInterface;

import fiji.jni.LibraryLoader;
import tracing.PathResult;

public class MultiScaleTubularityMeasure extends LibraryLoader {

    public native int GetPath(float [] jba);
    public native int FindPath(float [] Input, int [] pt1, int [] pt2 , float [] Path, int w, int h, int slices,int scales,  double pixw, double pixh, double pixd);
    public native int OrientedFlux(byte [] imageIn,float [] imageOut, int type, int width, int height, int Nslice, double pixwidth, double pixheight, double pixdepth, double sigmaMin, double sigmaMax, int scales, float [] pt1, float[] pt2);

    public native void getPathResult(String tubularityFilename,
                                     float [] p1,
                                     float [] p2,
                                     PathResult result);

}
