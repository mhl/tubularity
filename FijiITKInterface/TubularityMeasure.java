package FijiITKInterface;

import fiji.jni.LibraryLoader;

public class TubularityMeasure extends LibraryLoader {

    public native int OrientedFlux(byte [] imageIn,float [] imageOut, int type, int width, int height, int Nslice, double pixwidth, double pixheight, double pixdepth, double sigmaMin, double sigmaMax, int scales, String outputFilename);
}

