.PHONY : all clean test
.PRECIOUS : build/linux64/libTubularGeodesics.so \
	build/linux64/libOOFTubularityMeasure.so \
	build/linux/libTubularGeodesics.so \
	build/linux/libOOFTubularityMeasure.so \
	FijiITKInterface/FijiITKInterface_TubularGeodesics.h \
	FijiITKInterface/FijiITKInterface_OOFTubularityMeasure.h

# Ensure the JDK_INCLUDE_PATH is defined:
#
# Examples:
# On Linux it may well be:
#   /usr/lib/jvm/java-6-openjdk/
# On Mac OS:
#   /usr/lib/jvm/java-1.6.0-openjdk/include/
# On Windows:
#   [not sure yet]

ifndef JDK_HOME
$(error The environment variable JDK_HOME must be defined.)
endif

# You must also set the environment variable FIJI_BINARY to point
# to a working instance of the Fiji launcher.

ifndef FIJI_LAUNCHER
$(info The environment variable FIJI_LAUNCHER must be defined.)
$(error FIJI_LAUNCHER undefined)
endif

# Finally, you have to set up the environment variable ITK to point to
# the root of an ITK installation.  This expects libraries to then be
# under $ITK/bin.  If that is incorrect, you can override it by
# setting the environment variable ITK_LIBS.

ifndef ITK
$(info The environment variable ITK must be defined, and point)
$(info to the root of an ITK installation, in this case where)
$(info itkConfigure.h is found.  The ITK libraries are expected)
$(info to be in $ITK/bin - if that is not the case, you should)
$(info also define ITK_LIBS to point to their location.)
$(error ITK undefined)
endif

ifndef ITK_LIBS
ITK_LIBS=$(ITK)/bin/
endif

ifndef FFTW_PREFIX
FFTW_PREFIX=/usr/lib/
endif

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')

ifeq ($(ARCH),)
ifeq ($(uname_S),Linux)
JAVA_ARCH_NAME=linux
	LIBRARY_EXTENSION=so
ifeq ($(uname_M),x86_64)
         ARCH=linux64
else
         ARCH=linux
endif
endif
# FIXME: doesn't detect win64 at the moment
ifneq (,$(findstring MINGW,$(uname_S)))
        ARCH=win32
endif
ifeq ($(uname_S),Darwin)
        ARCH=macosx
	LIBRARY_EXTENSION=dylib
endif
endif

ifndef ARCH
$(info The architecture could not be automatically discovered.)
$(info Please report this as a bug to mark-imagej@longair.net,)
$(info including the following output:)
$(info uname_S is: $(uname_S))
$(info uname_M is: $(uname_M))
$(error)
endif


## 	-luuid \ TODO

# -framework CoreFoundation : works for macOS but not linux
	
LINK_LIBRARIES_ITK=-L$(ITK_LIBS) \
 $(ITK_LIBS)libitksys-4.1.a \
 $(ITK_LIBS)libitkvnl_algo-4.1.a \
 $(ITK_LIBS)libitkvnl-4.1.a \
 $(ITK_LIBS)libitkv3p_netlib-4.1.a \
 $(ITK_LIBS)libITKCommon-4.1.a \
 $(ITK_LIBS)libitkNetlibSlatec-4.1.a \
 $(ITK_LIBS)libITKStatistics-4.1.a \
 $(ITK_LIBS)libITKIOImageBase-4.1.a \
 $(ITK_LIBS)libITKMesh-4.1.a \
 $(ITK_LIBS)libitkzlib-4.1.a \
 $(ITK_LIBS)libITKMetaIO-4.1.a \
 $(ITK_LIBS)libITKSpatialObjects-4.1.a \
 $(ITK_LIBS)libITKPath-4.1.a \
 $(ITK_LIBS)libITKLabelMap-4.1.a \
 $(ITK_LIBS)libITKQuadEdgeMesh-4.1.a \
 $(ITK_LIBS)libITKOptimizers-4.1.a \
 $(ITK_LIBS)libITKPolynomials-4.1.a \
 $(ITK_LIBS)libITKBiasCorrection-4.1.a \
 $(ITK_LIBS)libITKBioCell-4.1.a \
 $(ITK_LIBS)libITKFFT-4.1.a \
 $(ITK_LIBS)libITKDICOMParser-4.1.a \
 $(ITK_LIBS)libITKEXPAT-4.1.a \
 $(ITK_LIBS)libITKIOXML-4.1.a \
 $(ITK_LIBS)libITKIOSpatialObjects-4.1.a \
 $(ITK_LIBS)libITKFEM-4.1.a \
 $(ITK_LIBS)libITKIOBMP-4.1.a \
 $(ITK_LIBS)libITKIOBioRad-4.1.a \
 $(ITK_LIBS)libitkopenjpeg-4.1.a \
 $(ITK_LIBS)libITKIOIPL-4.1.a \
 $(ITK_LIBS)libITKIOGE-4.1.a \
 $(ITK_LIBS)libITKIOGIPL-4.1.a \
 $(ITK_LIBS)libitkjpeg-4.1.a \
 $(ITK_LIBS)libITKIOJPEG-4.1.a \
 $(ITK_LIBS)libitktiff-4.1.a \
 $(ITK_LIBS)libITKIOTIFF-4.1.a \
 $(ITK_LIBS)libITKIOLSM-4.1.a \
 $(ITK_LIBS)libITKIOMeta-4.1.a \
 $(ITK_LIBS)libITKznz-4.1.a \
 $(ITK_LIBS)libITKniftiio-4.1.a \
 $(ITK_LIBS)libITKIONIFTI-4.1.a \
 $(ITK_LIBS)libITKNrrdIO-4.1.a \
 $(ITK_LIBS)libITKIONRRD-4.1.a \
 $(ITK_LIBS)libitkpng-4.1.a \
 $(ITK_LIBS)libITKIOPNG-4.1.a \
 $(ITK_LIBS)libITKIOSiemens-4.1.a \
 $(ITK_LIBS)libITKIOStimulate-4.1.a \
 $(ITK_LIBS)libITKIOVTK-4.1.a \
 $(ITK_LIBS)libITKKLMRegionGrowing-4.1.a \
 $(ITK_LIBS)libITKVTK-4.1.a \
 $(ITK_LIBS)libITKWatersheds-4.1.a \
 $(ITK_LIBS)libITKDeprecated-4.1.a \
 $(ITK_LIBS)libITKgiftiio-4.1.a \
 $(ITK_LIBS)libitkhdf5_cpp-4.1.a \
 $(ITK_LIBS)libitkhdf5-4.1.a \
 $(ITK_LIBS)libITKIOCSV-4.1.a \
 $(ITK_LIBS)libITKIOHDF5-4.1.a \
 $(ITK_LIBS)libITKIOMesh-4.1.a \
 $(ITK_LIBS)libITKIOTransformBase-4.1.a \
 $(ITK_LIBS)libITKIOTransformHDF5-4.1.a \
 $(ITK_LIBS)libITKIOTransformInsightLegacy-4.1.a \
 $(ITK_LIBS)libITKIOTransformMatlab-4.1.a \
 $(ITK_LIBS)libITKOptimizersv4-4.1.a \
 $(ITK_LIBS)libITKReview-4.1.a \
 $(ITK_LIBS)libITKVideoCore-4.1.a \
 $(ITK_LIBS)libITKVideoIO-4.1.a \
 -lfftw3 \
 $(ITK_LIBS)libITKDICOMParser-4.1.a \
 $(ITK_LIBS)libITKgiftiio-4.1.a \
 $(ITK_LIBS)libITKLabelMap-4.1.a \
 $(ITK_LIBS)libITKQuadEdgeMesh-4.1.a \
 $(ITK_LIBS)libITKBiasCorrection-4.1.a \
 $(ITK_LIBS)libITKPolynomials-4.1.a \
 $(ITK_LIBS)libITKBioCell-4.1.a \
 $(ITK_LIBS)libITKFFT-4.1.a \
 -lfftw3 -lfftw3_threads -lfftw3f -lfftw3f_threads \
 $(ITK_LIBS)libITKIOSpatialObjects-4.1.a \
 $(ITK_LIBS)libITKIOXML-4.1.a \
 $(ITK_LIBS)libITKFEM-4.1.a \
 $(ITK_LIBS)libITKOptimizers-4.1.a \
 $(ITK_LIBS)libITKIOBMP-4.1.a \
 $(ITK_LIBS)libITKIOBioRad-4.1.a \
 $(ITK_LIBS)libitkopenjpeg-4.1.a \
 $(ITK_LIBS)libITKEXPAT-4.1.a \
 $(ITK_LIBS)libITKIOGIPL-4.1.a \
 $(ITK_LIBS)libITKIOJPEG-4.1.a \
 $(ITK_LIBS)libITKIOLSM-4.1.a \
 $(ITK_LIBS)libITKIOTIFF-4.1.a \
 $(ITK_LIBS)libitktiff-4.1.a \
 $(ITK_LIBS)libitkjpeg-4.1.a \
 $(ITK_LIBS)libITKIOMeta-4.1.a \
 $(ITK_LIBS)libITKMetaIO-4.1.a \
 $(ITK_LIBS)libITKIONIFTI-4.1.a \
 $(ITK_LIBS)libITKniftiio-4.1.a \
 $(ITK_LIBS)libITKznz-4.1.a \
 $(ITK_LIBS)libITKIONRRD-4.1.a \
 $(ITK_LIBS)libITKNrrdIO-4.1.a \
 $(ITK_LIBS)libITKIOPNG-4.1.a \
 $(ITK_LIBS)libitkpng-4.1.a \
 $(ITK_LIBS)libITKIOSiemens-4.1.a \
 $(ITK_LIBS)libITKIOGE-4.1.a \
 $(ITK_LIBS)libITKIOIPL-4.1.a \
 $(ITK_LIBS)libITKIOStimulate-4.1.a \
 $(ITK_LIBS)libITKIOVTK-4.1.a \
 $(ITK_LIBS)libITKKLMRegionGrowing-4.1.a \
 $(ITK_LIBS)libITKVTK-4.1.a \
 $(ITK_LIBS)libITKWatersheds-4.1.a \
 $(ITK_LIBS)libITKSpatialObjects-4.1.a \
 $(ITK_LIBS)libITKMesh-4.1.a \
 $(ITK_LIBS)libITKPath-4.1.a \
 $(ITK_LIBS)libITKIOTransformHDF5-4.1.a \
 $(ITK_LIBS)libitkhdf5_cpp-4.1.a \
 $(ITK_LIBS)libitkhdf5-4.1.a \
 $(ITK_LIBS)libitkzlib-4.1.a \
 $(ITK_LIBS)libITKIOTransformInsightLegacy-4.1.a \
 $(ITK_LIBS)libITKIOTransformMatlab-4.1.a \
 $(ITK_LIBS)libITKIOTransformBase-4.1.a \
 $(ITK_LIBS)libITKStatistics-4.1.a \
 $(ITK_LIBS)libitkNetlibSlatec-4.1.a \
 $(ITK_LIBS)libITKIOImageBase-4.1.a \
 $(ITK_LIBS)libITKVideoCore-4.1.a \
 $(ITK_LIBS)libITKCommon-4.1.a \
 $(ITK_LIBS)libitksys-4.1.a \
 $(ITK_LIBS)libITKVNLInstantiation-4.1.a \
 $(ITK_LIBS)libitkvnl_algo-4.1.a \
 $(ITK_LIBS)libitkv3p_lsqr-4.1.a \
 $(ITK_LIBS)libitkvnl-4.1.a \
 $(ITK_LIBS)libitkvcl-4.1.a \
 $(ITK_LIBS)libitkv3p_netlib-4.1.a \
 -lpthread \
 -lm \
 -ldl	

FFTW_INCLUDE=$(FFTW_PREFIX)/include/
FFTW_LIB=$(FFTW_PREFIX)/

LINK_LIBRARIES_FFTW=-L$(FFTW_LIB) \
	 -lfftw3 \
	 -lfftw3f \
	 -lfftw3_threads \
	 -lfftw3f_threads

INCLUDE_ITK=-ftemplate-depth-50 -Wall -Wno-deprecated -msse2\
  -IITKIOFactoryRegistration \
	-I$(ITK) \
	-I$(ITK)/vnl \
	-IitkCVLab 


all : plugins/OOFTubularityMeasure_Plugin.jar \
	plugins/TubularGeodesics_Plugin.jar

plugins/%_Plugin.jar : FijiITKInterface/%.class \
			FijiITKInterface/%_Plugin.class \
			%.config \
			fiji/jni/LibraryLoader.class \
			build/$(ARCH)/lib%.$(LIBRARY_EXTENSION)
	mkdir -p plugins
	cp $*.config plugins.config
	jar cvf plugins/$*_Plugin.jar \
		FijiITKInterface/$*.class \
		FijiITKInterface/$*_Plugin.class \
		plugins.config \
		fiji/jni/LibraryLoader.class \
		-C build \
		$(ARCH)/lib$*.$(LIBRARY_EXTENSION)
	rm plugins.config

test :
	java -jar ij.jar -eval 'run("Bridge (174K)"); run("Tubularity Measure Plugin");'

clean :
	rm -fv *.class */*.class *.o *.$(LIBRARY_EXTENSION)
	rm -rf build/$(ARCH)/*
	rm -rf plugins

superclean: clean
	rm -rf ij.jar

build/$(ARCH)/lib%.$(LIBRARY_EXTENSION) : FijiITKInterface/FijiITKInterface_%.h c++/%JNIImplementation.cpp
	mkdir -p build/$(ARCH)/
	g++ -Wall -O3 -DWITH_JAVA -DITK_IO_FACTORY_REGISTER_MANAGER -o $@ -I$(FFTW_INCLUDE)  -I../c++ -fopenmp -lgomp c++/$*JNIImplementation.cpp -fPIC -shared  -I$(JDK_HOME)/include/ -I$(JDK_HOME)/Headers/ -I$(JDK_HOME)/include/$(JAVA_ARCH_NAME)/ -lstdc++ -I./FijiITKInterface/ $(INCLUDE_ITK) $(LINK_LIBRARIES_ITK) $(LINK_LIBRARIES_FFTW)

FijiITKInterface/FijiITKInterface_%.h : FijiITKInterface/%.class
	$(FIJI_LAUNCHER) --javah --class-path=.:$(JDK_HOME)/lib/tools.jar -jni -d FijiITKInterface FijiITKInterface.$*

%.class : %.java
	$(FIJI_LAUNCHER) --javac --class-path=. $<
