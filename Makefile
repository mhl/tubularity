.PHONY : all clean test

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')

ifeq ($(ARCH),)
ifeq ($(uname_S),Linux)
JAVA_ARCH_NAME=linux
	JDK_INCLUDE_PATH=/usr/lib/jvm/java-6-openjdk/include
	LIBRARY_EXTENSION=so
ifeq ($(uname_M),x86_64)
         ARCH=linux64
else
         ARCH=linux
endif
endif
ifneq (,$(findstring MINGW,$(uname_S)))
	# FIXME: doesn't detect win64 at the moment
        ARCH=win32
endif
ifeq ($(uname_S),Darwin)
        ARCH=macosx
p		JDK_INCLUDE_PATH= /usr/lib/jvm/java-1.6.0-openjdk/include/
		LIBRARY_EXTENSION=dylib
endif
endif

CLASSPATH=.:ij.jar:/home/mark/fiji/plugins/3D_Viewer.jar:/usr/share/java/j3dutils.jar:/usr/share/java/j3dcore.jar:/usr/share/java/vecmath.jar

ITK_HEADER_PREFIX= /home/mark/cvlab-work/ITK
ITK_LIB_LOCATION=  $(ITK_HEADER_PREFIX)/bin

LINK_LIBRARIES_ITK= -lITKIO -lITKAlgorithms -litkopenjpeg -litkpng -litktiff -litkzlib -lITKniftiio -lITKSpatialObject -lITKMetaIO -lITKDICOMParser -lITKEXPAT -lITKznz -lITKNrrdIO -lITKStatistics -litkNetlibSlatec -lITKNumerics -litkv3p_netlib -litkvnl -litkvcl -litkvnl_algo -litkvnl_inst -lITKCommon -litksys -litkgdcm -lpthread -lm -L$(ITK_LIB_LOCATION) -ldl 


INCLUDE_ITK=-ftemplate-depth-50 -Wall -Wno-deprecated -msse2 -I$(ITK_HEADER_PREFIX) -I$(ITK_HEADER_PREFIX)/Code/Algorithms -I$(ITK_HEADER_PREFIX)/Code/BasicFilters -I$(ITK_HEADER_PREFIX)/Code/Common -I$(ITK_HEADER_PREFIX)/Code/Numerics -I$(ITK_HEADER_PREFIX)/Code/IO -I$(ITK_HEADER_PREFIX)/Code/Numerics/FEM -I$(ITK_HEADER_PREFIX)/Code/Numerics/NeuralNetworks -I$(ITK_HEADER_PREFIX)/Code/SpatialObject -I$(ITK_HEADER_PREFIX)/Code/Review -I$(ITK_HEADER_PREFIX)/Utilities/MetaIO -I$(ITK_HEADER_PREFIX)/Utilities/NrrdIO -I$(ITK_HEADER_PREFIX)/bin/Utilities/NrrdIO -I$(ITK_HEADER_PREFIX)/Utilities/DICOMParser -I$(ITK_HEADER_PREFIX)/bin/Utilities/DICOMParser -I$(ITK_HEADER_PREFIX)/bin/Utilities/expat -I$(ITK_HEADER_PREFIX)/Utilities/expat -I$(ITK_HEADER_PREFIX)/Utilities/nifti/niftilib -I$(ITK_HEADER_PREFIX)/Utilities/nifti/znzlib -I$(ITK_HEADER_PREFIX)/Utilities/itkExtHdrs -I$(ITK_HEADER_PREFIX)/bin/Utilities -I$(ITK_HEADER_PREFIX)/Utilities -I$(ITK_HEADER_PREFIX)/Code/Numerics/Statistics -I$(ITK_HEADER_PREFIX)/Utilities/vxl/v3p/netlib -I$(ITK_HEADER_PREFIX)/Utilities/vxl/vcl -I$(ITK_HEADER_PREFIX)/Utilities/vxl/core -I$(ITK_HEADER_PREFIX)/bin/Utilities/vxl/v3p/netlib -I$(ITK_HEADER_PREFIX)/bin/Utilities/vxl/vcl -I$(ITK_HEADER_PREFIX)/bin/Utilities/vxl/core -I$(ITK_HEADER_PREFIX)/bin/Utilities/gdcm -I$(ITK_HEADER_PREFIX)/Utilities/gdcm/src  -I$ /home/mark/cvlab-work/latest-src/src/itkCVLab

all : MultiScaleTubularityMeasure_Plugin.jar

test :
	java -jar ij.jar -eval 'run("Bridge (174K)"); run("Tubularity Measure Plugin");'

clean :
	rm -fv *.class */*.class *.o *.$(LIBRARY_EXTENSION)
	rm -rf build/$(ARCH)/*
	rm -rf plugins

superclean: clean
	rm -rf ij.jar

ij.jar:
	wget http://rsb.info.nih.gov/ij/upgrade/ij.jar

build/$(ARCH)/libMultiScaleTubularityMeasure.$(LIBRARY_EXTENSION) : libMultiScaleTubularityMeasure.$(LIBRARY_EXTENSION)

MultiScaleTubularityMeasure_Plugin.jar : FijiITKInterface/MultiScaleTubularityMeasure.class FijiITKInterface/MultiScaleTubularityMeasure_Plugin.class plugins.config fiji/jni/LibraryLoader.class build/$(ARCH)/libMultiScaleTubularityMeasure.$(LIBRARY_EXTENSION) build/macosx/libTubularityMeasure.dylib 
	mkdir -p plugins
	jar cvf plugins/MultiScaleTubularityMeasure_Plugin.jar FijiITKInterface/MultiScaleTubularityMeasure.class FijiITKInterface/MultiScaleTubularityMeasure_Plugin.class plugins.config fiji/jni/LibraryLoader.class -C build $(ARCH)/libMultiScaleTubularityMeasure.$(LIBRARY_EXTENSION) build/macosx/libTubularityMeasure.dylib

libMultiScaleTubularityMeasure.$(LIBRARY_EXTENSION) : FijiITKInterface/FijiITKInterface_MultiScaleTubularityMeasure.h c++/MultiScaleTubularityJNIImplementation.cpp
	mkdir -p build/$(ARCH)/
	g++ -Wall -O3 -o build/$(ARCH)/$@ -I../c++ c++/MultiScaleTubularityJNIImplementation.cpp -fPIC -shared  -I$(JDK_INCLUDE_PATH) -I$(JDK_INCLUDE_PATH)/$(JAVA_ARCH_NAME)/ -lstdc++ -I./FijiITKInterface/  $(INCLUDE_ITK) $(LINK_LIBRARIES_ITK)

FijiITKInterface/FijiITKInterface_MultiScaleTubularityMeasure.h : FijiITKInterface/MultiScaleTubularityMeasure.class
	javah -classpath $(CLASSPATH) -jni -d FijiITKInterface FijiITKInterface.MultiScaleTubularityMeasure

%.class : %.java ij.jar
	javac -cp $(CLASSPATH)  $<

