.PHONY : all clean test

# Ensure the JDK_INCLUDE_PATH is defined:
#
# Examples:
# On Linux it may well be:
#   /usr/lib/jvm/java-6-openjdk/include/
# On Mac OS:
#   /usr/lib/jvm/java-1.6.0-openjdk/include/
# On Windows:
#   [not sure yet]

ifndef JDK_INCLUDE_PATH
$(error The environment variable JDK_INCLUDE_PATH must be defined.)
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
ITK_LIBS=$(ITK)/bin
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

LINK_LIBRARIES_ITK=-L$(ITK_LIBS) \
	-lITKIO \
	-lITKAlgorithms \
	-litkopenjpeg \
	-litkpng \
	-litktiff \
	-litkzlib \
	-lITKniftiio \
	-lITKSpatialObject \
	-lITKMetaIO \
	-lITKDICOMParser \
	-lITKEXPAT \
	-lITKznz \
	-lITKNrrdIO \
	-lITKStatistics \
	-litkNetlibSlatec \
	-lITKNumerics \
	-litkv3p_netlib \
	-litkvnl \
	-litkvcl \
	-litkvnl_algo \
	-litkvnl_inst \
	-lITKCommon \
	-litksys \
	-lpthread \
	-lm \
	-ldlx

INCLUDE_ITK=-ftemplate-depth-50 -Wall -Wno-deprecated -msse2 -I$(ITK) \
	-I$(ITK)/Code/Algorithms \
	-I$(ITK)/Code/BasicFilters \
	-I$(ITK)/Code/Common \
	-I$(ITK)/Code/Numerics \
	-I$(ITK)/Code/IO \
	-I$(ITK)/Code/Numerics/FEM \
	-I$(ITK)/Code/Numerics/NeuralNetworks \
	-I$(ITK)/Code/SpatialObject \
	-I$(ITK)/Code/Review \
	-I$(ITK)/Utilities/MetaIO \
	-I$(ITK)/Utilities/NrrdIO \
	-I$(ITK)/bin/Utilities/NrrdIO \
	-I$(ITK)/Utilities/DICOMParser \
	-I$(ITK)/bin/Utilities/DICOMParser \
	-I$(ITK)/bin/Utilities/expat \
	-I$(ITK)/Utilities/expat \
	-I$(ITK)/Utilities/nifti/niftilib \
	-I$(ITK)/Utilities/nifti/znzlib \
	-I$(ITK)/Utilities/itkExtHdrs \
	-I$(ITK)/bin/Utilities \
	-I$(ITK)/Utilities \
	-I$(ITK)/Code/Numerics/Statistics \
	-I$(ITK)/Utilities/vxl/v3p/netlib \
	-I$(ITK)/Utilities/vxl/vcl \
	-I$(ITK)/Utilities/vxl/core \
	-I$(ITK)/bin/Utilities/vxl/v3p/netlib \
	-I$(ITK)/bin/Utilities/vxl/vcl \
	-I$(ITK)/bin/Utilities/vxl/core \
	-I$(ITK)/bin/Utilities/gdcm \
	-I$(ITK)/Utilities/gdcm/src \
	-IitkCVLab

all : plugins/MultiScaleTubularityMeasure_Plugin.jar

test :
	java -jar ij.jar -eval 'run("Bridge (174K)"); run("MultiScaleTubularity Measure Plugin");'

clean :
	rm -fv *.class */*.class *.o *.$(LIBRARY_EXTENSION)
	rm -rf build/$(ARCH)/*
	rm -rf plugins

superclean: clean
	rm -rf ij.jar

plugins/TubularityMeasure_Plugin.jar : FijiITKInterface/TubularityMeasure.class FijiITKInterface/TubularityMeasure_Plugin.class plugins.config fiji/jni/LibraryLoader.class build/$(ARCH)/libTubularityMeasure.$(LIBRARY_EXTENSION)
	mkdir -p plugins
	jar cvf plugins/TubularityMeasure_Plugin.jar FijiITKInterface/TubularityMeasure.class FijiITKInterface/TubularityMeasure_Plugin.class plugins.config fiji/jni/LibraryLoader.class -C build $(ARCH)/libTubularityMeasure.$(LIBRARY_EXTENSION)

build/$(ARCH)/libTubularityMeasure.$(LIBRARY_EXTENSION) : FijiITKInterface/FijiITKInterface_TubularityMeasure.h c++/TubularityJNIImplementation.cpp
	mkdir -p build/$(ARCH)/
	g++ -Wall -O3 -o $@ -I../c++ c++/TubularityJNIImplementation.cpp -fPIC -shared  -I$(JDK_INCLUDE_PATH) -I$(JDK_INCLUDE_PATH)/$(JAVA_ARCH_NAME)/ -lstdc++ -I./FijiITKInterface/  $(INCLUDE_ITK) $(LINK_LIBRARIES_ITK)

FijiITKInterface/FijiITKInterface_TubularityMeasure.h : FijiITKInterface/TubularityMeasure.class
	$(FIJI_LAUNCHER) --javah --class-path=. -jni -d FijiITKInterface FijiITKInterface.TubularityMeasure

%.class : %.java ij.jar
	$(FIJI_LAUNCHER) --javac --class-path=. $<
