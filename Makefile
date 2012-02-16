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
ITK_LIBS=$(ITK)/bin
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

LINK_LIBRARIES_ITK=-L$(ITK_LIBS) \
	-luuid \
	-lITKBasicFilters \
	-lITKNumerics \
	-lITKIO \
	-lITKNrrdIO \
	-litkgdcm \
	-lITKAlgorithms \
	-litkopenjpeg \
	-litkpng \
	-litktiff \
	-litkjpeg8 \
	-litkjpeg16 \
	-litkjpeg12 \
	-litkzlib \
	-lITKniftiio \
	-lITKSpatialObject \
	-lITKMetaIO \
	-lITKDICOMParser \
	-lITKEXPAT \
	-lITKznz \
	-lITKStatistics \
	-litkNetlibSlatec \
	-lITKCommon \
	-litkvnl_inst \
	-litkvnl_algo \
	-litkv3p_netlib \
	-litkvnl \
	-litkvcl \
	-litksys \
	-lpthread \
	-lm \
	-ldl

FFTW_INCLUDE=$(FFTW_PREFIX)/include/
FFTW_LIB=$(FFTW_PREFIX)/lib/

LINK_LIBRARIES_FFTW=-L$(FFTW_LIB) \
	 -lfftw3 \
	 -lfftw3f \
	 -lfftw3_threads \
	 -lfftw3f_threads

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
	g++ -Wall -O3 -o $@ -I$(FFTW_INCLUDE) -I../c++ c++/$*JNIImplementation.cpp itkCVLab/itkFFTWLock.cxx -fPIC -shared  -I$(JDK_HOME)/include/ -I$(JDK_HOME)/Headers/ -I$(JDK_HOME)/include/$(JAVA_ARCH_NAME)/ -lstdc++ -I./FijiITKInterface/ $(INCLUDE_ITK) $(LINK_LIBRARIES_ITK) $(LINK_LIBRARIES_FFTW)

FijiITKInterface/FijiITKInterface_%.h : FijiITKInterface/%.class
	$(FIJI_LAUNCHER) --javah --class-path=.:$(JDK_HOME)/lib/tools.jar -jni -d FijiITKInterface FijiITKInterface.$*

%.class : %.java
	$(FIJI_LAUNCHER) --javac --class-path=. $<
