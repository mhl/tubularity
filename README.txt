// Contact: F. Benmansour, fethallah@gmail.com
            M. Longair


The required libraries with the appropriate options are listed here:
1- FFTW 3.3.1
2- ITK 4.1

For building on Linux or MacOS systems, use a short shell script
that sets up some environment variables and invokes make, like:

  export JDK_HOME=/usr/lib/jvm/java-6-openjdk/
  export FIJI_LAUNCHER=/home/mark/fiji/fiji
  export ITK=/home/mark/cvlab-work/ITK
  export ITK_LIBS=/usr/local/lib/
  make "$@"

The challenges in building the plugins are chiefly that:

 - The libraries that are loaded via JNI should be built
  statically.

 - The ITK code uses some classes from the "review" section of
  ITK, which typically isn't built for packaged versions.

 - ITK should be built to use fftw, so you have to build fftw
  before ITK.


------------------------------------------------------------------------
-------------------              FFTW                -------------------
need to install fftw, so follow the instructions here:

 http://www.fftw.org/install/mac.html
 http://www.fftw.org/install/linux.html

downloaded fftw-3.3.tar.gz, unpacked. configured with:

$ ./configure  --enable-threads  --enable-sse2 
$ make
$ make install
$ make clean

(the --enable-portable-binary option seems to be unknown)

Then configure for float - same configure but with --enable-float:

$ ./configure --enable-threads --enable-sse2 --enable-float 
$ make
$ make install
$ make clean

------------------------------------------------------------------------
-------------------              CMake               -------------------
downloaded CMake from http://www.cmake.org/cmake/resources/software.html
 (select the option to add it to all users’ paths)


------------------------------------------------------------------------
-------------------               ITK                -------------------
Download ITK4.1, www.itk.org


cd ITK

start cmake
select the ITK directory as the source code
select the an other ITK directory for the “where to build the binaries” option
then click “configure”

configuration should go through - a few options will be
highlighted in red.

check “advanced” and turn off “build testing”, “build examples”.
Turn on USE_FFTWD USE_FFTWF and USE_REVIEW.
Turn off BUILD_SHARED_LIBS.
CMAKE_CXX_FLAGS                -fPIC
CMAKE_C_FLAGS                  -fPIC

click configure - it should give errors to say that FFTWD_LIB
and four other paths couldn’t be found.  Set them as, for
example:

 FFTWD_LIB                        /usr/local/lib/libfftw3.a                                                                                                                                                                                                                  
 FFTWD_THREADS_LIB                /usr/local/lib/libfftw3_threads.a                                                                                                                                                                                                          
 FFTWF_LIB                        /usr/local/lib/libfftw3f.a                                                                                                                                                                                                                 
 FFTWF_THREADS_LIB                /usr/local/lib/libfftw3f_threads.a                                                                                                                                                                                                         
 FFTW_INCLUDE_PATH                /usr/local/include        

Turn on ITK_USE_REVIEW

configure again

Click "Generate" then exit

make; make install;


------------------------------------------------------------------------
-------------------             PLUGINS              -------------------
Setting properly the environment variables as explained above,  
then run the script !!
