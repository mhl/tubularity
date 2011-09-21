/* Copyright Johannes Schindelin 2010

   This file is taken from the branch jni-helper of fiji.git:

      http://pacific.mpi-cbg.de/cgi-bin/gitweb.cgi?p=fiji.git;a=shortlog;h=refs/heads/jni-helper

   ... so this version may get rapidly out of date.  Hopefully this
   will be merged into the master branch soon - this copy is just
   for experimentation.
*/

package fiji.jni;

import ij.IJ;
import ij.ImageJ;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import java.util.*;
import java.io.*;
import java.lang.reflect.Field;

import java.net.URL;

/**
 * A simple library loader for JNI.
 *
 * Extend this class, name the new class like your library, and it will be
 * loaded automatically upon construction.
 */

public class LibraryLoader {
	protected static File libraryDirectory;
	protected static boolean libraryLoaded;
	protected String baseURL;

	public static void addDir(String s) throws IOException {
   		 try {
      			  // This enables the java.library.path to be modified at runtime
      			  // From a Sun engineer at http://forums.sun.com/thread.jspa?threadID=707176
      			  //
      			 Field field = ClassLoader.class.getDeclaredField("usr_paths");
       			 field.setAccessible(true);
       			 String[] paths = (String[])field.get(null);
       			 for (int i = 0; i < paths.length; i++) {
        			    if (s.equals(paths[i])) {
        			        return;
        			    }
     			   }
      			  String[] tmp = new String[paths.length+1];
      			  System.arraycopy(paths,0,tmp,0,paths.length);
      			  tmp[paths.length] = s;
      			  field.set(null,tmp);
      			  System.setProperty("java.library.path", System.getProperty("java.library.path") + File.pathSeparator + s);
    			} catch (IllegalAccessException e) {
       			 throw new IOException("Failed to get permissions to set library path");
   			 } catch (NoSuchFieldException e) {
      			  throw new IOException("Failed to get field handle to set library path");
   			 }
	}

	protected LibraryLoader() {

		if (baseURL == null) {
			String classFile = getClass().getName().replace('.', '/') + ".class";
			URL url = getClass().getResource("/" + classFile);
			if (url != null) {
				String string = url.toString();
				if (string.endsWith(classFile))
					baseURL = string.substring(0, string.length() - classFile.length());
			}
		}

		if (!libraryLoaded) {
			//Trick to set java.library.path from the current directory
			String currentDir = System.getProperty("user.dir");
            		String platform = getPlatform();
            
			try{
				System.setProperty("java.library.path", "./plugins/ITK/" + platform + ":" + currentDir + "/plugins/ITK/" + platform);
            		   	String path = System.getProperty("java.library.path");

				//The variable sys_paths is re-initialized if it set to null
				Field fieldSysPath = ClassLoader.class.getDeclaredField( "sys_paths" );
				fieldSysPath.setAccessible( true );
				fieldSysPath.set( null, null );

				//Need to call loadLibrary to re-initilaize the sys_paths
				System.loadLibrary("itkzlib"); System.loadLibrary("ITKNrrdIO");
				System.loadLibrary("itkopenjpeg"); 

				if(platform.equals("linux64")){
					System.loadLibrary("gdcmCommon");System.loadLibrary("gdcmDSED"); System.loadLibrary("gdcmIOD");System.loadLibrary("gdcmDICT"); System.loadLibrary("gdcmjpeg8"); 
					System.loadLibrary("gdcmjpeg12"); System.loadLibrary("gdcmjpeg16"); System.loadLibrary("gdcmcharls");
				}
				
				System.loadLibrary("itkpng"); System.loadLibrary("itkjpeg8"); System.loadLibrary("itkjpeg12"); System.loadLibrary("itkjpeg16"); System.loadLibrary("itktiff"); 
				System.loadLibrary("itkv3p_netlib"); System.loadLibrary("itkvcl"); System.loadLibrary("itkvnl"); System.loadLibrary("itkvnl_algo");
               			System.loadLibrary("itkvnl_inst"); 
				System.loadLibrary("itksys"); System.loadLibrary("ITKCommon");
				System.loadLibrary("ITKMetaIO");System.loadLibrary("ITKSpatialObject");
				System.loadLibrary("ITKDICOMParser"); System.loadLibrary("ITKEXPAT"); 
				System.loadLibrary("ITKznz"); System.loadLibrary("ITKniftiio");

				System.loadLibrary("itkNetlibSlatec"); 

				System.loadLibrary("ITKNumerics");
				System.loadLibrary("ITKStatistics"); 
				
				System.loadLibrary("itkgdcm");System.loadLibrary("ITKIO"); System.loadLibrary("ITKAlgorithms");

		
			} catch (Exception e){
				throw new RuntimeException("Could not extract : " + e);
			}

			loadLibrary();
			libraryLoaded = true;
		}
	}

	protected static String getPlatform() {
		return IJ.isMacOSX() ? "macosx" :
			(IJ.isWindows() ? "win"
				+ (IJ.is64Bit() ? "64" : "32")
			 : "linux" + (IJ.is64Bit() ? "64" : ""));
	}

	protected static String getLibraryName(String name) {
		return (IJ.isWindows() ? "" : "lib") + name
			+ "." + (IJ.isMacOSX() ? "dylib" :
				IJ.isWindows() ? "dll" : "so");
	}

	// return $IMAGEJ_ROOT/lib/$PLATFORM/
	protected static File getLibraryDirectory() {
		if (ImageJ.VERSION.compareTo("1.43d") < 0)
			// Our best guess for ImageJ < 1.43d
			return new File(IJ.getDirectory("plugins")).getParentFile();
		else
			return new File(IJ.getDirectory("imagej"), "lib"
				+ File.separator + getPlatform());
	}

	protected static File getTempLibraryDirectory(String name) {
		try {
			File tmp = File.createTempFile(name, "");
			if (!tmp.delete() || !tmp.mkdirs())
				return null;
			tmp.deleteOnExit();
			return tmp;
		} catch (IOException e) {
			return null;
		}
	}

	protected void loadLibrary() {
		
		String className = getClass().getName();
		loadLibrary(className.substring(className.lastIndexOf('.') + 1));
	}

	protected void loadLibrary(String name) {
	String fileName = getLibraryName(name);
	
        String currentDir = System.getProperty("user.dir");
        String platform = getPlatform();

	File newlibraryDirectory = new File(currentDir + "/plugins/ITK/" + platform);
		
       /*if (libraryDirectory == null){
            libraryDirectory = new File(currentDir + "/plugins/ITK/" + platform);
        }*/

	File newbaseURL = new File(currentDir);

	File file = new File(newlibraryDirectory, fileName);
	if (!file.exists()) {
		// Try to write to the 
		if (baseURL == null)
			throw new RuntimeException("Could not determine .jar");
			try {
				URL url = new URL(newbaseURL + "/plugins/ITK/" + getPlatform() + "/" + fileName);
				try {
					copy(url, file);
				} catch (IOException e) {
					// Try again with temporary directory
					libraryDirectory = getTempLibraryDirectory(name);
					file = new File(newlibraryDirectory, fileName);
					copy(url, file);
				}
			} catch (Exception e) {
				throw new RuntimeException("Could not extract " + fileName + ": " + e);
			}
		}
        
		System.load(file.getAbsolutePath());

	}

	protected static void copy(URL source, File target) throws IOException {
		InputStream in = source.openStream();
		target.deleteOnExit();
		OutputStream out = new FileOutputStream(target);
		byte[] buffer = new byte[1<<16];
		for (;;) {
			int len = in.read(buffer);
			if (len < 0)
				break;
			out.write(buffer, 0, len);
		}
		in.close();
		out.close();
	}
}

