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
		if (libraryDirectory == null)
			libraryDirectory = getLibraryDirectory();

		File file = new File(libraryDirectory, fileName);
		if (!file.exists()) {
			// Try to write to the 
			if (baseURL == null)
				throw new RuntimeException("Could not determine .jar");
			try {
				URL url = new URL(baseURL + getPlatform() + "/" + fileName);
				try {
					copy(url, file);
				} catch (IOException e) {
					// Try again with temporary directory
					libraryDirectory = getTempLibraryDirectory(name);
					file = new File(libraryDirectory, fileName);
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

