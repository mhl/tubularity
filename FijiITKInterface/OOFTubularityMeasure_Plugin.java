package FijiITKInterface;

import ij.IJ;
import ij.plugin.PlugIn;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.io.FileInfo;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.gui.ImageCanvas;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;


import ij.measure.Calibration;

import ij3d.*;

public class OOFTubularityMeasure_Plugin implements MouseListener, MouseMotionListener, PlugIn {

    OOFTubularityMeasure ti;
	boolean ROIpt_1_2 = true;
	boolean position_checked = false;
    final static boolean MULTICOLORED = true;
    JPanel mainPane;
	Point p, ROI_p1, ROI_p2;
	int Slice1, Slice2;
	ImagePlus imagePlus;

    public OOFTubularityMeasure_Plugin() {
        ti = new OOFTubularityMeasure();
    }

    protected String getSavePath( FileInfo info ) {
	String extension = ".oof.nrrd";
	SaveDialog sd;

	if( info == null ) {

	    sd = new SaveDialog("Save output as...",
				"image",
				extension);

	} else {

	    String fileName = info.fileName;
	    String directory = info.directory;

	    String suggestedSaveFilename;

	    suggestedSaveFilename = fileName;

	    sd = new SaveDialog("Save output as...",
				directory,
				suggestedSaveFilename,
				extension);
	}

	String savePath = null;
	if(sd.getFileName()==null) {
	    IJ.error("You must choose a filename - exiting");
	    return null;
	} else {
	    savePath = sd.getDirectory()+sd.getFileName();
	}

	/*
	File file = new File(savePath);
	if ((file!=null)&&file.exists()) {
	    if (!IJ.showMessageWithCancel(
					  "Save output file...", "The file "+
					  savePath+" already exists.\n"+
					  "Do you want to replace it?"))
		return;
	}
	*/

	return savePath;
    }

    public void run( String ignored ) {

		///////initialisation////////////
        imagePlus = IJ.getImage();
		ImageStack stack = imagePlus.getImageStack();
		ImageCanvas imagecanvas = imagePlus.getCanvas();

        if( imagePlus == null ) {
            IJ.error("No image was open...");
            return;
        }

		//TODO: possible types are : GRAY8 GRAY16 and GRAY32.
        int imageType = imagePlus.getType();
        if( ! (imageType == ImagePlus.GRAY8 || imageType == ImagePlus.COLOR_256) ) {
            IJ.error("Not an 8-bit image");
            return;
        }

        ByteProcessor bp = (ByteProcessor)imagePlus.getProcessor();
        byte [] pixelData = (byte [])bp.getPixels();

		ROI_p1 = new Point(); ROI_p2 = new Point();
		ROI_p1.x = 0; ROI_p1.y = 0; ROI_p2.x = imagePlus.getWidth(); ROI_p2.y = imagePlus.getHeight();
		Slice1 = 1; Slice2 = stack.getSize();

		imagecanvas.addMouseListener(this);
   		imagecanvas.addMouseMotionListener(this);

		FileInfo Info = imagePlus.getOriginalFileInfo();
		Calibration Calib = imagePlus.getCalibration();

		int NSlices = stack.getSize(); int width = imagePlus.getWidth(); int height = imagePlus.getHeight(); 
	
		ByteProcessor bpStack;
        byte [] SlicepixelData = new byte[width*height];
		byte [] StackpixelData = new byte[width*height*NSlices];

		float [] StackPixelDataOut = new float[width*height*NSlices];

		///////////////////////////////////

		for(int i=0;i<NSlices;i++){
			bpStack = (ByteProcessor)stack.getProcessor(i+1);
			SlicepixelData = (byte [])bpStack.getPixels();
			for(int j=0;j<width*height;j++)
				StackpixelData[width*height*i + j] = SlicepixelData[j];
		}

		//IDE 
	 	GenericDialog gd = new GenericDialog("Optimally Oriented Flux Options");

		double pixelWidth = Math.abs(Calib.pixelWidth);

        gd.addNumericField("Number of scales", 1, 0);
	 	gd.addNumericField("Minimum scale", pixelWidth, 6);
 	 	gd.addNumericField("Maximum scale", pixelWidth, 6);
		gd.addCheckbox("Show Gaussian smoothed images:",false);
		gd.addCheckbox("Show filtered images at each scale:",false);
		gd.addCheckbox("Show which scales were used at each point:",false);

		gd.showDialog();	
		
		if( gd.wasCanceled() ) return;

 		int scales = (int)Math.round(gd.getNextNumber());

		if( scales < 1 ) {
		 	IJ.error("The minimum number of scales to try is 1");
			return;
        }

    	double minimumScale = gd.getNextNumber();
        double maximumScale = gd.getNextNumber();

        if( maximumScale < minimumScale ) {
		 IJ.error("The maximum scale cannot be less than the minimum scale");
			return;
		}

        boolean showGaussianImages = gd.getNextBoolean();
        boolean showFilteredImages = gd.getNextBoolean();
        boolean showWhichScales = gd.getNextBoolean();

	String outputFilename = getSavePath( Info );
	System.out.println("writing to outputFilename:"+outputFilename);

        ti.OrientedFlux(StackpixelData, StackPixelDataOut, imageType, width, height, NSlices, Calib.pixelWidth, Calib.pixelHeight, Calib.pixelDepth, minimumScale, maximumScale, scales, outputFilename);

		ImageStack newstack = new ImageStack(width, height);
	
		for(int i=0;i<NSlices;i++){
		
			float[] pix = new float[width*height]; // get your bytes somehow
			float max = -10, min = 10;
			for(int j=0;j<width*height;j++){
				pix[j] = StackPixelDataOut[width*height*i + j];
				if(pix[j] > max ) max = pix[j];
				if(pix[j] < min) min = pix[j];
			}
				
			FloatProcessor proc = new FloatProcessor(width, height, pix, null);
			newstack.addSlice("", proc);
		}	
	
		ImagePlus imp = new ImagePlus("Tubularity", newstack);
		imp.setFileInfo(Info);
		imp.setCalibration(imagePlus.getCalibration());

        imp.show();
    }

	/////////////////////////////////////////

	public void mouseClicked(MouseEvent me) {
  	}

  	public void mouseEntered(MouseEvent me) {
  	}

  	public void mouseExited(MouseEvent me) {
  	}

  	public void mousePressed(MouseEvent me) {
   	 	p = me.getPoint();
		if(ROIpt_1_2){
	 		Slice1 = imagePlus.getCurrentSlice();
			ROI_p1.x = p.x; ROI_p1.y = p.y;
			ROIpt_1_2  = false;
	 		IJ.error("Position 1: " + p.x + " " + p.y + " " + Slice1);
		}else{
			Slice2 = imagePlus.getCurrentSlice();
			ROI_p2.x = p.x; ROI_p2.y = p.y;
			ROIpt_1_2  = true;
			position_checked = true;
			IJ.error("Position 2: " + p.x + " " + p.y + " " + Slice2);
		}
  	}

  	public void mouseReleased(MouseEvent me) {
    	p = null;
  	}

  	public void mouseDragged(MouseEvent me) {
    	p = me.getPoint();
  	}

 	public void mouseMoved(MouseEvent me) {
  	}
	/////////////////////////////////////////
}



/*class MouseMotionEvents extends ImagePlus implements MouseListener,
    MouseMotionListener {
  Point p;

  public MouseMotionEvents() {
    addMouseListener(this);
    addMouseMotionListener(this);
  }

  public void mouseClicked(MouseEvent me) {
  }

  public void mouseEntered(MouseEvent me) {
  }

  public void mouseExited(MouseEvent me) {
  }

  public void mousePressed(MouseEvent me) {
    p = me.getPoint();
    repaint();
  }

  public void mouseReleased(MouseEvent me) {
    p = null;
    repaint();
  }

  public void mouseDragged(MouseEvent me) {
    p = me.getPoint();
    repaint();
  }

  public void mouseMoved(MouseEvent me) {
  }

  public void paint(Graphics g) {
    if (p != null) {
      Dimension d = getSize();
      int xc = d.width / 2;
      int yc = d.height / 2;
      g.drawLine(xc, yc, p.x, p.y);
    }
  }
}*/

