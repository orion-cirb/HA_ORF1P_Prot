package HA_ORF1P_Tools;


import ij.*;
import ij.process.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.measure.*;

/** Select focused slices from a Z stack. Based on the autofocus algorithm "Normalized Variance" (Groen et al., 1985; Yeo et
 * al., 1993). However, the images are first treated by a sobel edge filter. This provided a better result for fluorescent bead images.
 * Code modified from the "Select Frames With Best Edges" plugin from Jennifer West (http://www.umanitoba.ca/faculties/science/astronomy/jwest/plugins.html)
 * First version 2009-9-27
 * Second version 2010-11-27
 * Third version 2011-2-15
 * Forth version 2011-3-2
 * Fifth version 2020-2-17 output to result table
 * By TSENG Qingzong; qztseng at gmail.com
 */
 
public class Find_focused_slices implements PlugInFilter, Measurements {

    double percent, vThr;
    boolean consecutive, edge;
    
    
    public void setParams(double percent, double vThr, boolean edge, boolean consecutive) {
        this.percent = percent;
        this.vThr = vThr;
        this.edge = edge;
        this.consecutive = consecutive;
    }
   

    public ImagePlus run(ImagePlus imp) {       
        ImageStack stack = imp.getStack();
        int width = imp.getWidth();
        int height = imp.getHeight();
        String name = imp.getTitle();
        ImageStack stack2 = new ImageStack(width, height, imp.getProcessor().getColorModel());
        int fS = 0;

        int size = stack.getSize();
        if (size == 1){
            IJ.error("Stack required.");
            return(null);
        }

        double vMax = 0;
        double[] varA = new double[size];

        for (int slice = 1; slice <= size; slice++) {
            imp.setSlice(slice);
            ImageProcessor ip = imp.getProcessor();
            varA[slice - 1] = calVar(ip);
            if (varA[slice - 1] > vMax) {
                vMax = varA[slice - 1];
                fS = slice;
            }

        }
        if (vMax < vThr) {
            IJ.error("All slices are below the variance threshold value");
            return(null);
        }
		
	int nn = 0; 
	//go through the slices before the best focus slice
        boolean con = true;
        for (int slice = fS-1; slice >0; slice--) {
            if (varA[slice - 1] / vMax >= percent / 100 && varA[slice - 1] > vThr && con == true) {
                imp.setSlice(slice);
                ImageProcessor ip = imp.getProcessor();
                ip.resetRoi();
                ip = ip.crop();
                String label = stack.getSliceLabel(slice);
                if (label == null) {
                    label = "Z";
                }
                stack2.addSlice(label + "_" + slice, ip,0);
            }else{
            	if(consecutive)	con = false;	
            }
        }
	//go through the slices after the best focus slice 
	con = true;             
        for (int slice = fS; slice <= size; slice++) {
            if (varA[slice - 1] / vMax >= percent / 100 && varA[slice - 1] > vThr && con == true) {
                imp.setSlice(slice);
                ImageProcessor ip = imp.getProcessor();
                ip.resetRoi();
                ip = ip.crop();
                String label = stack.getSliceLabel(slice);
                if (label == null) {
                    label = "Z";
                }
                stack2.addSlice(label + "_" + slice, ip, nn);
                nn++;
            } else {
            	if(consecutive)	con = false;	
            }
        }
		
        ImagePlus focusstack = imp.createImagePlus();
        focusstack.setStack("Focused slices of " + name + "_" + percent + "%", stack2);
        focusstack.setCalibration(imp.getCalibration());
        if (focusstack.getStackSize() == 1) {
            focusstack.setProperty("Label", "Z_" + fS);
        }
        return(focusstack);
    }

    double calVar(ImageProcessor ip) {
        double variance = 0;
        int W = ip.getWidth();
        int H = ip.getHeight();

        Rectangle r = ip.getRoi();
        if (r == null) {
            r.x = 0;
            r.y = 0;
            r.height = H;
            r.width = W;
        }
        ImageProcessor edged = ip.duplicate();
        if (edge) edged.findEdges();
        double mean = ImageStatistics.getStatistics(edged, MEAN, null).mean;
        double a = 0;
        for (int y = r.y; y < (r.y + r.height); y++) {
            for (int x = r.x; x < (r.x + r.width); x++) {
                a += Math.pow(edged.getPixel(x, y) - mean, 2);
            }
        }
        variance = (1 / (W * H * mean)) * a;
        return variance;
    }
    
    
    public void run(ImageProcessor arg0) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    
    public int setup(String arg0, ImagePlus arg1) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
