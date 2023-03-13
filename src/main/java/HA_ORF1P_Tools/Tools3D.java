package HA_ORF1P_Tools;

import HA_ORF1P_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import HA_ORF1P_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiEnlarger;
import ij.plugin.filter.ThresholdToSelection;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DIntLabelImage;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureCompactness;
import mcib3d.geom2.measurements.MeasureEllipsoid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.processing.BinaryMorpho;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clijx.plugins.Skeletonize;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.SkeletonResult;



/**
 * @author phm
 */
public class Tools3D {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public Calibration cal = new Calibration();
    public float pixArea = 0;
    public double minNucArea= 50;
    public double maxNucArea = 3000;
    public float innerNucDil = 1;
    public float outerNucDil = 1;
    public double minCellArea= 100;
    public double maxCellArea = 5000;
    
  
    
    // Cellpose
    public int cellPoseCellsDiameter = 100;
    public int cellPoseNucDiameter = 80;
    public String cellPoseModel = "cyto2";
    public String cellPoseEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
      
    private final CLIJ2 clij2 = CLIJ2.getInstance(); 
    private final Find_focused_slices focus = new Find_focused_slices();
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
   
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
 /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelID(0, n) == null) ? channels[n] = Integer.toString(n) : meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelID(0, n) == null) ? channels[n] = Integer.toString(n) : meta.getChannelExcitationWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {     
        
        String[] channelsName = {"DAPI", "HA-ORF1P", "lamin", "Prot"}; 
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
          
        gd.addMessage("Channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        int index = 0;
        for (String chName: channelsName) {
            gd.addChoice(chName + ": ", channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus area (µm2):", minNucArea);
        gd.addNumericField("Max nucleus area (µm2):", maxNucArea);   
        
        gd.addMessage("Nuclei donut", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Nucleus outer ring (µm):", outerNucDil);
        gd.addNumericField("Nucleus inner ring (µm):", innerNucDil);
        
        gd.addMessage("Cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min cell area (µm2): ", minCellArea);
        gd.addNumericField("Max cell area (µm2): ", maxCellArea);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.addNumericField("Z calibration (µm):", cal.pixelDepth);
        gd.showDialog();
        
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if(gd.wasCanceled())
            ch = null;
       
        minNucArea = (float) gd.getNextNumber();
        maxNucArea = (float) gd.getNextNumber();
        outerNucDil = (float) gd.getNextNumber();
        innerNucDil = (float) gd.getNextNumber();
        minCellArea= (float) gd.getNextNumber();
        maxCellArea = (float) gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        pixArea = (float) (cal.pixelWidth*cal.pixelHeight);
        return(ch);
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    /**
     * Look for all 2D cells: 
     * - apply CellPose in 2D slice 
     * 
     * @param img
     * @return 
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, int diameter, Roi roi) throws IOException{
        ImagePlus imgDup = img.resize((int)(img.getWidth()*0.5), (int)(img.getHeight()*0.5), "none");
        ImagePlus imgMed = median_filter(img, 2);
        flush_close(imgDup);
        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellPoseModel, 1, diameter, cellPoseEnvDirPath);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgMed);
        ImagePlus imgOut = cellpose.run();
        flush_close(imgMed);
        ImagePlus imgLabels = imgOut.resize(img.getWidth(), img.getHeight(), "none");
        flush_close(imgOut);
        imgLabels.setCalibration(cal);
        if (roi != null)
            clearOutSide(imgLabels, roi);
        
        // Get cells as a population of objects
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));
        System.out.println(pop.getNbObjects() + " cell detections");
        popFilterSize(pop, minNucArea, maxNucArea);
        System.out.println(pop.getNbObjects() + " cells remaining after size filtering");       
        flush_close(imgLabels);
        return(pop);
    } 
        
    /**
     * Find best focus slice in stack
     */
    public ImagePlus findBestFocus(ImagePlus img) {
        focus.setParams(100, 0, false, false);
        ImagePlus imgFocus = focus.run(img);
        cal.pixelDepth = 1;
        return(imgFocus);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param img
     * @param thMed
     * @return 
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        imgBin.getProcessor().multiply(255);
        IJ.run(imgBin, "8-bit", "");
        return(imgBin);
    }
    
     /**
     * Clij2 skeletonize 2D
     */
    private ImagePlus clij_Skeletonize(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLSkel = clij2.create(imgCL);
        Skeletonize.skeletonize(clij2, imgCL, imgCLSkel);
        clij2.release(imgCL);
        ImagePlus imgSkel = clij2.pull(imgCLSkel);
        clij2.release(imgCLSkel);
        return(imgSkel);
    }
    
    /**
     * Find rois of object (contours)
     * @param obj
     * @return 
     */
    
    private Roi computeRoi(Object3DInt obj) {
        ImageHandler labelImageCrop = new Object3DIntLabelImage(obj).getCroppedLabelImage(255);
        // extract selection
        ImageByte binaryCrop = labelImageCrop.thresholdAboveExclusive(0);
        ByteProcessor mask = new ByteProcessor(binaryCrop.sizeX, binaryCrop.sizeY, (byte[]) binaryCrop.getArray1D(0));
        mask.setThreshold(1, 255, ImageProcessor.NO_LUT_UPDATE);
        ImagePlus maskPlus = new ImagePlus("mask", mask);
        ThresholdToSelection tts = new ThresholdToSelection();
        tts.setup("", maskPlus);
        tts.run(mask);
        Roi roi = maskPlus.getRoi();
        return(roi);
    }
    
    /**
     * Analyze skeleton
     * return branche and junction number
     */
    private double[] AnalyzeSkeleton(ImagePlus img, ImageHandler imgLab, Cell cell) {
        AnalyzeSkeleton_ anaSkel = new AnalyzeSkeleton_();
        anaSkel.setup("",img);
        SkeletonResult skelResult = anaSkel.run(AnalyzeSkeleton_.NONE, true, true, img, true, false);
        ImagePlus imgLabeledSkeletons = new ImagePlus("", anaSkel.getLabeledSkeletons());
        ImageHandler segImage = ImageHandler.wrap(imgLabeledSkeletons);
        segImage.setOffset(imgLab.offsetX - 4, imgLab.offsetY - 4,imgLab.offsetZ);
        segImage.setCalibration(cal);
        cell.setLaminSkel(new Object3DInt(segImage));
        int junctions = skelResult.getListOfJunctionVoxels().size();
        int branches = (junctions == 0) ? 0 : skelResult.getBranches().length;
        return(new double[] {branches, junctions});
    }
    
    /**
     * Find lamin parameters
     * Get number of branches and junctions in lamin skeleton
     */
    public double[] findLaminParams(ImagePlus img, Cell cell) {
        ImageHandler labelImage = new Object3DIntLabelImage(cell.nucleus).getCroppedLabelImage();
        BoundingBox bbox = cell.nucleus.getBoundingBox();
        Roi roi = new Roi(bbox.xmin - 4, bbox.ymin - 4, (bbox.xmax - bbox.xmin) + 4, (bbox.ymax - bbox.ymin) + 4);
        img.setRoi(roi);
        ImagePlus imgCrop = new Duplicator().crop(img);
        imgCrop.deleteRoi();
        ImagePlus imgLOG = LOG_filter(imgCrop, 6);
        flush_close(imgCrop);
        ImagePlus imgBin = threshold(imgLOG, "Li");
        flush_close(imgLOG);
        // clearOutside nucleus
        roi = computeRoi(cell.nucleus);
        // enlarge to 1 µm
        Roi roiEnlarged = RoiEnlarger.enlarge(roi, 4);
        clearOutSide(imgBin, roiEnlarged);
        ImagePlus imgSkel = clij_Skeletonize(imgBin);
        flush_close(imgBin);
        double[] params = AnalyzeSkeleton(imgSkel, labelImage, cell);
        flush_close(imgSkel);
        return(params);
    } 
    
     /**
     * Find coloc between pop1 and pop2
     * set label of colocalized object in Id object
     * @param cellsPop
     * @param nucleiPop
     * @param pourc
     * @return number of pop objects colocalized with pop1
     * @throws java.io.IOException
     */
    public ArrayList<Cell> findColocPop (Objects3DIntPopulation cellsPop, Objects3DIntPopulation nucleiPop, double pourc) throws IOException {
        ArrayList<Cell> colocPop = new ArrayList<Cell>();
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucleiPop, cellsPop);
        AtomicInteger ai = new AtomicInteger(0);
        nucleiPop.getObjects3DInt().forEach(obj1 -> {
            List<PairObjects3DInt> list = coloc.getPairsObject1(obj1.getLabel(), true);
            if (!list.isEmpty()) {
                PairObjects3DInt p = list.get(list.size() - 1);
                Object3DInt obj2 = p.getObject3D2();
                if (p.getPairValue() > obj1.size()*pourc) {
                    colocPop.add(new Cell(obj2, obj1));
                    ai.incrementAndGet();
                }
            }
            else
               colocPop.add(new Cell(null, obj1));
        });
        System.out.println(ai.get()+" cells colocalized with nuclei");
        return(colocPop);
    } 
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        roi.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.updateAndDraw();
    }
    
    
    /*
     * Reset labels of cells in population
     */
    public void resetLabels(ArrayList<Cell> cellPop) {
        float label = 1;
        for (Cell cell: cellPop) {
            if (cell.cell != null)
                cell.cell.setLabel(label);
            cell.nucleus.setLabel(label);
            label++;
        }
    }
    
     /**
     * 2D Laplace of Gaussian filter using CLIJ2
     * @param img
     * @param sizeXYZ
     * @return 
     */ 
    public ImagePlus  LOG_filter(ImagePlus img, double sizeXYZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLDOG = clij2.create(imgCL);
       clij2.gaussianBlur(imgCL, imgCLDOG, sizeXYZ, sizeXYZ);
       clij2.release(imgCL);
       ClearCLBuffer imgCLLOG = clij2.create(imgCL);
       clij2.laplaceBox(imgCLDOG, imgCLLOG);
       clij2.release(imgCLDOG);
       ImagePlus imgLOG = clij2.pull(imgCLLOG);
       clij2.release(imgCLLOG);
       return(imgLOG);
    }  
    
    /**
     * Max filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    private ImagePlus max_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DBox(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       return(clij2.pull(imgCLMax));
    } 
    
    /**
     * Min filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    private ImagePlus min_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMin = clij2.create(imgCL);
       clij2.minimum3DBox(imgCL, imgCLMin, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       return(clij2.pull(imgCLMin));
    } 
    
     /**
     * Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus median_filter(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median2DBox(imgCL, imgCLMed, sizeXY, sizeXY);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
        clij2.release(imgCLMed);
       return(imgMed);
    } 
    
    
    /**
     * Morphological operator
     * @param obj
     * @param op
     * @param rad (pixels)
     * @return 
     */
    
    private Object3DInt getMorphologicalObject3D(Object3DInt obj, ImagePlus img, int op, int rad) {
        int ext = (op == BinaryMorpho.MORPHO_DILATE) ? rad + 1 : 1;
        ImageHandler labelImage = new Object3DIntLabelImage(obj).getCroppedLabelImage(ext, ext, 0, 1, false);
        ImagePlus imgCrop = labelImage.getImagePlus();
        ImagePlus imgSeg = null;
        switch (op) {
            case BinaryMorpho.MORPHO_DILATE :
                imgSeg = max_filter(imgCrop, rad, 0);
                break;
            case BinaryMorpho.MORPHO_ERODE :
                imgSeg = min_filter(imgCrop, rad, 0);
                break;
        }
        ImageHandler segImage2 = ImageHandler.wrap(imgSeg);
        segImage2.setOffset(labelImage);
        segImage2.setCalibration(cal);
        Object3DInt objMorpho = new Object3DInt(segImage2);
        if ((op == BinaryMorpho.MORPHO_DILATE) && (new Object3DComputation(objMorpho).touchBorders(ImageHandler.wrap(img), false)))
            objMorpho = null;
        else
            objMorpho.setLabel(obj.getLabel());
        return objMorpho;
    }
    
    
    /**
     * Set the inner/outer ring of nuclei in cell population
     * @param cellsPop
     * @param img
     * @param dilCoef
     * @param dil
     */
    public void setNucleiRing(ArrayList<Cell> cellsPop, ImagePlus img, float dilCoef, boolean dil) {
        int dilCoefPixels  = (int)Math.ceil(dilCoef/cal.pixelWidth);
        for (Cell cell: cellsPop) {
            if (dil) {
                Object3DInt nucDil =  getMorphologicalObject3D(cell.nucleus, img,  BinaryMorpho.MORPHO_DILATE, dilCoefPixels);
                if (nucDil != null) { 
                    Object3DComputation objComputation = new Object3DComputation​(nucDil);
                    Object3DInt donut = objComputation.getObjectSubtracted(cell.nucleus);
                    donut.setLabel(cell.nucleus.getLabel());
                    donut.setVoxelSizeXY(cal.pixelWidth);
                    donut.setVoxelSizeZ(cal.pixelDepth);
                    cell.setOuterRing(donut);
                }
                else
                    cell.nucleus = null;
            } else {
                Object3DInt nucErod = getMorphologicalObject3D(cell.nucleus, img,  BinaryMorpho.MORPHO_ERODE, dilCoefPixels);
                nucErod.setLabel(cell.nucleus.getLabel());
                cell.setInnerNucleus(nucErod);
                Object3DComputation objComputation = new Object3DComputation​(cell.nucleus);
                Object3DInt donut = objComputation.getObjectSubtracted(nucErod);
                donut.setLabel(cell.nucleus.getLabel());
                donut.setVoxelSizeXY(cal.pixelWidth);
                donut.setVoxelSizeZ(cal.pixelDepth);
                cell.setInnerRing(donut);
            }
        } 
        // Remove nucleus touching border after dilation
        cellsPop.removeIf(p-> (p.nucleus == null));
    }
           
    /**
     * Find background image intensity:
     * Z projection over min intensity + read mean intensity
     * @param img
     */
    public double findBackground(ImagePlus img, Roi roi) {
      ImageProcessor imp = img.getProcessor();
      if (roi != null) {
          roi.setLocation(0, 0);
          imp.setRoi(roi);
      }
      double bg = imp.getStatistics().median;
      System.out.println("Background = " + bg);
      return(bg);
    }
      
    
    /**
     * Tag cell with parameters....
     * @param img
     * @param cellsPop
     * @param dotsPop
     */
    public void tagCells(ImagePlus imgProt, ImagePlus imgLamin, ArrayList<Cell> cellsPop) {
        ImageHandler imh = (imgProt == null) ? null : ImageHandler.wrap(imgProt);
        System.out.println("Measuring nucleus parameters ...");
        for (Cell cell: cellsPop) {
           double nucInt = 0;
           double innerNucArea = 0;
           double innerNucInt = 0;
           double innerRingNucArea = 0;
           double innerRingNucInt = 0;
           double outerRingArea = 0;
           double outerRingInt = 0;
           double nucArea = new MeasureVolume(cell.nucleus).getVolumeUnit();
           double nucComp = new MeasureCompactness(cell.nucleus).getValueMeasurement(MeasureCompactness.COMP_CORRECTED);
           double nucCirc = new MeasureCompactness(cell.nucleus).getValueMeasurement(MeasureCompactness.SPHER_CORRECTED);
           double nucElongation = new MeasureEllipsoid(cell.nucleus).getValueMeasurement(MeasureEllipsoid.ELL_ELONGATION);
           if (imgProt != null) {
                // Get nucleus intensity parameter
                nucInt = new MeasureIntensity(cell.nucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);

                // Get inner nucleus parameters
                innerNucArea = new MeasureVolume(cell.innerNucleus).getVolumeUnit();
                innerNucInt = new MeasureIntensity(cell.innerNucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);

                // Get inner ring parameters
                innerRingNucArea = new MeasureVolume(cell.innerRing).getVolumeUnit();
                innerRingNucInt = new MeasureIntensity(cell.innerRing, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);

                // Get outer ring parameters
                outerRingArea = new MeasureVolume(cell.outerRing).getVolumeUnit();
                outerRingInt = new MeasureIntensity(cell.outerRing, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
           }

            // Get lamin parameters
            double[] laminParams = findLaminParams(imgLamin, cell);
            
            // Add all parameters to cell
            cell.setParams(cell.nucleus.getLabel(), nucArea, nucComp, nucCirc, nucElongation, nucInt, innerNucArea, innerNucInt, innerRingNucArea, innerRingNucInt,
                    outerRingArea, outerRingInt, laminParams[0], laminParams[1]);
        }
        if (imgProt != null)
            imh.closeImagePlus();
    }
        
    /**
     * Label object
     * @param obj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int z = (int)(bb.zmin + 0.5*(bb.zmax - bb.zmin)); //bb.zmin;
        int x = (int)(bb.xmin + 0.5*(bb.xmax - bb.xmin)); //bb.xmin - 1;
        int y = (int)(bb.ymin + 0.5*(bb.ymax - bb.ymin)); //bb.ymin - 1;
        img.setSlice(z); //z+1
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
    
    
    /**
     * Save image objects
     * @param cellsPop
     * @param imgDAPI
     * @param imgORF1P
     * @param imgName 
     * @param outDir
     * @param fontSize 
     */
    public void drawResults(ArrayList<Cell> cellsPop, ImagePlus imgDAPI, ImagePlus imgORF1P, ImagePlus imgLamin, String imgName, String outDir, int fontSize, boolean nucRings) {
        ImageHandler imgObj1 = ImageHandler.wrap(imgDAPI).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.createSameDimensions();
        ImageHandler imgObj3 = imgObj1.createSameDimensions();
        ImageHandler imgObj4 = (nucRings) ? imgObj1.createSameDimensions() : null;
        ImageHandler imgObj5 = (nucRings) ? imgObj1.createSameDimensions() : null;
        ImageHandler imgObj6 = imgObj1.createSameDimensions();
        if (cellsPop.size() > 0) {
            for (Cell cell: cellsPop) {
                cell.nucleus.drawObject(imgObj1, 255);
                if (cell.cell != null)
                    cell.cell.drawObject(imgObj2, 255);
                labelObject(cell.nucleus, imgObj3.getImagePlus(), fontSize);
                if (nucRings) {
                    cell.innerRing.drawObject(imgObj4, 255);
                    cell.outerRing.drawObject(imgObj5, 255);
                }
                cell.laminSkel.drawObject(imgObj6, 255);
            }
        }
        
        System.out.println("Saving cells images ...");
        ImagePlus[] imgColors1 = {imgObj3.getImagePlus(), imgObj2.getImagePlus(), imgObj1.getImagePlus(), imgORF1P};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imgName + "_cells.tif"); 
        flush_close(imgObjects);
        
        if (nucRings) {
            System.out.println("Saving rings cell images ...");
            ImagePlus[] imgColors2 = {imgObj3.getImagePlus(), null, null, imgDAPI, imgObj4.getImagePlus(), null, imgObj5.getImagePlus()};
            imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors2, true);
            ImgObjectsFile = new FileSaver(imgObjects);
            ImgObjectsFile.saveAsTiff(outDir + imgName + "_rings.tif");
            flush_close(imgObjects);
        }
        
        System.out.println("Saving lamin skeleton cell images ...");
        ImagePlus[] imgColors3 = {imgObj6.getImagePlus(), imgObj3.getImagePlus(), null, imgLamin};
        imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors3, true);
        imgObjects.setCalibration(cal);
        ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imgName + "_skeletons.tif");
        flush_close(imgObjects);
        
        imgObj1.closeImagePlus();
        imgObj2.closeImagePlus();
        imgObj3.closeImagePlus();
        if (nucRings) {
            imgObj4.closeImagePlus(); 
            imgObj5.closeImagePlus();
        }
    }
}