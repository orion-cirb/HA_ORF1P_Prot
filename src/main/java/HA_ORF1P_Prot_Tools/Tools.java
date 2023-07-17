package HA_ORF1P_Prot_Tools;

import HA_ORF1P_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import HA_ORF1P_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
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
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Object3DIntLabelImage;
import mcib3d.geom2.Object3DPlane;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureCompactness;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.processing.BinaryMorpho;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;


/**
 * @author ORION-CIRB
 */
public class Tools {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/HA_ORF1P_Prot.git";
    
    private final CLIJ2 clij2 = CLIJ2.getInstance(); 
    private final Find_focused_slices focus = new Find_focused_slices();
    
    public String[] channelsName = {"DAPI (mandatory)", "HA-ORF1P (optional)", "Protein (mandatory)"}; 
    public Calibration cal = new Calibration();
    public float pixArea;
   
    // Cellpose
    public String cellposeEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    public String cellposeNucModel = "cyto2";
    public int cellposeNucDiameter = 120;
    public double minNucArea = 20;
    public double maxNucArea = 200;
    public String cellposeCellsModel = "cyto2";
    public int cellposeCellsDiameter = 200;
    public double minCellArea = 100;
    public double maxCellArea = 400;
    
    // Nuclei rings
    public float outerNucDil = 1;
    public float innerNucDil = 1;
    
      

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
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
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
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
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
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 60, 0);
        gd.addImage(icon);
          
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName: channelsName) {
            gd.addChoice(chName + ": ", channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min area (µm2):", minNucArea);
        gd.addNumericField("Max area (µm2):", maxNucArea);   
        
        gd.addMessage("Nuclei rings", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Outer ring (µm):", outerNucDil);
        gd.addNumericField("Inner ring (µm):", innerNucDil);
        
        gd.addMessage("HA-ORF1p cells detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min area (µm2): ", minCellArea);
        gd.addNumericField("Max area (µm2): ", maxCellArea);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
       
        minNucArea = (float) gd.getNextNumber();
        maxNucArea = (float) gd.getNextNumber();
        outerNucDil = (float) gd.getNextNumber();
        innerNucDil = (float) gd.getNextNumber();
        minCellArea= (float) gd.getNextNumber();
        maxCellArea = (float) gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixArea = (float) (cal.pixelWidth*cal.pixelHeight);
        
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImg(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Find best focuseds slice in stack
     */
    public ImagePlus findBestFocus(ImagePlus img) {
        focus.setParams(100, 0, false, false);
        ImagePlus imgFocus = focus.run(img);
        return(imgFocus);
    }
    
    
    /**
     * Apply 2D CellPose
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, Roi roi, String model, int diameter, double minArea, double maxArea) throws IOException{
        // Duplicate image
        ImagePlus imgDup = img.duplicate();
        
        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(model, 1, diameter, cellposeEnvDirPath);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgDup);
        ImagePlus imgOut = cellpose.run();
        imgOut.setCalibration(cal);
        clearOutside(imgOut, roi);
        
        // Get cells as a population of objects
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        System.out.println(pop.getNbObjects() + " Cellpose detections");
        popFilterSize(pop, minArea, maxArea);
        System.out.println(pop.getNbObjects() + " detections remaining after size filtering");
        
        closeImg(imgDup);
        closeImg(imgOut);
        return(pop);
    } 
    
    
    /**
     * Median filter using CLIJ2
     */ 
    public ImagePlus medianFilter(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median2DBox(imgCL, imgCLMed, sizeXY, sizeXY);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCL);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
    
    /**
     * Clear outside ROI
     */
    public void clearOutside(ImagePlus img, Roi roi) {
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
    
    
    /**
     * Remove object with size < min and size > max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    /**
     * Compute colocalization between cells and nuclei
     * @throws java.io.IOException
     */
    public ArrayList<Cell> findColocPop(Objects3DIntPopulation cellsPop, Objects3DIntPopulation nucleiPop, double perc) throws IOException {
        AtomicInteger ai = new AtomicInteger(0);
        ArrayList<Cell> colocPop = new ArrayList<Cell>();
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucleiPop, cellsPop);
        nucleiPop.getObjects3DInt().forEach(nucleus -> {
            List<PairObjects3DInt> list = coloc.getPairsObject1(nucleus.getLabel(), true);
            if (!list.isEmpty()) {
                PairObjects3DInt p = list.get(list.size() - 1);
                Object3DInt cell = p.getObject3D2();
                if (p.getPairValue() > nucleus.size()*perc) {
                    colocPop.add(new Cell(cell, nucleus));
                    ai.incrementAndGet();
                }
            } else {
               colocPop.add(new Cell(null, nucleus));
            }
        });
        System.out.println(ai.get()+" cells colocalized with a nucleus");
        return(colocPop);
    } 
    

    /*
     * Reset labels of cells in population
     */
    public void resetLabels(ArrayList<Cell> cellPop) {
        float label = 1;
        for (Cell cell: cellPop) {
            cell.nucleus.setLabel(label);
            cell.innerNucleus.setLabel(label);
            cell.innerRing.setLabel(label);
            cell.outerRing.setLabel(label);
            if (cell.cell != null) cell.cell.setLabel(label);
            label++;
        }
    }

    
    /**
     * Compute the inner/outer ring of nuclei in cell population
     */
    public void setNucleiRing(ArrayList<Cell> cellsPop, ImagePlus img, float dilCoef, boolean dil) {
        int dilCoefPix  = (int) Math.ceil(dilCoef / cal.pixelWidth);
        for (Cell cell: cellsPop) {
            if (dil) {
                Object3DInt nucDil =  getMorphologicalObject3D(cell.nucleus, img, BinaryMorpho.MORPHO_DILATE, dilCoefPix);
                if (nucDil != null) { 
                    Object3DComputation objComputation = new Object3DComputation​(nucDil);
                    Object3DInt ring = objComputation.getObjectSubtracted(cell.nucleus);
                    ring.setVoxelSizeXY(cal.pixelWidth);
                    ring.setVoxelSizeZ(cal.pixelDepth);
                    cell.setOuterRing(ring);
                } else {
                    cell.nucleus = null;
                }
            } else {
                Object3DInt nucErod = getMorphologicalObject3D(cell.nucleus, img,  BinaryMorpho.MORPHO_ERODE, dilCoefPix);
                cell.setInnerNucleus(nucErod);
                
                Object3DComputation objComputation = new Object3DComputation​(cell.nucleus);
                Object3DInt ring = objComputation.getObjectSubtracted(nucErod);
                ring.setVoxelSizeXY(cal.pixelWidth);
                ring.setVoxelSizeZ(cal.pixelDepth);
                cell.setInnerRing(ring);
            }
        } 
        // Remove nuclei touching border after dilation
        cellsPop.removeIf(p-> (p.nucleus == null));
    }
    
    
    /**
     * Morphological operator
     */
    private Object3DInt getMorphologicalObject3D(Object3DInt obj, ImagePlus img, int op, int rad) {
        int ext = (op == BinaryMorpho.MORPHO_DILATE) ? rad + 1 : 1;
        ImageHandler labelImage = new Object3DIntLabelImage(obj).getCroppedLabelImage(ext, ext, 0, 1, false);
        ImagePlus imgCrop = labelImage.getImagePlus();
        ImagePlus imgSeg = null;
        switch (op) {
            case BinaryMorpho.MORPHO_DILATE :
                imgSeg = maxFilter(imgCrop, rad, 0);
                break;
            case BinaryMorpho.MORPHO_ERODE :
                imgSeg = minFilter(imgCrop, rad, 0);
                break;
        }
        ImageHandler imhSeg = ImageHandler.wrap(imgSeg);
        imhSeg.setOffset(labelImage);
        imhSeg.setCalibration(cal);
        
        Object3DInt objMorpho = new Object3DInt(imhSeg);
        if ((op == BinaryMorpho.MORPHO_DILATE) && (new Object3DComputation(objMorpho).touchBorders(ImageHandler.wrap(img), false)))
            objMorpho = null;
        
        return objMorpho;
    }
    
    
    /**
     * Max filter using CLIJ2
     */ 
    private ImagePlus maxFilter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DBox(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);
       ImagePlus imgMax = clij2.pull(imgCLMax);
       clij2.release(imgCL);
       clij2.release(imgCLMax);
       return(imgMax);
    } 
    
    
    /**
     * Min filter using CLIJ2
     */ 
    private ImagePlus minFilter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMin = clij2.create(imgCL);
       clij2.minimum3DBox(imgCL, imgCLMin, sizeXY, sizeXY, sizeZ);
       ImagePlus imgMin = clij2.pull(imgCLMin);
       clij2.release(imgCL);
       clij2.release(imgCLMin);
       return(imgMin);
    } 
    

    /**
     * Find background image intensity:
     * Z projection over min intensity + read mean intensity
     * @param img
     */
    public double findBackground(ImagePlus img, Roi roi) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      roi.setLocation(0, 0);
      imp.setRoi(roi);
      
      double bg = imp.getStatistics().median;
      System.out.println("Background (median intensity of the min projection) = " + bg);
      closeImg(imgProj);
      return(bg);
    }
    
    
    /**
     * Do Z projection
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
          
    /**
     * Tag cells with parameters
     */
    public void tagCells(ImagePlus imgProt, ArrayList<Cell> cellsPop, double bg) {
        ImageHandler imh = ImageHandler.wrap(imgProt);
        
        for (Cell cell: cellsPop) {
            // Get nucleus parameters
            double nucArea = new MeasureVolume(cell.nucleus).getVolumeUnit();
            double nucCircV1 = new MeasureCompactness(cell.nucleus).getValueMeasurement(MeasureCompactness.SPHER_CORRECTED);
            double nucCircV2 = computeNucleusCircularity(cell.nucleus, imgProt);
            double nucInt = new MeasureIntensity(cell.nucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM)-bg*nucArea/pixArea;

            // Get inner nucleus parameters
            double innerNucArea = new MeasureVolume(cell.innerNucleus).getVolumeUnit();
            double innerNucInt = new MeasureIntensity(cell.innerNucleus, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM)-bg*innerNucArea/pixArea;

            // Get inner ring parameters
            double innerRingArea = new MeasureVolume(cell.innerRing).getVolumeUnit();
            double innerRingInt = new MeasureIntensity(cell.innerRing, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM)-bg*innerRingArea/pixArea;

            // Get outer ring parameters
            double outerRingArea = new MeasureVolume(cell.outerRing).getVolumeUnit();
            double outerRingInt = new MeasureIntensity(cell.outerRing, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM)-bg*outerRingArea/pixArea;

            // Cell
            double cellArea = (cell.cell != null) ? new MeasureVolume(cell.cell).getVolumeUnit() : 0;
                    
            // Save all parameters
            cell.setParams(cell.nucleus.getLabel(), nucArea, nucCircV1, nucCircV2, nucInt, innerNucArea, innerNucInt, 
                            innerRingArea, innerRingInt, outerRingArea, outerRingInt, cellArea);
        }
        imh.closeImagePlus();
    }
       
    
    public double computeNucleusCircularity(Object3DInt nucleus, ImagePlus img) {
        Object3DPlane plane = nucleus.getObject3DPlanes().get(0);
        ImageHandler planeImh = ImageHandler.wrap(IJ.createImage("", "8-bit black", img.getWidth(), img.getHeight(), img.getNSlices()));
        plane.drawObject(planeImh, 255);
        ImagePlus planeImg = planeImh.getImagePlus();
        planeImg.setZ(plane.getZPlane()+1);
        IJ.setAutoThreshold(planeImg, "Default dark no-reset");
        
        ResultsTable rt = new ResultsTable();
        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.CLEAR_WORKSHEET+ParticleAnalyzer.LIMIT, Measurements.SHAPE_DESCRIPTORS, rt, 0, Double.MAX_VALUE);
        pa.analyze(planeImg​);
        double circ = rt.getValue("Circ.", 0);
        
        closeImg(planeImg);
        return(circ);
    }
        
    
    /**
     * Compute ROI area
     */
    public double computeRoiArea(Roi roi, ImagePlus img) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        img.setRoi(poly);
        img.setCalibration(cal);
        
        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(img, Analyzer.AREA, rt);
        analyzer.measure();
        return(rt.getValue("Area", 0));
    } 
    
    
    /**
     * Draw results
     */
    public void drawResults(ArrayList<Cell> cellsPop, ImagePlus imgDAPI, ImagePlus imgORF1P, String imgName, String outDir) {
        ImageHandler imgObj1 = ImageHandler.wrap(imgDAPI).createSameDimensions();
        ImageHandler imgObj2 = imgObj1.createSameDimensions();
        ImageHandler imgObj3 = imgObj1.createSameDimensions();
        ImageHandler imgObj4 = imgObj1.createSameDimensions();
        
        if (cellsPop.size() > 0) {
            for (Cell cell: cellsPop) {
                int label = (int)((double)cell.params.get("label"));
                cell.nucleus.drawObject(imgObj1, label);
                if (cell.cell != null) cell.cell.drawObject(imgObj2, label);
                cell.innerRing.drawObject(imgObj3, label);
                cell.outerRing.drawObject(imgObj4, label);
            }
        }
        
        ImagePlus[] imgColors = {null, imgObj2.getImagePlus(), imgObj1.getImagePlus(), imgORF1P};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imgName + "_cells.tif"); 
        closeImg(imgObjects);
        
        ImagePlus[] imgColors1 = {null, null, null, imgDAPI, imgObj3.getImagePlus(), null, imgObj4.getImagePlus()};
        imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
        ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(outDir + imgName + "_rings.tif");
        closeImg(imgObjects);
        
        imgObj1.closeImagePlus();
        imgObj2.closeImagePlus();
        imgObj3.closeImagePlus();
        imgObj4.closeImagePlus();
    }
}