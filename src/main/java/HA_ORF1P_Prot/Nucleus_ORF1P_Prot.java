package HA_ORF1P_Prot;

import HA_ORF1P_Tools.Tools3D;
import HA_ORF1P_Tools.Cell;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.FileWriter;
import java.util.ArrayList;
import loci.common.Region;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/*
 * Find ORF1P cells and their corresponding DAPI nuclei
 * Measure nucleus and cytoplasm intensity in ORF1P channel (488)
 */

/**
 * @author phm
 */
public class Nucleus_ORF1P_Prot implements PlugIn {
    
    Tools3D tools = new Tools3D();
    private final boolean canceled = false;
    private String imageDir = "";
    public  String outDirResults = "";
    public  String rootName = "";
    public BufferedWriter results;
    
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage("Plugin canceled");
                return;
            }
            
                        
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            // Find images with nd extension
            ArrayList<String> imageFile = tools.findImages(imageDir, "nd");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header in results file
            String header = "Image name\tImage background\tNucleus ID\tIs HA-ORF1P\tNucleus area (µm3)\tNucleus compactness\tNucleus circularity\tNucleus elongation"
                    + "\tNucleus cor. intensity\tNucleus inner area (µm2)\tNucleus inner cor. intensity\tNucleus inner ring area (µm2)\tNucleus inner ring cor. intensity"
                    + "\tNucleus outer ring area (µm2)\tNucleus outer ring cor. intensity\tCell area (µm2)\tCell cor. intensty\tLamin branches\tLamin junctions\n";
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFile.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFile.get(0), meta, reader);

            // Channels dialog
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showStatus("Plugin cancelled or Stardist model not found");
                return;
            }
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                // Find ROI file
                String roi_file = imageDir+rootName+".zip";
                RoiManager rm = new RoiManager(false);
                if (new File(roi_file).exists()) {
                    rm.runCommand("Open", roi_file);
                }
                else {
                    roi_file = imageDir+rootName+".roi";
                    if (new File(roi_file).exists()) {
                        rm.runCommand("Open", roi_file);   
                    }
                    else {
                        IJ.showStatus("No ROI file found !");
                        rm.add(new Roi(0, 0, reader.getSizeX()-1, reader.getSizeY()-1), 0);
                    }
                }
                    
                // for each roi open image and crop
                for (Roi roi : rm.getRoisAsArray()) {
                    ImporterOptions options = new ImporterOptions();
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setQuiet(true);
                    options.setCrop(true);
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setCrop(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                    
                    // Open DAPI channel
                    tools.print("- Analyzing " + chs[0] + " channel -");
                    int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                    ImagePlus imgNucleus = BF.openImagePlus(options)[indexCh];
                    ImagePlus imgNucFocus = tools.findBestFocus(imgNucleus);
                    tools.flush_close(imgNucleus);
                    
                    // Find DAPI nuclei
                    System.out.println("Finding " + chs[0] + " nuclei....");
                    Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNucFocus, tools.cellPoseNucDiameter, roi);
                    
                    // Open HA-ORF1P channel
                    ImagePlus imgORF1PFocus = null;
                    ArrayList<Cell> colocPop = new ArrayList<>();
                    if (!chs[1].equals("None")) {
                        tools.print("- Analyzing " + chs[1] + " channel -");
                        indexCh = ArrayUtils.indexOf(channels, chs[1]);
                        ImagePlus imgORF1P = BF.openImagePlus(options)[indexCh];              
                        imgORF1PFocus =tools.findBestFocus(imgORF1P);
                        tools.flush_close(imgORF1P);

                        // Find cells
                        System.out.println("Finding " + chs[1] + " cells....");
                        Objects3DIntPopulation cellPop = tools.cellposeDetection(imgORF1PFocus, tools.cellPoseCellsDiameter, roi);

                        // Colocalization
                        System.out.println("Finding " + chs[1] + " cells colocalizing with a " + chs[0] + " nuclei...");
                        colocPop = tools.findColocPop(cellPop, nucPop, 0.02);
                        tools.resetLabels(colocPop);
                    }
                    
                    // No HA-ORF1P
                    if (colocPop.size() == 0)
                        for (Object3DInt nucObj : nucPop.getObjects3DInt())
                                colocPop.add(new Cell(null, nucObj));
                    
                    
                    // Find nuclei outer ring
                    // No Prot
                    if (!chs[3].equals("None")) {
                        // Find nuclei outer ring
                        System.out.println("Finding outer rings...");
                        tools.setNucleiRing(colocPop, imgNucFocus, tools.outerNucDil, true);
                        // Find nuclei inner ring and inner nucleus
                        System.out.println("Finding inner rings and inner nuclei...");
                        tools.setNucleiRing(colocPop, imgNucFocus, tools.innerNucDil, false);
                    }
                    
                    // Find lamin parameters
                    // Open lamin channel
                    tools.print("- Analyzing " + chs[2] + " channel -");
                    indexCh = ArrayUtils.indexOf(channels, chs[2]);
                    ImagePlus imgLamin = BF.openImagePlus(options)[indexCh];
                    ImagePlus imgLaminFocus =tools.findBestFocus(imgLamin);
                    tools.flush_close(imgLamin);
                    
                    // Open protein channel
                    ImagePlus imgProtFocus = null;
                    double bgProt = 0;
                    // No prot
                    if (!chs[3].equals("None")) {
                        tools.print("- Analyzing " + chs[3] + " channel -");
                        indexCh = ArrayUtils.indexOf(channels, chs[3]);
                        ImagePlus imgProt = BF.openImagePlus(options)[indexCh];
                        imgProtFocus = tools.findBestFocus(imgProt);
                        tools.flush_close(imgProt);
                        // Find background 
                        bgProt = tools.findBackground(imgProtFocus, roi);
                    }

                    // Tag nuclei with parameters
                    tools.print("- Measuring cells parameters -");
                    tools.tagCells(imgProtFocus, imgLaminFocus, colocPop);
                    if (imgProtFocus != null)
                        tools.flush_close(imgProtFocus);
                    
                    // Write results
                    for (Cell cell: colocPop) {
                        HashMap<String, Double> params = cell.params;
                        String isHAORF1P = (cell.cell == null) ? "No" : "Yes";
                        results.write(rootName+"\t"+bgProt+"\t"+(int)((double)params.get("index"))+"\t"+isHAORF1P+"\t"+params.get("nucArea")+"\t"+params.get("nucComp")+"\t"
                                +params.get("nucCirc")+"\t"+params.get("nucEllElong")+"\t"+(params.get("nucInt")-bgProt*params.get("nucArea"))+"\t"
                                +params.get("innerNucArea")+"\t"+(params.get("innerNucInt")-bgProt*params.get("innerNucArea"))+"\t"+params.get("innerRingArea")+"\t"
                                +(params.get("innerRingInt")-bgProt*params.get("innerRingArea"))+"\t"+params.get("outerRingArea")+"\t"
                                +(params.get("outerRingInt")-bgProt*params.get("outerRingArea"))+"\t"+params.get("cellArea")+"\t"
                                        +(params.get("cellInt")-bgProt*params.get("cellArea"))+"\t"+params.get("nucBranches")+"\t"+params.get("nucJunctions")+"\n");
                        results.flush();
                    }
                    
                    // Save image objects
                    tools.print("- Saving results -");
                    tools.drawResults(colocPop, imgNucFocus, imgORF1PFocus, imgLaminFocus, rootName, outDirResults, 40, !chs[3].equals("None"));
                    tools.flush_close(imgNucFocus);
                    tools.flush_close(imgLaminFocus);
                    
                    if (imgORF1PFocus != null)
                        tools.flush_close(imgORF1PFocus);
                }
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(Nucleus_ORF1P_Prot.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        tools.print("--- All done! ---");
    }
}
