
import HA_ORF1P_Prot_Tools.Tools;
import HA_ORF1P_Prot_Tools.Cell;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.Duplicator;
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


/**
 * Detect nuclei and HA-ORF1p cells in 2D
 * Compute their colocalization, distinguish nuclei being HA-ORF1p+ and HA-ORF1p-
 * Measure intensity of protein in different nuclear compartments
 * @author ORION-CIRB
 */
public class HA_ORF1P_Prot implements PlugIn {
    
    Tools tools = new Tools();
    
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            }             
            
            // Get input folder
            String imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with nd extension
            ArrayList<String> imageFiles = tools.findImages(imageDir, "nd");
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            
            // Create output folder
            String outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write header in results file
            String header = "Image name\tROI area (µm2)\tFocused slice\tProtein background\tNucleus ID\tNucleus area (µm2)"
                    + "\tNucleus circularity (v1)\tNucleus circularity (v2)\tNucleus cor. intensity\tNucleus inner area (µm2)"
                    + "\tNucleus inner cor. intensity\tNucleus inner ring area (µm2)\tNucleus inner ring cor. intensity"
                    + "\tNucleus outer ring area (µm2)\tNucleus outer ring cor. intensity\tIs HA-ORF1P?\tCell area (µm2)\n";
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            BufferedWriter results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);

            // Channels dialog
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                // Find ROI file
                RoiManager rm = new RoiManager(false);
                String roiName = imageDir+rootName+".roi";
                if (new File(roiName).exists()) {
                    rm.runCommand("Open", roiName);
                } else {
                    rm.add(new Roi(0, 0, reader.getSizeX(), reader.getSizeY()), 0);
                    System.out.println("No ROI file found, entire image is analyzed");
                }
                    
                // For each ROI, open image, crop it and analyze it
                for (Roi roi : rm.getRoisAsArray()) {
                    ImporterOptions options = new ImporterOptions();
                    options.setId(f);
                    options.setCrop(true);
                    options.setSplitChannels(true);                    
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setQuiet(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                    
                    // Open DAPI channel
                    tools.print("- Analyzing " + chs[0] + " channel -");
                    int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                    ImagePlus stackNuclei = BF.openImagePlus(options)[indexCh];
                    ImagePlus imgNuclei = tools.findBestFocus(stackNuclei);
                    int focusedSlice =  Integer.valueOf(imgNuclei.getProp("Label"));                   
                    tools.closeImg(stackNuclei);
                    
                    // Detect nuclei
                    System.out.println("Detecting nuclei...");
                    Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNuclei, roi, tools.cellposeNucModel, tools.cellposeNucDiameter, tools.minNucArea, tools.maxNucArea);
                    
                    // If provided, open HA-ORF1P channel
                    ImagePlus imgHAORF1P = null;
                    ArrayList<Cell> colocPop = new ArrayList<>();
                    if (!chs[1].equals("None")) {
                        tools.print("- Analyzing " + chs[1] + " channel -");
                        indexCh = ArrayUtils.indexOf(channels, chs[1]);
                        ImagePlus stackHAORF1P = BF.openImagePlus(options)[indexCh];              
                        imgHAORF1P = tools.findBestFocus(stackHAORF1P);
                        tools.closeImg(stackHAORF1P);

                        // Detect cells
                        System.out.println("Detecting HA-ORF1p nuclei...");
                        Objects3DIntPopulation cellPop = tools.cellposeDetection(imgHAORF1P, roi, tools.cellposeCellsModel, tools.cellposeCellsDiameter, tools.minCellArea, tools.maxCellArea);

                        // Colocalize cells with nuclei
                        System.out.println("Finding HA-ORF1p cells colocalizing with a nucleus...");
                        colocPop = tools.findColocPop(cellPop, nucPop, 0.02);
                    } else {
                        for (Object3DInt nucleus: nucPop.getObjects3DInt())
                                colocPop.add(new Cell(null, nucleus));
                    }
                    
                    // Find nuclei outer and inner ring
                    System.out.println("Finding nuclei outer ring...");
                    tools.setNucleiRing(colocPop, imgNuclei, tools.outerNucDil, true);
                    // Find nuclei inner ring and inner nucleus
                    System.out.println("Finding nuclei inner ring and inner...");
                    tools.setNucleiRing(colocPop, imgNuclei, tools.innerNucDil, false);
                    tools.resetLabels(colocPop);
                    
                    // Open protein channel
                    tools.print("- Analyzing " + chs[2] + " channel -");
                    indexCh = ArrayUtils.indexOf(channels, chs[2]);
                    ImagePlus stackProt = BF.openImagePlus(options)[indexCh];
                    ImagePlus imgProt = new Duplicator().run​(stackProt, focusedSlice, focusedSlice);
                    
                    // Compute protein background 
                    double bgProt = tools.findBackground(stackProt, roi);
                    tools.closeImg(stackProt);

                    // Tag nuclei with parameters
                    tools.print("- Measuring cells parameters -");
                    tools.tagCells(imgProt, colocPop, bgProt);                        
                    
                    // Write results
                    tools.print("- Saving results -");
                    for (Cell cell: colocPop) {
                        double roiArea = tools.computeRoiArea(roi, imgNuclei);
                        HashMap<String, Double> params = cell.params;
                        String isHAORF1P = (cell.cell == null) ? "No" : "Yes";
                        results.write(rootName+"\t"+roiArea+"\t"+focusedSlice+"\t"+bgProt+"\t"+(int)((double)params.get("label"))+
                                "\t"+params.get("nucArea")+"\t"+params.get("nucCircV1")+"\t"+params.get("nucCircV2")+
                                "\t"+params.get("nucInt")+"\t"+params.get("innerNucArea")+"\t"+params.get("innerNucInt")+
                                "\t"+params.get("innerRingArea")+"\t"+params.get("innerRingInt")+"\t"+params.get("outerRingArea")+
                                "\t"+params.get("outerRingInt")+"\t"+isHAORF1P+"\t"+params.get("cellArea")+"\n");
                        results.flush();
                    }
                    
                    // Draw results
                    tools.drawResults(colocPop, imgNuclei, imgHAORF1P, rootName, outDirResults);
                    
                    tools.closeImg(imgNuclei);
                    if (imgHAORF1P != null) tools.closeImg(imgHAORF1P);
                    tools.closeImg(imgProt);
                }
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(HA_ORF1P_Prot.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        tools.print("--- All done! ---");
    }
}
