
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
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import ij.plugin.PlugIn;
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
            
            // Get input directory
            String imgDir = IJ.getDirectory("Select images directory");
            if (imgDir == null) {
                return;
            }
            
            // Find extension of first image in input folder
            String fileExt = tools.findImageType(new File(imgDir));
            // Find all images with corresponding extension in folder
            ArrayList<String> imageFiles = tools.findImages(imgDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("ERROR", "No image found with " + fileExt + " extension in " + imgDir + " folder");
                return;
            }
            
            // Instantiate metadata and reader
            IMetadata meta = MetadataTools.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] chMeta = tools.findChannels(imageFiles.get(0), meta, reader);

            // Generate dialog box
            String[] chOrder = tools.dialog(chMeta);
            if (chOrder == null) {
                return;
            } else if(chOrder[0].equals("None") || chOrder[2].equals("None")) {
                IJ.showMessage("ERROR", "Nuclei or protein channel not defined");
                return;
            }
            
            // Create output directory
            String outDir = imgDir + File.separator + "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write header in results file
            FileWriter fwResults = new FileWriter(outDir + "results.csv", false);
            BufferedWriter results = new BufferedWriter(fwResults);
            results.write("Image name\tROI area (µm2)\tFocused slice\tProtein background\tNucleus ID\tNucleus area (µm2)"
                    + "\tNucleus circularity (v1)\tNucleus circularity (v2)\tNucleus cor. intensity\tNucleus inner area (µm2)"
                    + "\tNucleus inner cor. intensity\tNucleus inner ring area (µm2)\tNucleus inner ring cor. intensity"
                    + "\tNucleus outer ring area (µm2)\tNucleus outer ring cor. intensity\tIs HA-ORF1P?\tCell area (µm2)\n");
            results.flush();
            
            for (String f : imageFiles) {
                reader.setId(f);
                String imgName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + imgName + " ---");
                
                // Load ROIs, if any provided
                tools.print("- Loading ROIs -");
                List<Roi> rois = tools.loadRois(imgDir, imgName, reader);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setCrop(true);
                    
                // For each ROI, open image, crop it and analyze it
                for (Roi roi: rois) {
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                    
                    // Open nuclei channel
                    tools.print("- Opening nuclei channel -");
                    int index = ArrayUtils.indexOf(chMeta, chOrder[0]);
                    ImagePlus stackNuclei = BF.openImagePlus(options)[index];                   
                    ImagePlus imgNuclei = tools.findBestFocus(stackNuclei);
                    int focusedSlice =  Integer.valueOf(imgNuclei.getProp("Label"));
                    tools.closeImage(stackNuclei);
                    
                    // Detect nuclei
                    System.out.println("- Detecting nuclei -");
                    Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNuclei, roi, tools.cellposeNucModel, tools.cellposeNucDiameter, tools.minNucArea, tools.maxNucArea);
                    
                    // If provided, open, crop and analyze HA-ORF1P channel
                    ImagePlus imgHAORF1P = null;
                    ArrayList<Cell> colocPop = new ArrayList<>();
                    if (!chOrder[1].equals("None")) {
                        tools.print("- Opening HA-ORF1p channel -");
                        index = ArrayUtils.indexOf(chMeta, chOrder[1]);
                        ImagePlus stackHAORF1P = BF.openImagePlus(options)[index];
                        imgHAORF1P = tools.findBestFocus(stackHAORF1P);
                        tools.closeImage(stackHAORF1P);

                        // Detect cells
                        System.out.println("- Detecting HA-ORF1p cells -");
                        Objects3DIntPopulation cellPop = tools.cellposeDetection(imgHAORF1P, roi, tools.cellposeCellsModel, tools.cellposeCellsDiameter, tools.minCellArea, tools.maxCellArea);

                        // Colocalize cells with nuclei
                        System.out.println("- Finding HA-ORF1p cells colocalizing with a nucleus -");
                        colocPop = tools.findColocPop(cellPop, nucPop, 0.02);
                    } else {
                        for (Object3DInt nucleus: nucPop.getObjects3DInt())
                                colocPop.add(new Cell(null, nucleus));
                    }
                    
                    // Find nuclei outer and inner ring
                    System.out.println("- Computing nuclei inner and outer rings -");
                    tools.setNucleiRing(colocPop, imgNuclei, tools.outerNucDil, true);
                    // Find nuclei inner ring and inner nucleus
                    tools.setNucleiRing(colocPop, imgNuclei, tools.innerNucDil, false);
                    tools.resetLabels(colocPop);
                    
                    // Open protein channel
                    tools.print("- Opening protein channel -");
                    index = ArrayUtils.indexOf(chMeta, chOrder[2]);
                    ImagePlus stackProt = BF.openImagePlus(options)[index];
                    ImagePlus imgProt = new Duplicator().run​(stackProt, focusedSlice, focusedSlice);
                    
                    // Compute protein background
                    tools.print("- Computing protein channel background noise -");
                    double bgProt = tools.findBackground(stackProt, roi);
                    tools.closeImage(stackProt);

                    // Tag nuclei with parameters
                    tools.print("- Measuring cells parameters -");
                    tools.tagCells(imgProt, colocPop, bgProt);                        
                    
                    // Write results
                    tools.print("- Saving results -");
                    for (Cell cell: colocPop) {
                        double roiArea = tools.computeRoiArea(roi, imgNuclei);
                        HashMap<String, Double> params = cell.params;
                        String isHAORF1P = (cell.cell == null) ? "No" : "Yes";
                        results.write(imgName+"\t"+roiArea+"\t"+focusedSlice+"\t"+bgProt+"\t"+(int)((double)params.get("label"))+
                                "\t"+params.get("nucArea")+"\t"+params.get("nucCircV1")+"\t"+params.get("nucCircV2")+
                                "\t"+params.get("nucInt")+"\t"+params.get("innerNucArea")+"\t"+params.get("innerNucInt")+
                                "\t"+params.get("innerRingArea")+"\t"+params.get("innerRingInt")+"\t"+params.get("outerRingArea")+
                                "\t"+params.get("outerRingInt")+"\t"+isHAORF1P+"\t"+params.get("cellArea")+"\n");
                        results.flush();
                    }
                    
                    // Draw results
                    tools.drawResults(colocPop, imgNuclei, imgHAORF1P, imgName, outDir);
                    
                    tools.closeImage(imgNuclei);
                    if (imgHAORF1P != null) tools.closeImage(imgHAORF1P);
                    tools.closeImage(imgProt);
                }
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(HA_ORF1P_Prot.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        tools.print("--- All done! ---");
    }
}
