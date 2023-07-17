package HA_ORF1P_Prot_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author ORION-CIRB
 */
public class Cell {
    
    public Object3DInt cell;
    public Object3DInt nucleus;
    public Object3DInt innerRing;
    public Object3DInt outerRing;
    public Object3DInt innerNucleus;
    public HashMap<String, Double> params;
    
    public Cell(Object3DInt cell, Object3DInt nucleus) {
        this.cell = cell;
        this.nucleus = nucleus;
        this.params = new HashMap<>();
    }
    
    public void setInnerNucleus(Object3DInt innerNucleus) {
        this.innerNucleus = innerNucleus;
    }
    
    public void setInnerRing(Object3DInt innerRing) {
        this.innerRing = innerRing;
    }
    
    public void setOuterRing(Object3DInt outerRing) {
        this.outerRing = outerRing;
    }
    
    public void setParams(double label, double nucArea, double nucCircV1, double nucCircV2, double nucInt, double innerNucArea, 
            double innerNucInt, double innerRingArea, double innerRingInt, double outerRingArea, double outerRingInt, double cellArea) {
        
        params.put("label", label);
        
        // Nucleus
        params.put("nucArea", nucArea);
        params.put("nucCircV1", nucCircV1);
        params.put("nucCircV2", nucCircV2);
        params.put("nucInt", nucInt);
        
        // Inner nucleus
        params.put("innerNucArea", innerNucArea);
        params.put("innerNucInt", innerNucInt);
        
        // Inner ring
        params.put("innerRingArea", innerRingArea);
        params.put("innerRingInt", innerRingInt);
        
        // Outer ring
        params.put("outerRingArea", outerRingArea);
        params.put("outerRingInt", outerRingInt);
        
        // Cell
        params.put("cellArea", cellArea);
    }
}
