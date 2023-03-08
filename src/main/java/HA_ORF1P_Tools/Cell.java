package HA_ORF1P_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author hm
 */
public class Cell {
    
    public Object3DInt cell;
    public Object3DInt nucleus;
    public Object3DInt innerRing;
    public Object3DInt outerRing;
    public Object3DInt innerNucleus;
    public double branches;
    public double junctions;
    public HashMap<String, Double> params;
    
    public Cell(Object3DInt cell, Object3DInt nucleus) {
        this.cell = cell;
        this.nucleus = nucleus;
        this.params = new HashMap<>();
    }
    
    public void setInnerRing(Object3DInt innerRing) {
        this.innerRing = innerRing;
    }
    
    public void setOuterRing(Object3DInt outerRing) {
        this.outerRing = outerRing;
    }
    
    public void setInnerNucleus(Object3DInt innerNucleus) {
        this.innerNucleus = innerNucleus;
    }
    
    
    public void setParams(double index, double nucVol, double nucComp, double nucSph, double nucEllElong, double nucEllFlat, double nucInt, 
                double innerNucVol, double innerNucInt, double innerRingVol, double innerRingInt, double outerRingVol, double outerRingInt,
                double nucBranches, double nucJunctions) {
        
        params.put("index", index);
        
        // Nucleus
        params.put("nucVol", nucVol);
        params.put("nucComp", nucComp);
        params.put("nucSph", nucSph);
        params.put("nucEllElong", nucEllElong);
        params.put("nucEllFlat", nucEllFlat);
        params.put("nucInt", nucInt);
        params.put("nucBranches", nucBranches);
        params.put("nucJunctions", nucJunctions);
        
        // Inner nucleus
        params.put("innerNucVol", innerNucVol);
        params.put("innerNucInt", innerNucInt);
        
        // Inner ring
        params.put("innerRingVol", innerRingVol);
        params.put("innerRingInt", innerRingInt);
        
        // Outer ring
        params.put("outerRingVol", outerRingVol);
        params.put("outerRingInt", outerRingInt);
        
    }
}
