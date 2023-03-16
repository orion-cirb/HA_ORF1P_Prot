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
    public Object3DInt laminSkel;
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
    
    public void setLaminSkel(Object3DInt laminSkel) {
        this.laminSkel = laminSkel;
    }
    
    
    public void setParams(double index, double nucArea, double nucComp, double nucCirc, double nucEllElong, double nucInt, 
                double innerNucArea, double innerNucInt, double innerRingArea, double innerRingInt, double outerRingArea, double outerRingInt,
                double cellArea, double cellInt, double nucBranches, double nucJunctions) {
        
        params.put("index", index);
        
        // Nucleus
        params.put("nucArea", nucArea);
        params.put("nucComp", nucComp);
        params.put("nucCirc", nucCirc);
        params.put("nucEllElong", nucEllElong);
        params.put("nucInt", nucInt);
        params.put("nucBranches", nucBranches);
        params.put("nucJunctions", nucJunctions);
        
        // Cell
        params.put("cellArea", cellArea);
        params.put("cellInt", cellInt);
        
        // Inner nucleus
        params.put("innerNucArea", innerNucArea);
        params.put("innerNucInt", innerNucInt);
        
        // Inner ring
        params.put("innerRingArea", innerRingArea);
        params.put("innerRingInt", innerRingInt);
        
        // Outer ring
        params.put("outerRingArea", outerRingArea);
        params.put("outerRingInt", outerRingInt);
        
    }
}
