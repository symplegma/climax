/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geom;

import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Vertex {
    private static int numVerteces = 0;
    protected int id;
    protected Map<Integer,Curve> theCurves = new TreeMap<Integer,Curve>();
    
    // constructor
    public Vertex(){
    }
    
    public Vertex(int id){
        this.id=id;
        numVerteces++;
    }
    
    public void putCurve(Curve aLine){
        this.theCurves.put(aLine.getID(), aLine);
    }
    
    public Curve getCurve(int id){
        return this.theCurves.get(id);
    }

    public Map<Integer,Curve> getCurves(){
        return this.theCurves;
    }

    public int getID() {
        return this.id;
    }
}
