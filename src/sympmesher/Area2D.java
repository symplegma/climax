/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sympmesher;

import geom.Point;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Area2D {
    private Map<Integer,Patch> thePatches = new TreeMap<Integer,Patch>();
    private Map<Integer,Curve3P> theCurves = new TreeMap<Integer,Curve3P>();
    private Map<Integer,Point> thePoints = new TreeMap<Integer,Point>();
    private Map<Integer,Point> theBoundaryPoints = new TreeMap<Integer,Point>();
    private List<quad> thequads= new ArrayList<quad>();
    private static int globalIndex=0;
    
    public Area2D(){};

    public void putPatch(Patch aPatch) {
        thePatches.put(aPatch.getID(), aPatch);
        aPatch.setArea2D(this);
    }

    public void putCurve3P(Curve3P aCurve3P) {
        theCurves.put(aCurve3P.getID(), aCurve3P);
    }

    public void putPoint(Point aPoint) {
        globalIndex++;
        aPoint.setID(globalIndex);
        thePoints.put(aPoint.getID(), aPoint);
    }
    
    public void putBoundaryPoint(Point aPoint) {
        aPoint.setID(globalIndex);
        theBoundaryPoints.put(aPoint.getID(), aPoint);
    }
    
    public Curve3P getCurve3P(int id){
        return theCurves.get(id);
    }
    
    public Patch getPatch(int id){
        return thePatches.get(id);
    }
    
    public void putquads(quad aquad){this.thequads.add(aquad);}
    
    public Map<Integer,Curve3P> getCurve3Ps(){
        return theCurves;
    }
    
    public Map<Integer,Patch> getPatches(){
        return thePatches;
    }
    
    public Map<Integer,Point> getPoints(){
        return thePoints;
    }
    
    public Map<Integer,Point> getBoundaryPoints(){
        return theBoundaryPoints;
    }
    
    public int ExistBoundaryLocation(double x, double y, double z){
        int pid=0;
        for (Map.Entry<Integer,Point> entry : theBoundaryPoints.entrySet()) {
            Point tp = new Point(x,y, z);
            double tol=0.00000001;
            if(tp.getDist(entry.getValue())<=tol){
                pid=entry.getKey();
                break;
            }
        }
        return pid;
    }
    
    public void dothemesh(){
        // some check() operation before start the meshing
        for(Iterator<Patch> it=this.thePatches.values().iterator(); it.hasNext();){
            Patch aPatch = it.next();
            aPatch.meshit();
        }
    }
    
    public void dothemeshtr(){
        dothemesh();
        List<quad> thetriangles= new ArrayList<quad>(); 
        for(quad q: thequads){
            thetriangles.add(new quad(q.getN1(),q.getN3(),q.getN4(),q.getN1()));
            thetriangles.add(new quad(q.getN1(),q.getN2(),q.getN3(),q.getN1()));
        }
        thequads= new ArrayList<quad>();
        for(quad q: thetriangles){
            thequads.add(q);
        }
    }
    
    public List<quad> getquads(){return this.thequads;}
    
    public static void initglobalIndex(){globalIndex=0;}
}
