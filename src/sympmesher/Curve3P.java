/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sympmesher;

import geom.Point;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Curve3P {
    private final int id;
    private final Point startPoint,midPoint,endPoint; // these 3 points are not considered in the map thePoints
    private Map<Integer,Patch> thePatches = new TreeMap<Integer,Patch>();
    private int mesh=0;
    private double ratio=0.0;
    
    public Curve3P(int id, Point start, Point mid, Point end){
        this.id=id;
        this.startPoint=start;
        this.midPoint=mid;
        this.endPoint=end;
    }
    
    public Point getPointStart(){return startPoint;}
    public Point getPointMid(){return midPoint;}
    public Point getPointEnd(){return endPoint;}

    public int getID(){
        return this.id;
    }

    public void putPatch(Patch aPatch, int localID) {
        thePatches.put(localID, aPatch);
        double pratio = 0.0;
        int pmesh=0;
        switch(localID){
            case 1:pratio=aPatch.getratioXsi(); pmesh=aPatch.getmeshXsi(); break;
            case 3:pratio=aPatch.getratioXsi(); pmesh=aPatch.getmeshXsi(); break;
            case 2:pratio=aPatch.getratioEta(); pmesh=aPatch.getmeshEta(); break;
            case 4:pratio=aPatch.getratioEta(); pmesh=aPatch.getmeshEta(); break;
        }
        if(this.ratio==0.0){
            ratio=pratio;
        }else if(Math.abs(ratio-pratio)>Double.MIN_VALUE){
            System.err.println("Curve3P with id="+this.id+" found to have incompatible geometric progression ratio");
            for (Map.Entry<Integer,Patch> entry : thePatches.entrySet()) {
                
                switch(entry.getKey() ){
                    case 1:pratio=entry.getValue().getratioXsi();break;
                    case 3:pratio=entry.getValue().getratioXsi();break;
                    case 2:pratio=entry.getValue().getratioEta();break;
                    case 4:pratio=entry.getValue().getratioEta();break;
                }
                
                System.err.println("Local curve position=" + entry.getKey() + ", Patch id=" + entry.getValue().getID()+", ratio: "+pratio);
            }
        }
        if(this.mesh==0){
            mesh=pmesh;
        }else if(mesh!=pmesh){
            System.err.println("Curve3P with id="+this.id+" found to have incompatible mesh parameter");
            for (Map.Entry<Integer,Patch> entry : thePatches.entrySet()) {
                switch(entry.getKey() ){
                    case 1:pmesh=entry.getValue().getmeshXsi();break;
                    case 3:pmesh=entry.getValue().getmeshXsi();break;
                    case 2:pmesh=entry.getValue().getmeshEta();break;
                    case 4:pmesh=entry.getValue().getmeshEta();break;
                }
                
                System.err.println("Local curve position=" + entry.getKey() + ", Patch id=" + entry.getValue().getID()+", mesh: "+pmesh);
            }
        }
    }
    
    public void printOnPatches(){
        System.out.println(this.id+": "+this.mesh+", "+this.ratio);
        for (Map.Entry<Integer,Patch> entry : thePatches.entrySet()) {
            System.out.println("Local curve position=" + entry.getKey() + " for Patch id=" + entry.getValue().getID());
        }
    }
    
    public void print(){
        System.out.println(this.id+": "+this.mesh+", "+this.ratio);
        System.out.println("Start point "+this.startPoint.getCoordinates()[0]+" "+this.startPoint.getCoordinates()[1]);
        System.out.println("Mid point "+this.midPoint.getCoordinates()[0]+" "+this.midPoint.getCoordinates()[1]);
        System.out.println("End point "+this.endPoint.getCoordinates()[0]+" "+this.endPoint.getCoordinates()[1]);
        
    }
}
