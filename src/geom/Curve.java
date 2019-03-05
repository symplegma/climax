/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geom;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Curve {
    private static int numCurves = 0;
    protected int id;
    protected Map<Integer,Point> thePoints = new TreeMap<Integer,Point>();
    
    // constructor
    public Curve(){
    }
    
    public Curve(int id){
        this.id=id;
        numCurves++;
    }
    
    public Curve(int id, Point p1, Point p2){
        this.id=id;
        numCurves++;
        putPoint(p1);
        putPoint(p2);
    }
    
    public void putPoint(Point aPoint){
        this.thePoints.put(aPoint.getID(), aPoint);
    }
    
    public Point getPoint(int id){
        return this.thePoints.get(id);
    }

    public Map<Integer,Point> getPoints(){
        return this.thePoints;
    }

    public int getID() {
        return this.id;
    }
    
    public boolean contains(geom.Point P){
        boolean itis=false;
        geom.Point p1 = null;
        geom.Point p2 = null;
        int count=0;
        for(Iterator<geom.Point> it=this.getPoints().values().iterator(); it.hasNext();){
            if(count==0){p1 = it.next();}else{p2= it.next();}
            count++;
        }
        if(p1==null || p2==null){
            itis=false;
        }else{
            //if (Math.abs(p1.getDist(P)+P.getDist(p2)-p1.getDist(p2))<=Double.MIN_VALUE ){
            if (Math.abs(p1.getDist(P)+P.getDist(p2)-p1.getDist(p2))<=1.0e-10 ){
                itis = true;
            }
        }
        return itis;
    }
}
