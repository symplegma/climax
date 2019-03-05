/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geom;

import gendomain.Node;
import java.awt.Color;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
abstract public class Shape {
    protected int id;
    protected Map<Integer,Point> ShapePoints = new TreeMap<Integer,Point>();
    //protected Map<Integer,Integer> theNodesHierarchy = new TreeMap<Integer,Integer>();
    protected int numNodes=0;
    protected Color RestCol,atMotionCol;
    
    protected void putPoint(Point aPoint) {
        ShapePoints.put(++numNodes, aPoint);
        //ShapePoints.put(aPoint.getID(), aPoint);
        aPoint.putShape(this);
    }
    
    public int getID() {
        return id;
    }
    
    public Map<Integer,Point> getPoints(){
        return this.ShapePoints;
    }
    
    public double ShapeFunction(int wsf, double xsi, double eta){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public Node setCentroid(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getMinX(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getMinY(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getMaxX(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getMaxY(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public int getNumNodes(){
        return this.ShapePoints.size();
    }
    
    public Point getNodeHierarchy(int nodeHier) {
        return this.ShapePoints.get(nodeHier);
    }
}
