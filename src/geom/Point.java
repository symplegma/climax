/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package geom;

import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Point {
    private static int numPoints = 0;
    private static double zero=0.0;
    protected int id;
    protected double[] coordinates;
    private Map<Integer,Shape> ConnectedElements = new TreeMap<Integer,Shape>();
    
    // constructor
    public Point(){
    }
    
    public Point(double x, double y, double z){
        this.coordinates=new double[3];
        coordinates[0]=x;
        coordinates[1]=y;
        coordinates[2]=z;
    }
    
    public Point(double x, double y){
        this.coordinates=new double[3];
        coordinates[0]=x;
        coordinates[1]=y;
        coordinates[2]=zero;
    }
    
    public Point(int id, double x, double y){
        double[] coords = new double[2];
        coords[0]=x; coords[1]=y;
        if(id>0){
            ++numPoints;
            this.id=id;
            this.coordinates = new double[3];
            for(int i=0; i<coords.length; i++){
                this.coordinates[i]=coords[i];
            }
            if(coords.length<3)coordinates[2]=zero;
        }else{
            System.err.println("Can not create a point with id = "+id);
        }
    }
    
    public Point(int id, double[] coords){
        if(id>0){
            ++numPoints;
            this.id=id;
            this.coordinates = new double[3];
            for(int i=0; i<coords.length; i++){
                this.coordinates[i]=coords[i];
            }
            if(coords.length<3)coordinates[2]=zero;
        }else{
            System.err.println("Can not create a point with id = "+id);
        }
    }
    
    public Point(Point... P){
        ++numPoints;
        this.coordinates = new double[3];
        this.coordinates[0]=0.0;
        this.coordinates[1]=0.0;
        this.coordinates[2]=0.0;
        for(int i=0;i<P.length;i++){
            this.coordinates[0]+=P[i].X();
            if(P[i].coordinates.length>1)this.coordinates[1]+=P[i].Y();
            if(P[i].coordinates.length>2)this.coordinates[1]+=P[i].Z();
        }
        this.coordinates[0]/=P.length;
        this.coordinates[1]/=P.length;
        this.coordinates[2]/=P.length;
    }
    
    public void setID(int id){
        this.id=id;
    }
    
    // methods
    public double[] getCoordinates(){
        return this.coordinates;
    }
    
    
    public int getID(){
        return this.id;
    }
    
    public int getNumPoints(){
        return Point.numPoints;
    }
    
    public double getDist(Point anotherPoint){
        double dist=0.;
        int n=Math.min(this.coordinates.length, anotherPoint.getCoordinates().length);
        for(int i=0;i<n;i++){
            dist+=(this.coordinates[i]-anotherPoint.coordinates[i])*(this.coordinates[i]-anotherPoint.coordinates[i]);
        }
        dist=Math.sqrt(dist);
        return dist;
    }
    
    public void print(){
        System.out.print(id+" ");
        for(int i=0;i<this.getCoordinates().length;i++){
            System.out.print(getCoordinates()[i]+" ");
        }
        System.out.println();
    }
    
    public void putShape(Shape aShape) {
        this.ConnectedElements.put(aShape.getID(), aShape);
    }
    
    public double X(){return this.coordinates[0];}
    
    public double Y(){return this.coordinates[1];}
    
    public double Z(){return this.coordinates[2];}
}
