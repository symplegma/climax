/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geom;

import gendomain.Node;

/**
 *
 * @author pchr
 */
public class Triangle extends Shape{

    public Triangle(int id, Node[] thePoints) {
        this.id=id;
        for(int i=0;i<thePoints.length;i++){
            //if(!this.ShapePoints.containsKey(thePoints[i].getID())){
                this.putPoint(thePoints[i]);
                thePoints[i].putShape(this);
            //}
        }
    }

    public Triangle(int id, jbem.Node[] thePoints) {
        this.id=id;
        for(int i=0;i<thePoints.length;i++){
            if(!this.ShapePoints.containsKey(thePoints[i].getID())){
                this.putPoint(thePoints[i]);
                thePoints[i].putShape(this);
            }
        }
    }
    
    public Triangle(int id, Point[] thePoints){
        this.id=id;
        for(int i=0;i<thePoints.length;i++){
            if(!this.ShapePoints.containsKey(thePoints[i].getID())){
                this.putPoint(thePoints[i]);
                thePoints[i].putShape(this);
            }
        }
    }
    
    @Override
    public double ShapeFunction(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=1.0-xsi-eta; break;
            case 2: N=xsi; break;
            case 3: N=eta; break;
            default:N=0.0 ;  break;
        }
        return N;
    }
}
