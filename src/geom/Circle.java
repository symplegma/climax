/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geom;

import gendomain.Node;

/**
 *
 * @author antonia
 */
public class Circle extends Shape{
    private final double R;
    
    public Circle(Point Center, double R){
        this.R=R;
        this.putPoint(Center);
    }
    
    @Override
    public Node setCentroid(){
        return Node.parsePoint(this.ShapePoints.values().iterator().next());
    }
    
    public double getMinX(){
        return (this.ShapePoints.values().iterator().next().X()-R);
    }
    
    public double getMaxX(){
        return (this.ShapePoints.values().iterator().next().X()+R);
    }
    
    public double getMinY(){
        return (this.ShapePoints.values().iterator().next().Y()-R);
    }
    
    public double getMaxY(){
        return (this.ShapePoints.values().iterator().next().Y()+R);
    }
    
    public double getR(){return this.R;}
}
