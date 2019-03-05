/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
abstract public class Material {
    int id;
    protected int MatType;
    
    protected double angle1x=0.;
    
    // constructor
    public Material(){}
    
    // methods
    public int getID(){
        return this.id;
    }
    
    public void setAngle(double angle){this.angle1x=angle;}
    
    public double getAngle(){return this.angle1x;}

}
