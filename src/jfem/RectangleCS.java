/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jfem;

/**
 *
 * @author pchr
 */
public class RectangleCS extends CrossSection{
    private double h,w;
    
    public RectangleCS(int id, double width, double height ){
        this.id=id;
        ++numberOfCrossSection;
        w=width;
        h=height;
	A=w*h;
	As2=5./6.*A;
	As3=5./6.*A;
	J1=h*w*(w*w+h*h)/12.;
	J2=h*w*w*w/12.;
	J3=w*h*h*h/12.;
    }
    
    public double getHeight(){return h;}
    
    public double getWidth(){return w;}
    
    @Override
    public double DA_height(){
        return w;
    }
    
    @Override
    public double DJ3_height(){
        return w*h*h/4.;
    }
    
}
