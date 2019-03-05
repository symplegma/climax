/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import gendomain.Node;
import geom.Shape;
import mater.Material;

/**
 *
 * @author antonia
 */
public class Particle {
    private int id;
    private Node theCentroid;
    private Shape theShape;
    private double m,Jm;
    private Material theMaterial;
    private double Fx_o,Fy_o,Mz_o;
    private double Fx_n,Fy_n,Mz_n;
    
    public Particle(int id){
        this.id=id;
    }
    
    public void setShape(Shape theShape){
        this.theShape=theShape;
        
        theCentroid = theShape.setCentroid();
    }
    
    public Shape getShape(){
        return this.theShape;
    }

    public int getID() {
        return this.id;
    }
    
    public Node getCentroid(){
        return theCentroid;
    }
    
    public double getMinX(){
        return this.theShape.getMinX();
    }
    
    public double getMaxX(){
        return this.theShape.getMaxX();
    }
    
    public double getMinY(){
        return this.theShape.getMinY();
    }
    
    public double getMaxY(){
        return this.theShape.getMaxY();
    }
    
    public double getFx_o(){return this.Fx_o;}
    
    public double getFx_n(){return this.Fx_n;}
    
    public double getFy_o(){return this.Fy_o;}
    
    public double getFy_n(){return this.Fy_n;}
    
    public double getMz_o(){return this.Mz_o;}
    
    public double getMz_n(){return this.Mz_n;}
    
    public void setFx_o(double v){Fx_o=v;}
    
    public void setFx_n(double v){Fx_n=v;}
    
    public void setFy_o(double v){Fy_o=v;}
    
    public void setFy_n(double v){Fy_n=v;}
    
    public void setMz_o(double v){Mz_o=v;}
    
    public void setMz_n(double v){Mz_n=v;}
    
    public void addFx_o(double v){Fx_o+=v;}
    
    public void addFx_n(double v){Fx_n+=v;}
    
    public void addFy_o(double v){Fy_o+=v;}
    
    public void addFy_n(double v){Fy_n+=v;}
    
    public void addMz_o(double v){Mz_o+=v;}
    
    public void addMz_n(double v){Mz_n+=v;}
    
    public void timesFx_o(double v){Fx_o*=v;}
    
    public void timesFx_n(double v){Fx_n*=v;}
    
    public void timesFy_o(double v){Fy_o*=v;}
    
    public void timesFy_n(double v){Fy_n*=v;}
    
    public void timesMz_o(double v){Mz_o*=v;}
    
    public void timesMz_n(double v){Mz_n*=v;}
    
    public void setMass(double m, double Jm){
        this.m=m;
        this.Jm=Jm;
    }
    
    public double get_m(){return m;}
    
    public double get_Jm(){return Jm;}
    
    public void setMaterial(Material aMat){this.theMaterial=aMat;}
}
