/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Cell {
    private int id;
    private Map<Integer,Particle> theParticles;
    private static int Capacity=10; // rerfers to Particles
    private int Population=0; // rerfers to Particles
    
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    
    public Cell(int id){
        this.id=id;
        theParticles = new TreeMap<Integer,Particle>();
    }
    
    public static void setCellCapacity(int c){Capacity=c;}
    
    public static int getCellCapacity(){return Capacity;}
    
    public int getPopulation(){return this.Population;}
    
    public void putParticle(Particle aParticle){
        this.theParticles.put(aParticle.getID(), aParticle);
        Population++;

        if(aParticle.getMaxX()>maximum_X_coordinate)maximum_X_coordinate=aParticle.getMaxX();
        if(aParticle.getMaxY()>maximum_Y_coordinate)maximum_Y_coordinate=aParticle.getMaxY();
        if(aParticle.getMinX()<minimum_X_coordinate)minimum_X_coordinate=aParticle.getMinX();
        if(aParticle.getMinY()<minimum_Y_coordinate)minimum_Y_coordinate=aParticle.getMinY();
    }
    
    public Map<Integer,Particle> getParticles(){return this.theParticles;}
    
    public int getID(){
        return id;
    }
    
    public double getMaxX(){return this.maximum_X_coordinate;}
    
    public double getMaxY(){return this.maximum_Y_coordinate;}
    
    public double getMinX(){return this.minimum_X_coordinate;}
    
    public double getMinY(){return this.minimum_Y_coordinate;}
}
