/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Block {
    private int id;
    private Map<Integer,Cell> theCells;
    private static int Capacity=20; // rerfers to Particles
    
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    
    public Block(int id){
        this.id=id;
        theCells = new TreeMap<Integer,Cell>();
    }
    
    public static void setBlockCapacity(int c){Capacity=c;}
    
    public static int getBlockCapacity(){return Capacity;}
    
    public int getPopulation(){
        int pop=0;
        for(Iterator<Cell> cit=this.theCells.values().iterator(); cit.hasNext();){
            Cell theCell = cit.next();
            pop+=theCell.getPopulation();
        }
        return pop;
    }
    
    public void putCell(Cell aCell){
        this.theCells.put(aCell.getID(), aCell);

        if(aCell.getMaxX()>maximum_X_coordinate)maximum_X_coordinate=aCell.getMaxX();
        if(aCell.getMaxY()>maximum_Y_coordinate)maximum_Y_coordinate=aCell.getMaxY();
        if(aCell.getMinX()<minimum_X_coordinate)minimum_X_coordinate=aCell.getMinX();
        if(aCell.getMinY()<minimum_Y_coordinate)minimum_Y_coordinate=aCell.getMinY();
    }
    
    public Map<Integer,Cell> getCells(){return this.theCells;}
    
    public Map<Integer,Particle> getParticles(){
        Map<Integer, Particle> BlockParticles=new TreeMap<Integer, Particle>();
        for(Iterator<Cell> cit=this.theCells.values().iterator(); cit.hasNext();){
            Cell theCell = cit.next();
            BlockParticles.putAll(theCell.getParticles());
        }
        return BlockParticles;
    }
    
    public int getID(){
        return id;
    }
    
    public double getMaxX(){return this.maximum_X_coordinate;}
    
    public double getMaxY(){return this.maximum_Y_coordinate;}
    
    public double getMinX(){return this.minimum_X_coordinate;}
    
    public double getMinY(){return this.minimum_Y_coordinate;}
}
