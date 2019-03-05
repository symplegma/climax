/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class City {
    private int id;
    private Map<Integer,Block> theBlocks;
    private static int Capacity=100; // rerfers to Particles
    
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    
    public City(int id){
        this.id=id;
        theBlocks = new TreeMap<Integer,Block>();
    }
    
    public void putBlock(Block aBlock){
        this.theBlocks.put(aBlock.getID(), aBlock);

        if(aBlock.getMaxX()>maximum_X_coordinate)maximum_X_coordinate=aBlock.getMaxX();
        if(aBlock.getMaxY()>maximum_Y_coordinate)maximum_Y_coordinate=aBlock.getMaxY();
        if(aBlock.getMinX()<minimum_X_coordinate)minimum_X_coordinate=aBlock.getMinX();
        if(aBlock.getMinY()<minimum_Y_coordinate)minimum_Y_coordinate=aBlock.getMinY();
        
        //if(this.getPopulation()>Capacity)System.err.println("Warning: Capacity for City ("+id+") has been exceeded");
    }
    
    public static void setCityCapacity(int c){Capacity=c;}
    
    public static int getCityCapacity(){return Capacity;}
    
    public int getPopulation(){
        int pop=0;
        for(Iterator<Block> cit=this.theBlocks.values().iterator(); cit.hasNext();){
            Block theBlock = cit.next();
            pop+=theBlock.getPopulation();
        }
        return pop;
    }
    
    public Map<Integer,Block> getBlocks(){return this.theBlocks;}
    
    public int getID(){
        return id;
    }
    
    public double getMaxX(){return this.maximum_X_coordinate;}
    
    public double getMaxY(){return this.maximum_Y_coordinate;}
    
    public double getMinX(){return this.minimum_X_coordinate;}
    
    public double getMinY(){return this.minimum_Y_coordinate;}
}
