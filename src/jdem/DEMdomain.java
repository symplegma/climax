/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import gendomain.Node;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import mater.Material;

/**
 *
 * @author antonia
 */
public class DEMdomain {
    private int id;
    protected Map<Integer,Particle> theParticles = new TreeMap<Integer,Particle>();
    protected Map<Integer,Node> theNodes = new TreeMap<Integer,Node>();
    private Map<Integer,InitialCondition> theInitialConditions = new TreeMap<Integer,InitialCondition>();
    private Map<Integer,City> theCities= new TreeMap<Integer,City>();
    private Map<Integer,Material> theMaterials = new TreeMap<Integer,Material>();
    private static int numberOfDomains = 0;
    
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Z_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Z_coordinate=Double.POSITIVE_INFINITY;
    
    public DEMdomain(){
        id = ++numberOfDomains;
    }
    
    public DEMdomain(int nCities, int nBlocks, int nCells){
        id = ++numberOfDomains;
    }

    public int getID(){
        return this.id;
    }
    
    public void clsNumberOfDomains(){numberOfDomains=0;}
    
    public void putNode(Node aNode){
        if(aNode.getCoordinates()[0]>maximum_X_coordinate)maximum_X_coordinate=aNode.getCoordinates()[0];
        if(aNode.getCoordinates()[0]<minimum_X_coordinate)minimum_X_coordinate=aNode.getCoordinates()[0];
        if(aNode.getCoordinates().length>1){
            if(aNode.getCoordinates()[1]>maximum_Y_coordinate)maximum_Y_coordinate=aNode.getCoordinates()[1];
            if(aNode.getCoordinates()[1]<minimum_Y_coordinate)minimum_Y_coordinate=aNode.getCoordinates()[1];
        }
        if(aNode.getCoordinates().length>2){
            if(aNode.getCoordinates()[2]>maximum_Z_coordinate)maximum_Z_coordinate=aNode.getCoordinates()[2];
            if(aNode.getCoordinates()[2]<minimum_Z_coordinate)minimum_Z_coordinate=aNode.getCoordinates()[2];
        }
        this.theNodes.put(aNode.getID(), aNode);
    }
    
    public void putCity(City Rethymno){
        this.theCities.put(Rethymno.getID(), Rethymno);
    }
    
    public void setParticles(){
        for(Iterator<City> city=this.theCities.values().iterator(); city.hasNext();){
            City theCity = city.next();
            for(Iterator<Block> bit=theCity.getBlocks().values().iterator(); bit.hasNext();){
                Block theBlock = bit.next();
                for(Iterator<Cell> cit=theBlock.getCells().values().iterator(); cit.hasNext();){
                    Cell theCell = cit.next();
                    for(Iterator<Particle> pit=theCell.getParticles().values().iterator(); pit.hasNext();){
                        Particle aParticle=pit.next();
                        this.theParticles.put(aParticle.getID(), aParticle);
                        this.putNode(aParticle.getCentroid());
                        if(aParticle.getMaxX()>maximum_X_coordinate)maximum_X_coordinate=aParticle.getMaxX();
                        if(aParticle.getMaxY()>maximum_Y_coordinate)maximum_Y_coordinate=aParticle.getMaxY();
                        if(aParticle.getMinX()<minimum_X_coordinate)minimum_X_coordinate=aParticle.getMinX();
                        if(aParticle.getMinY()<minimum_Y_coordinate)minimum_Y_coordinate=aParticle.getMinY();
                    }
                }
            }
        }
        
    }
    
    public Node getNode(int id){
        return this.theNodes.get(id);
    }

    public Map<Integer,Node> getNodes(){
        return this.theNodes;
    }
    
    public Particle getParticle(int id){
        return this.theParticles.get(id);
    }
    
    public Map<Integer,Particle> getParticles(){
        return this.theParticles;
    }
    
    public City getCity(int id){
        return this.theCities.get(id);
    }
    
    public Map<Integer,City> getCities(){
        return this.theCities;
    }
    
    public void putInitialCondition(InitialCondition aInitialCondition) {
        theInitialConditions.put(aInitialCondition.getID(), aInitialCondition);
    }
    
    public InitialCondition getInitialCondition(int id_){
        return theInitialConditions.get(id_);
    }
    
    public Map getInitialConditions(){
        return this.theInitialConditions;
    }
    
    public void DomainInitialConditions(){
        for(Iterator<InitialCondition> it=this.theInitialConditions.values().iterator(); it.hasNext();){
            InitialCondition theInitialCondition = it.next();
            switch(theInitialCondition.getType()){
                case 1: theInitialCondition.getNode().initInitialialConditionsTD(theInitialCondition.getInitialValue(), theInitialCondition.getDOF()-1); break;
                default: theInitialCondition.getNode().initInitialialConditions(theInitialCondition.getInitialValue(), theInitialCondition.getDOF()-1); break;
                        
            }
        }
    }
    
    public double getMaxX(){return this.maximum_X_coordinate;}
    
    public double getMaxY(){return this.maximum_Y_coordinate;}
    
    public double getMaxZ(){return this.maximum_Z_coordinate;}
    
    public double getMinX(){return this.minimum_X_coordinate;}
    
    public double getMinY(){return this.minimum_Y_coordinate;}
    
    public double getMinZ(){return this.minimum_Z_coordinate;}
    
    public double maximum_X_coordinate(){return getMaxX();}
    
    public double maximum_Y_coordinate(){return getMaxY();}
    
    public double minimum_X_coordinate(){return getMinX();}
    
    public double minimum_Y_coordinate(){return getMinY();}
    
    public void printCities(){
        System.out.println(this.getClass().toString()+" with ID= "+this.id+" ("+this.theParticles.size()+" particles)");
        for (City theCity : this.theCities.values()) {
            System.out.println("    City: "+theCity.getID()+" ("+theCity.getPopulation()+" particles)");
            for (Block theBlock : theCity.getBlocks().values()) {
                System.out.println("        Block: "+theBlock.getID()+" ("+theBlock.getPopulation()+" particles)");
                theBlock.getCells().values().stream().map((theCell) -> {
                    System.out.println("            Cell: "+theCell.getID()+" ("+theCell.getPopulation()+" particles)");
                    return theCell;
                }).forEach((theCell) -> {
                    theCell.getParticles().values().stream().forEach((theP) -> {
                        System.out.println("                Particle: "+theP.getID());
                    });
                });
            }
        }
    }
    
    public void printParticles(){
        System.out.println(this.getClass().toString()+" with ID= "+this.id);
        this.theParticles.values().stream().forEach((theP) -> {
            System.out.println("Particle: "+theP.getID());
        });
    }
    
    public void putMaterial(Material aMaterial) {
        theMaterials.put(aMaterial.getID(), aMaterial);
    }
    
    public Material getMaterial(int id_){
        return theMaterials.get(id_);
    }
}
