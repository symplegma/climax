/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

/**
 *
 * @author pchr
 */
abstract public class Analysis {
    protected int id;
    protected static int numOfAnalyses;
    protected DEMdomain theDomain;
    protected BCImposer theImposer;
    protected int step=0;
    protected double timeStep;
    protected int numSteps;
    
    public Analysis(){
        id = ++numOfAnalyses;
    }
    
    public int getID(){
        return this.id;
    }
    
    public void setImposer(BCImposer theImposer){
        this.theImposer=theImposer;
    }
    
    public int getStep(){
        return this.step;
    }
    
    public double getTimeStepLenght(){
        return this.timeStep;
    }
    
    public void setNumSteps(int num){this.numSteps=num;}
    
    public void setTimeStep(double dt){this.timeStep=dt;}
    
    abstract public int analyse();
}
