/*******************************************************************************
* Climax.                                                                      *
* Copyright (C) 2009-2017 C.G. Panagiotopoulos [http://www.symplegma.org]      *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): C.G. Panagiotopoulos (pchr76@gmail.com)
// *****************************************************************************
package jfem;
import java.util.Iterator;
import java.io.*;
import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;

import java.util.Map;
import java.util.TreeMap;
import javax.swing.JTextArea;
import jmat.AbstractMatrix;


/**
 *
 * @author pchr
 */
abstract public class Analysis {
    protected int id;
    protected static int numOfAnalyses;
    protected Domain theDomain;
    protected BCImposer theImposer;
    protected SOE theSOE;
    protected Algorithm theAlgorithm;
    protected PrintWriter AnalysisFile;
    protected int theType;
    protected int step=0;
    protected boolean RayleighDamping=false;
    protected double a0=0.0,a1=0.0;
    
    private double c_f1;
    private double c_f2;
    private int imp_f;
    
    protected double timeStep;
    protected int currentLC;
    protected static boolean  Analysis2File=false;
    
    protected Map<Integer,Tracker> theTrackers = new TreeMap<Integer,Tracker>();
    
    // constructor
    public Analysis(){
        c_f1=1.0;
        c_f2=0.0;
        try {
            id = ++numOfAnalyses;
            String filename="AnalysisFile_"+id+".jaf";
            if(Analysis2File)AnalysisFile = new PrintWriter(new FileWriter(filename));
        } catch (IOException ex) {
            Logger.getLogger(Analysis.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    // methods
    public void putTracker(Tracker aTracker) {
        theTrackers.put(aTracker.getID(), aTracker);
    }
    
    public void set_coefs(double f1, double f2){
        this.c_f1=f1;
        this.c_f2=f2;
    }
    
    public Domain getDomain(){
        return this.theDomain;
    }
    
    public SOE getSOE(){
        return this.theSOE;
    }
    
    public int getType(){
        return this.theType;
    }
    
    public int getStep(){
        return this.step;
    }
    
    public double getTimeStepLenght(){
        return this.timeStep;
    }
    
    public void setAlgorithm(Algorithm theAlgorithm){
        this.theAlgorithm=theAlgorithm;
    }
    
    public void setImposer(BCImposer theImposer){
        this.theImposer=theImposer;
    }
    
    protected void formStiffnessMatrix(){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getK();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToK(mat, aInt, 1.0);
        }
    }
    
    protected void formMassMatrix(){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToM(mat, aInt, 1.0);
        } 
    }
    
    protected void formFext(int LCase,int time, double coef){
        this.theSOE.clearFext();
        
        LoadCase LC =theDomain.getLoadCase(LCase);
        
        for(Iterator<Load> lit=LC.getLoads().values().iterator(); lit.hasNext();){
            Load theLoad = lit.next();
            int LoadTime = theLoad.getTime();
            if(theLoad.getClass().getName().compareTo(LoadDist.class.getName())==0){
                if(LoadTime==time){
                    int we = ((LoadDist)theLoad).getwelem();
                    AbstractMatrix vec = theDomain.getElement(we).getFequivalent((LoadDist) theLoad);
                    int[] aInt= theDomain.getElement(we).getFtable();
                    theSOE.addToFext(vec, aInt, this.c_f1);
                }
            }else{
                if(LoadTime==time){
                    double val = theLoad.getLoadValue();
                    int wn = theLoad.getwnode();
                    int wd = theLoad.getwdof();
                    AbstractMatrix vec = new AbstractMatrix(1,1);
                    int[] aInt= new int[1];
                    vec.putVal(0, 0, val*coef);
                    aInt[0]=theDomain.getNode(wn).getFtable()[wd-1];
                    theSOE.addToFext(vec, aInt, this.c_f1);
                }
            }
        }
        
        for(Iterator<Load> lit=LC.getLoads().values().iterator(); lit.hasNext();){
            Load theLoad = lit.next();
            int LoadTime = theLoad.getTime();
            if(theLoad.getClass().getName().compareTo(LoadDist.class.getName())==0){
                if(LoadTime==time-1){
                    int we = ((LoadDist)theLoad).getwelem();
                    AbstractMatrix vec = theDomain.getElement(we).getFequivalent((LoadDist) theLoad);
                    int[] aInt= theDomain.getElement(we).getFtable();
                    theSOE.addToFext(vec, aInt, this.c_f2);
                }
            }else{
                if(LoadTime==time-1){
                    double val = theLoad.getLoadValue();
                    int wn = theLoad.getwnode();
                    int wd = theLoad.getwdof();
                    AbstractMatrix vec = new AbstractMatrix(1,1);
                    int[] aInt= new int[1];
                    vec.putVal(0, 0, val*coef);
                    aInt[0]=theDomain.getNode(wn).getFtable()[wd-1];
                    theSOE.addToFext(vec, aInt, this.c_f2);
                }
            }
        }
        
        theImposer.imposeF_ext(this, time);
    }
    
    protected void formFint(){
        this.theSOE.clearFint();
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getF();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFint(mat, aInt, 1.0);
        }
        theImposer.imposeF_int(this);
    }
    
    public int getID(){
        return this.id;
    }
    
    public int analyse(){return analyse(1);}
    
    abstract public int analyse(int LC);
    abstract public void formLeft();
    abstract public void formRight();
    abstract public void run();
    abstract public void update();
    abstract public void commit();
    abstract public void printGenerals();
    //abstract public void endAnalysis();
    
    //protected void setNodes(int LC, int numSteps){
//        this.theDomain.
    //}
    
    public void endAnalysis() throws SQLException {
        if(Analysis2File)this.AnalysisFile.close();
        
        for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
        }
        
        
        for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
            theTracker.close();
        }
        
        Node aNode = new Node(); aNode.initDOFS();
        
//        this.theSOE.numDOFS=0;
//        this.theSOE=null;

    }
    
    public Map getTrackers(){
        return this.theTrackers;
    }
    
    public Tracker getTracker(int id_){
        return theTrackers.get(id_);
    }

    public int getNumTrackers(){
        return this.theTrackers.size();
    }
    
    public void setRayleighDamping(double a0, double a1){
        this.a0=a0;
        this.a1=a1;
        RayleighDamping=true;
    }
}
