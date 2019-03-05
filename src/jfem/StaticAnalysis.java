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


/**
 *
 * @author pchr
 */
public class StaticAnalysis extends Analysis{
    
    // constructor
    public StaticAnalysis(){
        this.theType=1;
    }
    
    public StaticAnalysis(Domain theDomain){
        this.theDomain = theDomain;
        this.theType=1;
        theAlgorithm = new Linear(this, 10E-12);
        theImposer = new StaticImposer();
    }
    
    public StaticAnalysis(int typeofBCImposer){
        theDomain = new Domain();
        this.theType=1;
        theAlgorithm = new Linear(this, 0.001);
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        
    }
    
    public StaticAnalysis(Domain theDomain, int typeofBCImposer, Algorithm theAlgorithm){
        this.theDomain = theDomain;
        this.theType=1;
        this.theAlgorithm=theAlgorithm;
        theAlgorithm.setAnalysis(this);
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        if(Analysis2File)this.printGenerals();
    }
    
    public StaticAnalysis(Domain theDomain, Algorithm theAlgorithm){
        this.theDomain = theDomain;
        this.theType=1;
        this.theAlgorithm=theAlgorithm;
        theAlgorithm.setAnalysis(this);
        theImposer = new StaticImposer();
        if(Analysis2File)this.printGenerals();
    }
    
    // methods 
    public int analyse(int LC) {
        int numSteps= theDomain.getLoadCase(LC).getNumOfIncrements();
        //this.setNodes(LC, numSteps);
        this.currentLC=LC;
        int log=0;
        this.step=0;
        int order=this.theDomain.getNdofs();
        double[] Res = new double[order];
        theSOE = new SOE(order);
        theDomain.clearDomain();
        double coef=1./numSteps;
        
        // in the case where algortihm is initialNR or Linear
        if(this.theAlgorithm.getType()==4||this.theAlgorithm.getType()==1){
            this.formLeft();
        }
//        this.theSOE.printK();
        if(Analysis2File){
            this.AnalysisFile.println();
            this.AnalysisFile.println("LoadCase with id: "+LC);
            this.AnalysisFile.println("______________________");
            this.AnalysisFile.println();
            this.AnalysisFile.println("Increment: "+0);
        }
        for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
//            theTracker.write("LoadCase with id: "+LC);
        }
        this.commit();
        // Increments
        this.formRight();
        for(int i=1; i<=numSteps; i++){
            this.step=i;
            this.formFext(LC,0, i*coef);
//            this.theSOE.printF();
            
            if(Analysis2File)this.AnalysisFile.println();
            if(Analysis2File)this.AnalysisFile.println("Increment: "+i);
            
            log=this.theAlgorithm.solve();
            if(Analysis2File){
                this.AnalysisFile.println();
                this.AnalysisFile.print("number of plastified elements: ");
                this.AnalysisFile.println(this.theDomain.getNumOfPlastifiedElems());
                this.AnalysisFile.println("log of increment "+i+" is "+log);
            }
            if(log==1){
                System.err.println("warning: log of increment "+i+" is "+log);
            }
            if(Analysis2File)this.AnalysisFile.println("iteration for this step: "+theAlgorithm.getIteraions());
            
        }
        if(Analysis2File)this.AnalysisFile.println();
        return log;
    }
    
    @Override
    public void run(){
        this.theSOE.StaticRun();
    }
    
    @Override
    public void update(){
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateDispsStep(theSOE.getX(),1.0);
        }
    }
    
    @Override
    public void commit(){
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitDisps();
            theNode.commitLoadCaseDisp(this.currentLC, step);
            if(Analysis2File)theNode.printDisps(this);
        }
        // commit sta elements
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.commit();
            //theElement.printAccls(this);
        }
        
        // print for trackers
        for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
            theTracker.write(this.step);
        }
    }

    @Override
    public void formLeft() {
        theSOE.clearK();
        this.formStiffnessMatrix();
        theImposer.imposeLeft(this);
    }

    @Override
    public void formRight() {
        this.formFint();
    }    
    
    
    public void printGenerals() {
        String aString="Analysis with id: "+this.id;
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        aString="Type of Analysis: "+"Static Analysis";
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        switch(this.theAlgorithm.getType()){
            case 1: aString="Type of Algorithm: "+"Linear"; break;
            case 2: aString="Type of Algorithm: "+"Full Newron-Raphson"; break;
            case 3: aString="Type of Algorithm: "+"Modified Newron-Raphson"; break;
            case 4: aString="Type of Algorithm: "+"Initial Newron-Raphson"; break;
            default: aString="Type of Algorithm: "+"Unknown"; break;
        }
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        switch(this.theImposer.getType()){
            case 1: aString="Type of Boundary Conditions Imposer: "+"Penalty"; break;
            default: aString="Type of Boundary Conditions Imposer: "+"Unknown"; break;
        }
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        aString="Domain of analysis id: "+this.theDomain.getID();
        if(Analysis2File)this.AnalysisFile.println(aString);
        aString="Number of degrees of freedom: "+this.theDomain.getNdofs();
        if(Analysis2File)this.AnalysisFile.println(aString);
        aString="Number of Load Cases: "+this.theDomain.getLoadCases().size();
        if(Analysis2File)this.AnalysisFile.println(aString);
        if(Analysis2File)this.AnalysisFile.println("===================================================");
    }

}
