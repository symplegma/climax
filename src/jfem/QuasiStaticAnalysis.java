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
public class QuasiStaticAnalysis extends Analysis{
    private double[] Fext;
    private double[] DeltaU;
    
    // constructor
    public QuasiStaticAnalysis(){
        this.theType=1;
    }
    
    public QuasiStaticAnalysis(int typeofBCImposer){
        theDomain = new Domain();
        this.theType=1;
        theAlgorithm = new Linear(this, 0.001);
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        
    }
    
    public QuasiStaticAnalysis(Domain theDomain, int typeofBCImposer, Algorithm theAlgorithm){
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
    
    public QuasiStaticAnalysis(Domain theDomain, Algorithm theAlgorithm){
        this.theDomain = theDomain;
        this.theType=1;
        this.theAlgorithm=theAlgorithm;
        theAlgorithm.setAnalysis(this);
        theImposer = new StaticImposer();
        if(Analysis2File)this.printGenerals();
        
    }
    
    public QuasiStaticAnalysis(Domain theDomain){
        this.theDomain = theDomain;
        this.theType=1;
        theAlgorithm = new Linear(this, 10E-12);
        theAlgorithm.setAnalysis(this);
        theImposer = new StaticImposer();
        if(Analysis2File)this.printGenerals();
    }
    
    // methods 
    public int analyse(int LC) {
        int numSteps= theDomain.getLoadCase(LC).getNumOfIncrements();
        this.currentLC=LC;
//        this.setNodes(LC, numSteps);
        int log=0;
        this.step=0;
        int order=this.theDomain.getNdofs();
        double[] Res = new double[order];
        theSOE = new SOE(order);
        theDomain.clearDomain();
        //double coef=1./numSteps;
        
        this.Fext = new double[order];
        this.DeltaU = new double[order];
        
        // in the case where algortihm is initialNR or Linear
        if(this.theAlgorithm.getType()==4||this.theAlgorithm.getType()==1){
            this.formLeft();
        }
        if(Analysis2File){
            this.AnalysisFile.println();
            this.AnalysisFile.println("LoadCase with id: "+LC);
            this.AnalysisFile.println("______________________");
            this.AnalysisFile.println();
            this.AnalysisFile.println("Increment: "+0);
        }
        System.out.println("Elastic and Kinetic energy evolution");
        double elastic=0.5*theDomain.getuKu();
        double kinetic=0.5*theDomain.getvMv();
        double damped=0.;//this.timeStep*theDomain.getuCv(a0,a1);
        double surfEnrg=0.;
        double surfKnrg=0.;
        double surfEPnrg=0.;
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element theElement = it.next();
            if(theElement.getClass()==OneDOF.class || theElement.getClass()==Crack2d.class){
                surfEnrg+=0.5*theElement.getuKu();
                surfKnrg+=0.5*theElement.getvMv();
                surfEPnrg+=theElement.PlasticNRG();
                theElement.checkForDamage();
            }
        }
        
        double diss=theDomain.getDissipation();
        double extwork=0.;
        for(int j=0;j<Fext.length;j++){
            Fext[j]=this.theSOE.getFext()[j]/2.;
            DeltaU[j]=theDomain.getDisplacements()[j];
            extwork+=DeltaU[j]*Fext[j];
        }
        System.out.println(0+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+damped+" "+extwork);
        
        
        // print for trackers
            for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
//            theTracker.write("LoadCase with id: "+LC);
            theTracker.write(0.);
            }
        this.commit();
        // Increments
        this.formRight();
        for(int i=1; i<=numSteps; i++){
            this.step=i;
            this.formFext(LC,i, 1.0);
            for(int j=0;j<Fext.length;j++)Fext[j]=theSOE.getFext()[j]/2.;
            if(Analysis2File)this.AnalysisFile.println();
            if(Analysis2File)this.AnalysisFile.println("Increment: "+i);
            
            log=this.theAlgorithm.solve();
            
            surfEnrg=0.;
            surfKnrg=0.;
            surfEPnrg=0;
            int[] ids = new int[100];
            for(int j=0;j<ids.length;j++)ids[j]=0;
            int counter=0;
            boolean bool;
            for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
                Element theElement = it.next();
                if(theElement.getClass()==OneDOF.class || theElement.getClass()==Crack2d.class){
                    bool=false;
                    if(theElement.getDamage())bool=theElement.checkForDamage();
                    if(bool){
                        ids[counter]=theElement.getID(); counter+=1;
                    }
                    surfEnrg+=0.5*theElement.getuKu();
                    surfKnrg+=0.5*theElement.getvMv();
                    surfEPnrg+=theElement.PlasticNRG();
                }
            }
            diss=theDomain.getDissipation();
            extwork=0.;
            for(int j=0;j<Fext.length;j++){
                DeltaU[j]=theDomain.getDisplacements()[j]-DeltaU[j];
                extwork+=DeltaU[j]*Fext[j];
                DeltaU[j]=theDomain.getDisplacements()[j];
            }
            for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
                Element theElement = it.next();
                if(theElement.getDamage())theElement.commitDamage();
            }
            elastic=0.5*theDomain.getuKu();
            kinetic=0.5*theDomain.getvMv();
            damped=theDomain.getuCv(a0,a1);//this.timeStep* : used in combination with velocity calc of damp energy
            
            System.out.print(i+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+damped+" "+extwork+" ");
            for(int j=0; j<ids.length;j++){
                if(ids[j]!=0){System.out.print(ids[j]+" ");}
            }
            System.out.println();
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
            // print for trackers
            for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
//            theTracker.write(i);
            if(i % theTracker.getStep() == 0)theTracker.write(i);
            }
            
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
        
        aString="Type of Analysis: "+"PseudoDynamic Static Analysis";
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        switch(this.theAlgorithm.getType()){
            case 1: aString="Type of Algorithm: "+"Linear"; break;
            case 2: aString="Type of Algorithm: "+"Full Newton-Raphson"; break;
            case 3: aString="Type of Algorithm: "+"Modified Newton-Raphson"; break;
            case 4: aString="Type of Algorithm: "+"Initial Newton-Raphson"; break;
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

