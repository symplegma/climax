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
import jmat.AbstractMatrix;
/**
 *
 * @author pchr
 */
public class SimoTransientAnalysis extends Analysis{
    private double[] Fext;
    private double[] DeltaU;
    // constructor
    public SimoTransientAnalysis(){
        this.theType=4;
    }
    
    public SimoTransientAnalysis(double timeStep){
        theDomain = new Domain();
        this.theType=4;
        this.timeStep=timeStep;
        theImposer = new SimoImposer();
    }
    
    public SimoTransientAnalysis(Domain theDomain, double timeStep){
        this.theDomain = theDomain;
        this.theType=4;
        this.timeStep=timeStep;
        theImposer = new SimoImposer();
    }
    
    public SimoTransientAnalysis(Domain theDomain, Algorithm theAlgorithm, 
                               double timeStep){
        this.theDomain = theDomain;
        this.theType=4;
        this.timeStep=timeStep;
        this.theAlgorithm=theAlgorithm;
        theAlgorithm.setAnalysis(this);
        theImposer = new SimoImposer();
    }

    @Override
    public int analyse(int LC) {
        int numSteps= theDomain.getLoadCase(LC).getNumOfIncrements();
        this.currentLC=LC;
//        this.setNodes(LC, numSteps);
        int log=0;
        int order=this.theDomain.getNdofs();
        double[] Res = new double[order];
        theSOE = new SOE(order,order);
        theDomain.clearDomain();
        theDomain.clearDeltaDisps();
        // configure initial conditions for 
        // displacements and velocities
        theDomain.DomainInitialConditions();
        this.set_coefs(1., 1.);
        this.formLeft();
        this.formFext(LC, 0, 1.0);
        this.formFint();
        
//        theSOE.TransientRunInit(); // to eixa comment
//        this.updateInitialAcceleration(); // to eixa comment
        
        this.Fext = new double[order];
        this.DeltaU = new double[order];
        if(Analysis2File){
            this.AnalysisFile.println();
            this.AnalysisFile.println("LoadCase with id: "+LC);
            this.AnalysisFile.println("______________________");
            this.AnalysisFile.println("time step: "+0);
        }
        this.commit2();
        if(Analysis2File)this.AnalysisFile.print("number of plastified elements: ");
        if(Analysis2File)this.AnalysisFile.println(this.theDomain.getNumOfPlastifiedElems());
        System.out.println("Elastic and Kinetic energy evolution");
        double elastic=0.5*theDomain.getuKu();
        double kinetic=0.5*theDomain.getvMv();
        double damped=this.timeStep*theDomain.getuCv(a0,a1);
        double surfEnrg=0.;
        double surfKnrg=0.;
        double constraintNRG=0.5*theDomain.getuKu_constraints();
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element theElement = it.next();
            if(theElement.getClass()==OneDOF.class || theElement.getClass()==Crack2d.class){
                surfEnrg+=0.5*theElement.getuKu();
                surfKnrg+=0.5*theElement.getvMv();
                theElement.checkForDamage();
            }
        }
        double diss=theDomain.getDissipation();
        double extwork=0.;
        for(int j=0;j<Fext.length;j++){
            Fext[j]=this.theSOE.getFext()[j];
            DeltaU[j]=theDomain.getDisplacements()[j];
            extwork+=DeltaU[j]*Fext[j];
        }
//        extwork = theDomain.getuKu_constraints();
        System.out.println(0+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+damped+" "+extwork+" "+constraintNRG);
//        this.printTextArea.append(0+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+damped+" "+extwork+" "+constraintNRG+'\n');
        // print for trackers
            for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
//            theTracker.write("LoadCase with id: "+LC);
            theTracker.write(0.);
            }
        
        // TimeSteps
        this.formRight();
        for(int i=1; i<=numSteps; i++){
            this.theDomain.clearDeltaDisps();
            this.step=i;
            this.ClearAccumulatedDips();
            this.formFext(LC,i, 1.0);
            for(int j=0;j<Fext.length;j++)Fext[j]=theSOE.getFext()[j]/2.;
            this.addM_u(4./(timeStep*timeStep));
            this.addM_v(4./timeStep);
            this.addK_u(-1.);
            
            if(this.RayleighDamping){
                this.addC_u(2./timeStep);
            }

            if(Analysis2File)this.AnalysisFile.println();
            if(Analysis2File)this.AnalysisFile.println("time step: "+i);

            log=this.theAlgorithm.solve();
            
            surfEnrg=0.;
            surfKnrg=0.;
            int[] ids = new int[100];
            for(int j=0;j<ids.length;j++)ids[j]=0;
            int counter=0;
            boolean bool;
            for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
                Element theElement = it.next();
                if(theElement.getClass()==OneDOF.class || theElement.getClass()==Crack2d.class){
                    surfEnrg+=0.5*theElement.getuKu();
                    surfKnrg+=0.5*theElement.getvMv();
                    bool=false;
                    if(theElement.getDamage())bool=theElement.checkForDamage();
                    if(bool){
                        ids[counter]=theElement.getID(); counter+=1;
                    }
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
            constraintNRG=0.5*theDomain.getuKu_constraints();
            
            System.out.print(i+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+damped+" "+extwork+" "+constraintNRG+" ");
//            this.printTextArea.append(i+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+damped+" "+extwork+" "+constraintNRG+" ");
            for(int j=0; j<ids.length;j++){
                if(ids[j]!=0){
                    System.out.print(ids[j]+" "); //this.printTextArea.append(ids[j]+" ");
                }
            }
            System.out.println();
//            this.printTextArea.append(" "+'\n');
            if(Analysis2File)this.AnalysisFile.print("number of plastified elements: ");
            if(Analysis2File)this.AnalysisFile.println(this.theDomain.getNumOfPlastifiedElems());
            if(Analysis2File)this.AnalysisFile.println("log of increment "+i+" is "+log);
            if(log==1){
                System.err.println("warning: log of increment "+i+" is "+log);
            }
            if(Analysis2File)this.AnalysisFile.println("iteration for this step: "+theAlgorithm.getIteraions());
            
            // print for trackers
            for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
                Tracker theTracker = it.next();
                if(i % theTracker.getStep() == 0)theTracker.write(i*timeStep);
            }
        }
        if(Analysis2File)this.AnalysisFile.println();
        return log;
    }
     
     protected void formVelocity(){
         theSOE.clearXv();
         for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
             Node theNode = it.next();
             AbstractMatrix mat ;
             mat=theNode.getVelcsTrial_mat();
             int[] aInt ;
             aInt=theNode.getFtable();
             theSOE.addToXv(mat, aInt,1.0);
         } 
    }
     
     protected void formAcceleration(){
         theSOE.clearXa();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            AbstractMatrix mat ;
            mat=theNode.getAcclsTrial_mat();
            int[] aInt ;
            aInt=theNode.getFtable();
            theSOE.addToXa(mat, aInt,1.0);
        } 
    }
     
     protected void formDisplacement(){
         theSOE.clearXu();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            AbstractMatrix mat ;
            //mat=theNode.getDispsConvg_mat();
            mat=theNode.getDispsAccumulated_mat();
            int[] aInt ;
            aInt=theNode.getFtable();
            theSOE.addToXu(mat, aInt,1.0);
        } 
    }
     
     public void ClearAccumulatedDips(){
         for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.clearAccumulatedDips();
        } 
     }

    @Override
    public void formLeft() {
        theSOE.clearK();
        theSOE.clearM();
        this.theImposer.setCoefficient_k(1.0);
        this.theImposer.setCoefficient_m(1.0);
        this.formStiffnessMatrix();
        this.formMassMatrix();
        theImposer.imposeLeft(this);
    }

    @Override
    public void formRight() {
        this.formFintEff();
    }
    
    protected void addM_v(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM_v();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getM_v();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            }
        }
        
    }
    
    protected void addK_u(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            //mat=elem.getK_u();
            mat=elem.getF();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                //mat=elem.getK_u();
                mat=elem.getF();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            } 
        }

    }
    
    protected void addM_u(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM_u();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getM_u();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            } 
        }
    }
    
    protected void addC_u(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getM_u();
            theSOE.addToFext(mat, aInt, coef*a0);
            mat=elem.getK_u();
            theSOE.addToFext(mat, aInt, coef*a1);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                int[] aInt ;
                aInt=elem.getFtable();
                mat=elem.getM_u();
                theSOE.addToFext(mat, aInt, coef*a0);
                mat=elem.getK_u();
                theSOE.addToFext(mat, aInt, coef*a1);
            } 
        }
    }

    @Override
    public void run() {
        if(this.RayleighDamping){
            this.theSOE.TransientRun(1.0, 4.0/( timeStep*timeStep ), 2.0/( timeStep ), a0, a1 );
        }else{
            this.theSOE.TransientRun(1.0, 4.0/( timeStep*timeStep ));
        }
    }

    @Override
    public void update(){
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateDispsStep(theSOE.getX(),1.0);
            theNode.accumulateDisps(theSOE.getX());
        }
    }
    
    //===================================================================
    public void updateVelc(){
        this.formVelocity();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateVelcsStep( theSOE.getXu(), 2./this.timeStep);
            theNode.updateVelcsStep( theSOE.getXv(), -2. );
            //theNode.updateVelcsStep( theSOE.getXa(), 1. );
        }
    }
    
    public void updateAccl(){
        this.formAcceleration();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateAcclsStep( theSOE.getXu(), 1.);
            theNode.updateAcclsStep( theSOE.getXv(), 1. );
            theNode.updateAcclsStep( theSOE.getXa(), -1. );
        }
    }
    //===================================================================
    

    @Override
    public void commit(){
        this.formDisplacement();
        this.updateVelc();
        this.updateAccl();
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitDisps();
            theNode.commitLoadCaseDisp(this.currentLC, step);
            if(Analysis2File)theNode.printDisps(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitVelcs();
            if(Analysis2File)theNode.printVelcs(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitAccls();
            if(Analysis2File)theNode.printAccls(this);
        }
        
        // commit sta elements
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.commit();
            //theElement.printAccls(this);
        }
    }
    
    public void commit2(){
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitDisps();
            theNode.commitLoadCaseDisp(this.currentLC, step);
            if(Analysis2File)theNode.printDisps(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitVelcs();
            if(Analysis2File)theNode.printVelcs(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitAccls();
            if(Analysis2File)theNode.printAccls(this);
        }
    }
    
    private void updateInitialAcceleration(){
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateAccls(theSOE.getX());
        }
    }
    
    @Override
    public void printGenerals() {
        String aString="Analysis with id: "+this.id;
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        aString="Type of Analysis: "+"b-Newmark Transient Analysis";
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
    
    protected void formFintEff(){
        this.theSOE.clearFint();
        
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getF();
            theSOE.addToFint(mat, aInt, 1.0);
            mat=elem.getM_u();
            mat=mat.times(4./(this.timeStep*this.timeStep));
            theSOE.addToFint(mat, aInt, 1.0);
            if(this.RayleighDamping){
                mat=elem.getM_u();
                mat=mat.times(2.*a0/this.timeStep);
                theSOE.addToFint(mat, aInt, 1.0);
                mat=elem.getK_u();
                mat=mat.times(2.*a1/this.timeStep);
                theSOE.addToFint(mat, aInt, 1.0);
            }
        }
        for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
            ConstraintElement elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getF();
            theSOE.addToFint(mat, aInt, 1.0);
            if (this.theImposer.TypeOfImposer!=1){
                mat=elem.getM_u();
                mat=mat.times(4./(this.timeStep*this.timeStep));
                theSOE.addToFint(mat, aInt, 1.0);
                if(this.RayleighDamping){
                    mat=elem.getM_u();
                    mat=mat.times(2.*a0/this.timeStep);
                    theSOE.addToFint(mat, aInt, 1.0);
                    mat=elem.getK_u();
                    mat=mat.times(2.*a1/this.timeStep);
                    theSOE.addToFint(mat, aInt, 1.0);
                }
            }
            
        }

    }

}
