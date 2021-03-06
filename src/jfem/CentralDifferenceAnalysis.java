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
public class CentralDifferenceAnalysis extends TimeIntegration{
    private double[] w;
    
    public CentralDifferenceAnalysis(){}
    
    public CentralDifferenceAnalysis(Domain theDomain, int typeofBCImposer, 
            Algorithm theAlgorithm, double timeStep){
        this.theDomain = theDomain;
        this.timeStep=timeStep;
        this.theAlgorithm=theAlgorithm;
        theAlgorithm.setAnalysis(this);
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            case 2: theImposer = new ClassicLargeMass(); break;
            case 3: theImposer = new GeneralizedLargeMass(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        w = new double[4];
        w[0]=1./(timeStep*timeStep);   w[1]=1.0/(2.0*timeStep); 
        w[2]=2.0*w[0];   w[3]=1.0/w[2];
    }

    @Override
    public int analyse(int LC) {
        int numSteps= theDomain.getLoadCase(LC).getNumOfIncrements();
//        this.setNodes(LC, numSteps);
        this.currentLC=LC;
        int log=0;
        int order=this.theDomain.getNdofs();
        double[] Res = new double[order];
        theSOE = new SOE(order,order);
        theDomain.clearDomain();
        // configure initial conditions for 
        // displacements and velocities
        theDomain.DomainInitialConditions();
        this.formLeft();
        this.formFext(LC, 0, 1.0);
        this.formFint();
        
//        theSOE.printK();
//        theSOE.printM();
        theSOE.TransientRunInit();
        this.updateInitialAcceleration();
        
        this.SpecialStartingProcedure();
        if(Analysis2File){        
            this.AnalysisFile.println();
            this.AnalysisFile.println("LoadCase with id: "+LC);
            this.AnalysisFile.println("______________________");
            this.AnalysisFile.println("time step: "+0);
        }
        this.commit2();
        if(Analysis2File){  
            this.AnalysisFile.print("number of plastified elements: ");
            this.AnalysisFile.println(this.theDomain.getNumOfPlastifiedElems());
        }
        double elastic,kinetic,damped,surfEnrg,surfKnrg,surfEPnrg,diss;
        
        if(printinfo){
            System.out.println("Elastic and Kinetic energy evolution");
            elastic=0.5*theDomain.getuKu();
            kinetic=0.5*theDomain.getvMv();
            damped=0.;//this.timeStep*theDomain.getuCv(a0,a1);
            surfEnrg=0.;
            surfKnrg=0.;
            surfEPnrg=0.;
            for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
                Element theElement = it.next();
                if(theElement.getClass()==OneDOF.class || theElement.getClass()==Crack2d.class){
                    surfEnrg+=0.5*theElement.getuKu();
                    surfKnrg+=0.5*theElement.getvMv();
                    surfEPnrg+=0.5*theElement.PlasticNRG();
                    theElement.checkForDamage();
                }
            }
            diss=theDomain.getDissipation();
    //        double[] u=this.theSOE.getX();
    //        double[] fext=this.theSOE.getFext();
    //        double extwork=0.;
    //        for(int j=0;j<fext.length;j++){extwork+=u[j]*fext[j];}
            System.out.println(0+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+surfEPnrg+" "+damped);
    //        System.out.println(0+" "+diss);
        }else{System.out.println("time step: "+0);}
        // print for trackers
        for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
            Tracker theTracker = it.next();
            theTracker.write(0.);
        }
        
        // TimeSteps
        this.formRight();
        for(int i=1; i<=numSteps; i++){
            this.step=i;
            this.ClearAccumulatedDips();
            this.formFext(LC,i, 1.0);
            this.addM_u(w[2]);
            this.addM_upre(-w[0]);            
            this.addK_u(-1.0);
            
            if(this.RayleighDamping){
                this.addC_upre(w[1]);
            }

            if(Analysis2File){  
                this.AnalysisFile.println();
                this.AnalysisFile.println("time step: "+i);
            }
            log=this.theAlgorithm.solve();
            if(printinfo){
                elastic=0.5*theDomain.getuKu();
                kinetic=0.5*theDomain.getvMv();
                damped=this.timeStep*theDomain.getuCv(a0,a1);

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
                        surfEnrg+=0.5*theElement.getuKu();
                        surfKnrg+=0.5*theElement.getvMv();
                        surfEPnrg+=0.5*theElement.PlasticNRG();
                        bool=false;
                        if(theElement.getDamage())bool=theElement.checkForDamage();
                        if(bool){
                            ids[counter]=theElement.getID(); counter+=1;
                        }
                    }
                }
                diss=theDomain.getDissipation();
    //            u=this.theSOE.getX();
    //            fext=this.theSOE.getFext();
    //            extwork=0.;
    //            for(int j=0;j<fext.length;j++){extwork+=u[j]*fext[j];}
                System.out.print(i+" "+elastic+" "+kinetic+" "+surfEnrg+" "+surfKnrg+" "+diss+" "+surfEPnrg+" "+damped+" ");
    //            System.out.print(i+" "+diss);
                for(int j=0; j<ids.length;j++){
                    if(ids[j]!=0)System.out.print(ids[j]+" ");
                }
                System.out.println();
                for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
                    Element theElement = it.next();
                    if(theElement.getDamage())theElement.commitDamage();
                }
            }else{System.out.println("time step: "+i);}
            if(Analysis2File){  
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
                if(i % theTracker.getStep() == 0)theTracker.write(i*timeStep);
            }
//            System.out.println(i+" "+0.5*theDomain.getuKu()+" "+0.5*theDomain.getvMv());
        }
        if(Analysis2File)this.AnalysisFile.println();
        
        return log;
    }
    
    @Override
    public void formRight() {
        this.formFintEff();
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
            mat=mat.times(w[2]);
            theSOE.addToFint(mat, aInt, 1.0);
            mat=elem.getM_upre();
            mat=mat.times(-w[0]);
            theSOE.addToFint(mat, aInt, 1.0);
            mat=elem.getK_u();
            mat=mat.times(-1.0);
            theSOE.addToFint(mat, aInt, 1.0);
            if(this.RayleighDamping){
                mat=elem.getM_upre();
                mat=mat.times(w[1]*a0);
                theSOE.addToFint(mat, aInt, 1.0);
                mat=elem.getK_upre();
                mat=mat.times(w[1]*a1);
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
                mat=mat.times(w[2]);
                theSOE.addToFint(mat, aInt, 1.0);
                mat=elem.getM_upre();
                mat=mat.times(-w[0]);
                theSOE.addToFint(mat, aInt, 1.0);
                mat=elem.getK_u();
                mat=mat.times(-1.0);
                theSOE.addToFint(mat, aInt, 1.0);
                if(this.RayleighDamping){
                    mat=elem.getM_upre();
                    mat=mat.times(w[1]*a0);
                    theSOE.addToFint(mat, aInt, 1.0);
                    mat=elem.getK_upre();
                    mat=mat.times(w[1]*a1);
                    theSOE.addToFint(mat, aInt, 1.0);
                }
            }
        }
    }

    @Override
    public void run() {
        if(this.RayleighDamping){
            this.theSOE.TransientRun(0.0, w[0], w[3], a0, a1 );
        }else{
            this.theSOE.TransientRun(0.0, w[0]);
        }
    }

    @Override
    public void commit() {
        this.formDisplacement();
//        this.updateVelc();
//        this.updateAccl();
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitDisps();
            theNode.commitLoadCaseDisp(this.currentLC, step);
            if(Analysis2File)theNode.printDisps(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
//        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
//            Node theNode = it.next();
//            theNode.commitVelcs();
//            theNode.printVelcs(this);
//        }
//        this.AnalysisFile.println();
//        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
//            Node theNode = it.next();
//            theNode.commitAccls();
//            theNode.printAccls(this);
//        }
        
        // commit sta elements
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.commit();
        }
    }

    @Override
    public void printGenerals() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    private void SpecialStartingProcedure(){
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.setDispsStepFromDerivatives(timeStep,-timeStep*timeStep/2.0);
        }
    }
    
}
