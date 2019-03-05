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
public class BetaNewmarkAnalysis extends TimeIntegration{
    private double beta;
    private double gamma;
    private double[] w;
    
    // constructor
    public BetaNewmarkAnalysis(){
        this.theType=3;
    }
    
    public BetaNewmarkAnalysis(int typeofBCImposer, double timeStep, double beta, double gamma){
        theDomain = new Domain();
        this.theType=3;
        this.timeStep=timeStep;
        this.beta=beta;
        this.gamma=gamma;
        //System.out.println("Ωcritical= "+Math.sqrt(gamma/2.-beta)/(gamma/2.-beta));
        theAlgorithm = new Linear(this, 0.001);
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            case 2: theImposer = new ClassicLargeMass(); break;
            case 3: theImposer = new GeneralizedLargeMass(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        
    }
    
    public BetaNewmarkAnalysis(Domain theDomain, double timeStep, double beta, double gamma){
        this.theDomain = theDomain;
        this.theType=3;
        this.timeStep=timeStep;
        this.beta=beta;
        this.gamma=gamma;
        theAlgorithm = new Linear(this, 0.001);
        theAlgorithm.setAnalysis(this);
        theImposer = new StaticImposer();
        w = new double[8];
        w[0]=1./(beta*timeStep*timeStep);   w[1]=gamma/(beta*timeStep); w[2]=1.0/(beta*timeStep);
        w[3]=1.0/(2*beta)-1.;               w[4]=gamma/beta-1.;         w[5]=timeStep*(gamma/beta-2.)/2.0;
        w[6]=timeStep*(1.-gamma);           w[7]=timeStep*gamma;
        if(Analysis2File)this.printGenerals();
    }
    
    public BetaNewmarkAnalysis(Domain theDomain, int typeofBCImposer, Algorithm theAlgorithm, 
                               double timeStep, double beta, double gamma){
        this.theDomain = theDomain;
        this.theType=3;
        this.timeStep=timeStep;
        this.beta=beta;
        this.gamma=gamma;
        this.theAlgorithm=theAlgorithm;
        theAlgorithm.setAnalysis(this);
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            case 2: theImposer = new ClassicLargeMass(); break;
            case 3: theImposer = new GeneralizedLargeMass(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        w = new double[8];
        w[0]=1./(beta*timeStep*timeStep);   w[1]=gamma/(beta*timeStep); w[2]=1.0/(beta*timeStep);
        w[3]=1.0/(2*beta)-1.;               w[4]=gamma/beta-1.;         w[5]=timeStep*(gamma/beta-2.)/2.0;
        w[6]=timeStep*(1.-gamma);           w[7]=timeStep*gamma;
        if(Analysis2File)this.printGenerals();
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
        updateInitialAcceleration();
        
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
        }else{
            //System.out.println("time step: "+0);
        }
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
            this.addM_u(w[0]);
            this.addM_v(w[2]);
            this.addM_a(w[3]);
            
            if(this.RayleighDamping){
                this.addC_u(w[1]);
                this.addC_v(w[4]);
                this.addC_a(w[5]);
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
            }else{
                //System.out.println("time step: "+i);
            }
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
    public void run() {
        if(this.RayleighDamping){
            this.theSOE.TransientRun(1.0, w[0], w[1], a0, a1 );
        }else{
            this.theSOE.TransientRun(1.0, w[0]);
        }
    }

    
    
    //===================================================================
    public void updateVelc(){
        this.formVelocity();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateVelcsStep( theSOE.getXu(), w[0]*w[7]);
            theNode.updateVelcsStep( theSOE.getXv(), -w[7]*w[2] );
            theNode.updateVelcsStep( theSOE.getXa(), w[6]-w[7]*w[3] );
        }
    }
    
    public void updateAccl(){
        this.formAcceleration();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateAcclsStep( theSOE.getXu(), w[0]);
            theNode.updateAcclsStep( theSOE.getXv(), -w[2] );
            theNode.updateAcclsStep( theSOE.getXa(), -w[3]-1. );
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
            mat=mat.times(w[0]);
            theSOE.addToFint(mat, aInt, 1.0);
            if(this.RayleighDamping){
                mat=elem.getM_u();
                mat=mat.times(w[1]*a0);
                theSOE.addToFint(mat, aInt, 1.0);
                mat=elem.getK_u();
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
                mat=mat.times(w[0]);
                theSOE.addToFint(mat, aInt, 1.0);
                if(this.RayleighDamping){
                    mat=elem.getM_u();
                    mat=mat.times(w[1]*a0);
                    theSOE.addToFint(mat, aInt, 1.0);
                    mat=elem.getK_u();
                    mat=mat.times(w[1]*a1);
                    theSOE.addToFint(mat, aInt, 1.0);
                }
            }
        }
    }
}
