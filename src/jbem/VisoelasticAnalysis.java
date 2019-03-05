/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

import java.util.Iterator;

/**
 *
 * @author pchr
 */
public class VisoelasticAnalysis extends StaticAnalysis{
    double tau=1.0;
    // constructor
    public VisoelasticAnalysis(Domain theDomain, int steps){
        this.putDomain(theDomain);
        this.steps=steps+1;
    }

    public VisoelasticAnalysis(int steps){
        this.steps=steps+1;
    }
    
//    @Override
//    public void init() {
//        super.init();
//        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
//            Domain theDomain = dit.next();
//            theDomain.setAuxiliaryVariables();
//        }
//    }
    
    @Override
    public void run() {
        for(int itime=0; itime<steps; itime++){
//            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            System.out.println("time step= "+itime+" from total "+steps);
//            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            run(itime,0);
            this.constructLHS=false;
        }
    }
    
    @Override
    public void run(int WhichStep, int Where2Save) {
        if(this.constructLHS){
            formLeft(0);
//            theSOE.getA().print(40, 40);
        }
        formRight(WhichStep);
        
//        theSOE.getB().print(20, 16);
        theSOE.solve();
        //theSOE.getX().print(20, 16);
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();
            if(Where2Save>=0){
                //theSOE.getA().print(20, 16);
                theDomain.updateNodes(theSOE.getX(),WhichStep,Where2Save);
                theDomain.updateRigidMultipliers(theSOE.getX(),WhichStep);
            }else{
                theDomain.updateNodesT(theSOE.getX());
            }
            theDomain.Aux2MainVariables_(WhichStep, tau);
        }
    }
    
}
