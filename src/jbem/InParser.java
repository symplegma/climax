/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;



/**
 *
 * @author pchr
 */
public class InParser extends Parser{
    
    // constructor
    public InParser(){
    }
    
    // methods
    
    public void Parse(){
        Domain theDomain = new Domain(1);
        
        ElasticMat ElasticMaterial = new ElasticMat(1);
        ElasticMaterial.setE_v(1000.,0.);
        ElasticMaterial.setDensity(1.);
        theDomain.addMaterial(ElasticMaterial);
        
        double dt=0.01;
        int timeSteps=3;
        
        double a=2.;
        double b=2.;
        double c=2.;
        double[] cord = new double[3];
        
        cord[0]=-a; cord[1]=-b; cord[2]=c;
        Node aNode = new Node(1,cord);
        theDomain.putNode(aNode);
        cord[0]=a; cord[1]=-b; cord[2]=c;
        aNode = new Node(2,cord);
        theDomain.putNode(aNode);
        
        cord[0]=a; cord[1]=b; cord[2]=c;
        aNode = new Node(3,cord);
        theDomain.putNode(aNode);
        cord[0]=-a; cord[1]=b; cord[2]=c;
        aNode = new Node(4,cord);
        theDomain.putNode(aNode);
        
        cord[0]=-a; cord[1]=-b; cord[2]=-c;
        aNode = new Node(5,cord);
        theDomain.putNode(aNode);
        cord[0]=a; cord[1]=-b; cord[2]=-c;
        aNode = new Node(6,cord);
        theDomain.putNode(aNode);
        
        cord[0]=-a; cord[1]=b; cord[2]=-c;
        aNode = new Node(8,cord);
        theDomain.putNode(aNode);
        cord[0]=a; cord[1]=b; cord[2]=-c;
        aNode = new Node(7,cord);
        theDomain.putNode(aNode);
        
        EQuad4 
        quad4Element = new EQuad4(1,theDomain.getNode(1),
                                    theDomain.getNode(2),
                                    theDomain.getNode(3),
                                    theDomain.getNode(4));
        theDomain.putElement(quad4Element);
        
        quad4Element = new EQuad4(2,theDomain.getNode(2),
                                    theDomain.getNode(6),
                                    theDomain.getNode(7),
                                    theDomain.getNode(3));
        theDomain.putElement(quad4Element);
        
        quad4Element = new EQuad4(3,theDomain.getNode(6),
                                    theDomain.getNode(5),
                                    theDomain.getNode(8),
                                    theDomain.getNode(7));
        theDomain.putElement(quad4Element);
        
        quad4Element = new EQuad4(4,theDomain.getNode(5),
                                    theDomain.getNode(1),
                                    theDomain.getNode(4),
                                    theDomain.getNode(8));
        theDomain.putElement(quad4Element);
        
        quad4Element = new EQuad4(5,theDomain.getNode(4),
                                    theDomain.getNode(3),
                                    theDomain.getNode(7),
                                    theDomain.getNode(8));
        theDomain.putElement(quad4Element);
        
        quad4Element = new EQuad4(6,theDomain.getNode(1),
                                    theDomain.getNode(5),
                                    theDomain.getNode(6),
                                    theDomain.getNode(2));
        theDomain.putElement(quad4Element);

        
        SpaceQuadIntegrator SI=new SpaceQuadIntegrator(6);
        //SI.setMesh(5, 5);
        SI.setMesh(7, 7);
        theDomain.setSpaceIntegrator(SI);
        
        ConstraintEquation aConstraintEquation;
        ConstraintTerm aConstraintTerm;
        double[] cvals= new double[timeSteps+1];
        int id=0;
        
        int uORv=1;
        for(int i=1; i<=3 ; i++){
            aConstraintEquation = new ConstraintEquation(++id);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(1), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

            aConstraintEquation = new ConstraintEquation(++id);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(4), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

            aConstraintEquation = new ConstraintEquation(++id);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(5), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

            aConstraintEquation = new ConstraintEquation(++id);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(8), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

        }
        
        for(int i=1; i<=3 ; i++){
            if(i==1){
                    //cvals[0]=0.;
                    for(int ii=0; ii<timeSteps+1; ii++){
                        
                        //if(ii*dt<0.02){cvals[ii]=0.;}
                        //else{cvals[ii]=-Math.cos(60.*(ii*dt-dt))/60.+1./60;}
                        
                        if(ii*dt<0.01){cvals[ii]=0.;}
                        else{cvals[ii]=-Math.sin(60.*(ii*dt-0.01));}
                        
                        //if(ii*dt<0.01){cvals[ii]=0.;}
                        //else{cvals[ii]=Math.cos(60.*(ii*dt-0.01))/60.-1./60.;}
                        
                    }
                }else{
                    for(int ii=0; ii<timeSteps+1; ii++){
                        cvals[ii]=0.;
                    }
                }
            aConstraintEquation = new ConstraintEquation(++id,cvals);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(2), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

            aConstraintEquation = new ConstraintEquation(++id,cvals);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(3), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

            aConstraintEquation = new ConstraintEquation(++id,cvals);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(6), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

            aConstraintEquation = new ConstraintEquation(++id,cvals);
            aConstraintTerm = new ConstraintTerm(theDomain.getNode(7), uORv, i, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            theDomain.putConstraintEquation(aConstraintEquation);

        }
        
        ConstraintTermElement aConstraintTermElement;
        
        for(int i=1; i<=4; i++){
            for(int j=1; j<=3; j++){
                aConstraintEquation = new ConstraintEquation(++id);
                aConstraintTermElement = new 
                        ConstraintTermElement(theDomain.getElement(1),
                        theDomain.getElement(1).getNodeHier(i).getID(), 
                        j, 1.0);
                aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                theDomain.putConstraintEquation(aConstraintEquation);
                
                aConstraintEquation = new ConstraintEquation(++id);
                aConstraintTermElement = new 
                        ConstraintTermElement(theDomain.getElement(3),
                        theDomain.getElement(3).getNodeHier(i).getID(), 
                        j, 1.0);
                aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                theDomain.putConstraintEquation(aConstraintEquation);
                
                aConstraintEquation = new ConstraintEquation(++id);
                aConstraintTermElement = new 
                        ConstraintTermElement(theDomain.getElement(5),
                        theDomain.getElement(5).getNodeHier(i).getID(), 
                        j, 1.0);
                aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                theDomain.putConstraintEquation(aConstraintEquation);
                
                aConstraintEquation = new ConstraintEquation(++id);
                aConstraintTermElement = new 
                        ConstraintTermElement(theDomain.getElement(6),
                        theDomain.getElement(6).getNodeHier(i).getID(), 
                        j, 1.0);
                aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                theDomain.putConstraintEquation(aConstraintEquation);

                /*
                if(j==1){
                    cvals[0]=0.;
                    for(int ii=0; ii<timeSteps+1; ii++){
                        //cvals[ii]=100.*Math.sin(ii*dt*60.);
                        cvals[ii]=100.;
                    }
                }else{
                    for(int ii=0; ii<timeSteps+1; ii++){
                        cvals[ii]=0.;
                    }
                }
                aConstraintEquation = new ConstraintEquation(++id,cvals);
                aConstraintTermElement = new 
                        ConstraintTermElement(theDomain.getElement(2),
                        theDomain.getElement(2).getNodeHier(i).getID(),
                        j, 1.0);
                aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                theDomain.putConstraintEquation(aConstraintEquation);
                */
            }
        }
        
        TransientElasticity3DFS fund = new TransientElasticity3DFS();
        TransientRodFS.theFSdata.setpTimeDist(1);
        TransientRodFS.theFSdata.setuTimeDist(1);
        TransientRodFS.theFSdata.setvTimeDist(1);
        theDomain.setFundamentalSolution(fund);
        theDomain.setDOFs(timeSteps);
        
        //TransientAnalysis1 theAnalysis=new TransientAnalysis1(theDomain,dt,timeSteps,"C:/Documents and Settings/pchr/Τα έγγραφά μου/NetBeansProjects/JBEMapplication");
        TransientAnalysis1 theAnalysis=new TransientAnalysis1(theDomain,dt,timeSteps);
        //TransientAnalysis_v1 theAnalysis=new TransientAnalysis_v1(theDomain,dt,timeSteps);
        theAnalysis.init();
        theAnalysis.run();
        theDomain.printElement(4, timeSteps-1, "el4.dat");
        theDomain.printElement(2, timeSteps-1, "el2.dat");

    }
}
