/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
import java.util.Iterator;

/**
 *
 * @author pchr
 */
public class SpaceQuadIntegrator extends SpaceIntegrator{
    
    // constructor
    public SpaceQuadIntegrator(){}
    
    public SpaceQuadIntegrator(int numofGP){
        this.setNumofGauss(numofGP);
    }
    
     @Override
    AbstractMatrix IntegrateDu(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        return IntegrateDu_Dv(colNode, theDomain, elem_, elemNode, 0);
    }

    @Override
    AbstractMatrix IntegrateDp(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        AbstractMatrix H;
        EQuad elem=(EQuad) elem_;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        boolean colonelem=false;
        if(elem.getHierOfNode(colNode.getID())!=0){colonelem=true;}
        boolean coidentity=false;
        if(colNode.getID()==elemNode.getID()){coidentity=true;}
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int xmesh=this.getXsi_mesh();
        int ymesh=this.getEta_mesh();
        Transformation_002 theTrans2;
        double transJac2;
        int[] subindex2 = new int[4];
        theTrans2 = new Transformation_002(subindex2);
        int m=FundamentalSolution.theFSdata.getMstep();
        int nn=FundamentalSolution.theFSdata.getNstep();
        
        if(!colonelem){
            for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                    theTrans2.setNodes(subindex2);
                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){

                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                            //xsi=CoordsOnElem[0];
                            //eta=CoordsOnElem[1];
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);
                            
                            R=elem.getDistance(colNode, xsi, eta);
                            Normal=elem.getNormal(xsi, eta);
                            FundamentalSolution.theFSdata.setR(R);
                            FundamentalSolution.theFSdata.setOutwardNormal(Normal);


                            val=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;
                            H=H.plus(
                                    theDomain.getFundamentalSolution().get_p_fund().times(val)
                                    );
                        }
                    }                    
                    
                }
            }
        }else if((!coidentity)&&colonelem){
            int colhierInelem=elem.getHierOfNode(colNode.getID());
            int subs=0;
            switch(colhierInelem){
                case 1: subs=2; break;
                case 2: subs=2; break;
                case 3: subs=2; break;
                case 4: subs=2; break;
                
                case 5: subs=3; break;
                case 6: subs=3; break;
                case 7: subs=3; break;
                case 8: subs=3; break;
                
                case 9: subs=4; break;
            }
            Transformation_001 theTrans;
            double transJac;
            int[] subindex = new int[3];
            for(int ii=1; ii<=subs; ii++){
                if(ii==1){
                    switch(colhierInelem){
                        case 1: subindex[0]=2; subindex[1]=3; subindex[2]=1; break;
                        case 2: subindex[0]=3; subindex[1]=4; subindex[2]=2; break;
                        case 3: subindex[0]=1; subindex[1]=2; subindex[2]=3; break;
                        case 4: subindex[0]=1; subindex[1]=2; subindex[2]=4; break;
                        
                        case 5: subindex[0]=4; subindex[1]=1; subindex[2]=5; break;
                        case 6: subindex[0]=1; subindex[1]=2; subindex[2]=6; break;
                        case 7: subindex[0]=4; subindex[1]=1; subindex[2]=7; break;
                        case 8: subindex[0]=1; subindex[1]=2; subindex[2]=8; break;
                        
                        case 9: subindex[0]=1; subindex[1]=2; subindex[2]=9; break;
                    }
                }else if(ii==2){
                    switch(colhierInelem){
                        case 1: subindex[0]=3; subindex[1]=4; subindex[2]=1; break;
                        case 2: subindex[0]=4; subindex[1]=1; subindex[2]=2; break;
                        case 3: subindex[0]=4; subindex[1]=1; subindex[2]=3; break;
                        case 4: subindex[0]=2; subindex[1]=3; subindex[2]=4; break;
                        
                        case 5: subindex[0]=2; subindex[1]=3; subindex[2]=5; break;
                        case 6: subindex[0]=3; subindex[1]=4; subindex[2]=6; break;
                        case 7: subindex[0]=2; subindex[1]=3; subindex[2]=7; break;
                        case 8: subindex[0]=3; subindex[1]=4; subindex[2]=8; break;
                        
                        case 9: subindex[0]=2; subindex[1]=3; subindex[2]=9; break;
                    }                    
                }else if(ii==3){
                    switch(colhierInelem){                        
                        case 5: subindex[0]=3; subindex[1]=4; subindex[2]=5; break;
                        case 6: subindex[0]=4; subindex[1]=1; subindex[2]=6; break;
                        case 7: subindex[0]=1; subindex[1]=2; subindex[2]=7; break;
                        case 8: subindex[0]=2; subindex[1]=3; subindex[2]=8; break;
                        
                        case 9: subindex[0]=3; subindex[1]=4; subindex[2]=9; break;
                    }
                }else if(ii==4){
                    switch(colhierInelem){
                        case 9: subindex[0]=4; subindex[1]=1; subindex[2]=9; break;
                    }
                }
                theTrans = new Transformation_001(subindex);
                
                for(int xi=1; xi<=xmesh; xi++){
                    for(int yi=1; yi<=ymesh; yi++){
                        subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                        theTrans2.setNodes(subindex2);
                        for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                            for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                                double val=1.0;
                                double xsi_,eta_;
                                CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                                CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                                weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                                weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                                xsi_=theTrans.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                                eta_=theTrans.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                                transJac=theTrans.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                                //xsi=CoordsOnElem[0];
                                //eta=CoordsOnElem[1];

                                xsi=theTrans2.getXsi(xsi_, eta_);
                                eta=theTrans2.getEta(xsi_, eta_);
                                transJac2=theTrans2.getJacobian(xsi_, eta_);

                                R=elem.getDistance(colNode, xsi, eta);
                                Normal=elem.getNormal(xsi, eta);
                                FundamentalSolution.theFSdata.setR(R);
                                FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                                val*=elem.getShapeFunction(hierInelem, xsi, eta);
                                val*=elem.getJacobian(xsi, eta);
                                val*=weight[0]*weight[1];
                                val*=transJac;
                                val*=transJac2;
                                H=H.plus(
                                        theDomain.getFundamentalSolution().get_p_fund().times(val)
                                        );
                            }
                        }
                    }
                }           
            }            
        }else if( coidentity && ((m-nn)>1) ){
            int colhierInelem=elem.getHierOfNode(colNode.getID());
            int subs=0;
            switch(colhierInelem){
                case 1: subs=2; break;
                case 2: subs=2; break;
                case 3: subs=2; break;
                case 4: subs=2; break;
                
                case 5: subs=3; break;
                case 6: subs=3; break;
                case 7: subs=3; break;
                case 8: subs=3; break;
                
                case 9: subs=4; break;
            }
            Transformation_001 theTrans;
            double transJac;
            int[] subindex = new int[3];
            for(int ii=1; ii<=subs; ii++){
                if(ii==1){
                    switch(colhierInelem){
                        case 1: subindex[0]=2; subindex[1]=3; subindex[2]=1; break;
                        case 2: subindex[0]=3; subindex[1]=4; subindex[2]=2; break;
                        case 3: subindex[0]=1; subindex[1]=2; subindex[2]=3; break;
                        case 4: subindex[0]=1; subindex[1]=2; subindex[2]=4; break;
                        
                        case 5: subindex[0]=4; subindex[1]=1; subindex[2]=5; break;
                        case 6: subindex[0]=1; subindex[1]=2; subindex[2]=6; break;
                        case 7: subindex[0]=4; subindex[1]=1; subindex[2]=7; break;
                        case 8: subindex[0]=1; subindex[1]=2; subindex[2]=8; break;
                        
                        case 9: subindex[0]=1; subindex[1]=2; subindex[2]=9; break;
                    }
                }else if(ii==2){
                    switch(colhierInelem){
                        case 1: subindex[0]=3; subindex[1]=4; subindex[2]=1; break;
                        case 2: subindex[0]=4; subindex[1]=1; subindex[2]=2; break;
                        case 3: subindex[0]=4; subindex[1]=1; subindex[2]=3; break;
                        case 4: subindex[0]=2; subindex[1]=3; subindex[2]=4; break;
                        
                        case 5: subindex[0]=2; subindex[1]=3; subindex[2]=5; break;
                        case 6: subindex[0]=3; subindex[1]=4; subindex[2]=6; break;
                        case 7: subindex[0]=2; subindex[1]=3; subindex[2]=7; break;
                        case 8: subindex[0]=3; subindex[1]=4; subindex[2]=8; break;
                        
                        case 9: subindex[0]=2; subindex[1]=3; subindex[2]=9; break;
                    }                    
                }else if(ii==3){
                    switch(colhierInelem){                        
                        case 5: subindex[0]=3; subindex[1]=4; subindex[2]=5; break;
                        case 6: subindex[0]=4; subindex[1]=1; subindex[2]=6; break;
                        case 7: subindex[0]=1; subindex[1]=2; subindex[2]=7; break;
                        case 8: subindex[0]=2; subindex[1]=3; subindex[2]=8; break;
                        
                        case 9: subindex[0]=3; subindex[1]=4; subindex[2]=9; break;
                    }
                }else if(ii==4){
                    switch(colhierInelem){
                        case 9: subindex[0]=4; subindex[1]=1; subindex[2]=9; break;
                    }
                }
                theTrans = new Transformation_001(subindex);
                
                for(int xi=1; xi<=xmesh; xi++){
                    for(int yi=1; yi<=ymesh; yi++){
                        subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                        theTrans2.setNodes(subindex2);
                        for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                            for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                                double val=1.0;
                                double xsi_,eta_;
                                CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                                CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                                weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                                weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                                xsi_=theTrans.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                                eta_=theTrans.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                                transJac=theTrans.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                                //xsi=CoordsOnElem[0];
                                //eta=CoordsOnElem[1];

                                xsi=theTrans2.getXsi(xsi_, eta_);
                                eta=theTrans2.getEta(xsi_, eta_);
                                transJac2=theTrans2.getJacobian(xsi_, eta_);

                                R=elem.getDistance(colNode, xsi, eta);
                                Normal=elem.getNormal(xsi, eta);
                                FundamentalSolution.theFSdata.setR(R);
                                FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                                val*=elem.getShapeFunction(hierInelem, xsi, eta);
                                val*=elem.getJacobian(xsi, eta);
                                val*=weight[0]*weight[1];
                                val*=transJac;
                                val*=transJac2;
                                H=H.plus(
                                        theDomain.getFundamentalSolution().get_p_fund().times(val)
                                        );
                            }
                        }
                    }
                }       
            }
        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDv(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        return IntegrateDu_Dv(colNode, theDomain, elem_, elemNode, 1);
    }
    
    private AbstractMatrix IntegrateDu_Dv(Node colNode, Domain theDomain, Element elem_, Node elemNode,int flag) {
        // flag 0 for u fundamental solution
        // flag 1 for v fundamental solution
        AbstractMatrix H;
        EQuad elem=(EQuad) elem_;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        boolean colonelem=false;
        if(elem.getHierOfNode(colNode.getID())!=0){colonelem=true;}
        boolean coidentity=false;
        if(colNode.getID()==elemNode.getID()){coidentity=true;}
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int xmesh=this.getXsi_mesh();
        int ymesh=this.getEta_mesh();
        Transformation_002 theTrans2;
        double transJac2;
        int[] subindex2 = new int[4];
        theTrans2 = new Transformation_002(subindex2);
        if(!coidentity){
            for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                    theTrans2.setNodes(subindex2);
                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){

                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                            //xsi=CoordsOnElem[0];
                            //eta=CoordsOnElem[1];
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);
                            
                            R=elem.getDistance(colNode, xsi, eta);
                            Normal=elem.getNormal(xsi, eta);
                            FundamentalSolution.theFSdata.setR(R);
                            FundamentalSolution.theFSdata.setOutwardNormal(Normal);


                            val*=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;

                            if(flag==0){
                                H=H.plus(
                                    theDomain.getFundamentalSolution().get_u_fund().times(val)
                                    );
                            } else {
                                TimeFundamentalSolution tFunf;
                                tFunf=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
                                H=H.plus(
                                    tFunf.get_v_fund().times(val)
                                    );
                            }
                        }
                    }                    
                    
                }
            }
        }else{
            int colhierInelem=elem.getHierOfNode(colNode.getID());
            int subs=0;
            switch(colhierInelem){
                case 1: subs=2; break;
                case 2: subs=2; break;
                case 3: subs=2; break;
                case 4: subs=2; break;
                
                case 5: subs=3; break;
                case 6: subs=3; break;
                case 7: subs=3; break;
                case 8: subs=3; break;
                
                case 9: subs=4; break;
            }
            Transformation_001 theTrans;
            double transJac;
            int[] subindex = new int[3];
            for(int ii=1; ii<=subs; ii++){
                if(ii==1){
                    switch(colhierInelem){
                        case 1: subindex[0]=2; subindex[1]=3; subindex[2]=1; break;
                        case 2: subindex[0]=3; subindex[1]=4; subindex[2]=2; break;
                        case 3: subindex[0]=1; subindex[1]=2; subindex[2]=3; break;
                        case 4: subindex[0]=1; subindex[1]=2; subindex[2]=4; break;
                        
                        case 5: subindex[0]=4; subindex[1]=1; subindex[2]=5; break;
                        case 6: subindex[0]=1; subindex[1]=2; subindex[2]=6; break;
                        case 7: subindex[0]=4; subindex[1]=1; subindex[2]=7; break;
                        case 8: subindex[0]=1; subindex[1]=2; subindex[2]=8; break;
                        
                        case 9: subindex[0]=1; subindex[1]=2; subindex[2]=9; break;
                    }
                }else if(ii==2){
                    switch(colhierInelem){
                        case 1: subindex[0]=3; subindex[1]=4; subindex[2]=1; break;
                        case 2: subindex[0]=4; subindex[1]=1; subindex[2]=2; break;
                        case 3: subindex[0]=4; subindex[1]=1; subindex[2]=3; break;
                        case 4: subindex[0]=2; subindex[1]=3; subindex[2]=4; break;
                        
                        case 5: subindex[0]=2; subindex[1]=3; subindex[2]=5; break;
                        case 6: subindex[0]=3; subindex[1]=4; subindex[2]=6; break;
                        case 7: subindex[0]=2; subindex[1]=3; subindex[2]=7; break;
                        case 8: subindex[0]=3; subindex[1]=4; subindex[2]=8; break;
                        
                        case 9: subindex[0]=2; subindex[1]=3; subindex[2]=9; break;
                    }                    
                }else if(ii==3){
                    switch(colhierInelem){
                        case 5: subindex[0]=3; subindex[1]=4; subindex[2]=5; break;
                        case 6: subindex[0]=4; subindex[1]=1; subindex[2]=6; break;
                        case 7: subindex[0]=1; subindex[1]=2; subindex[2]=7; break;
                        case 8: subindex[0]=2; subindex[1]=3; subindex[2]=8; break;
                        
                        case 9: subindex[0]=3; subindex[1]=4; subindex[2]=9; break;
                    }
                }else if(ii==4){
                    switch(colhierInelem){
                        case 9: subindex[0]=4; subindex[1]=1; subindex[2]=9; break;
                    }
                }
                theTrans = new Transformation_001(subindex);
                for(int xi=1; xi<=xmesh; xi++){
                    for(int yi=1; yi<=ymesh; yi++){
                        subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                        theTrans2.setNodes(subindex2);
                        for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                            for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                                double val=1.0;
                                double xsi_,eta_;
                                CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                                CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                                weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                                weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                                xsi_=theTrans.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                                eta_=theTrans.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                                transJac=theTrans.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                                //xsi=CoordsOnElem[0];
                                //eta=CoordsOnElem[1];

                                xsi=theTrans2.getXsi(xsi_, eta_);
                                eta=theTrans2.getEta(xsi_, eta_);
                                transJac2=theTrans2.getJacobian(xsi_, eta_);

                                R=elem.getDistance(colNode, xsi, eta);
                                Normal=elem.getNormal(xsi, eta);
                                FundamentalSolution.theFSdata.setR(R);
                                FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                                val*=elem.getShapeFunction(hierInelem, xsi, eta);
                                val*=elem.getJacobian(xsi, eta);
                                val*=weight[0]*weight[1];
                                val*=transJac;
                                val*=transJac2;

                                if(flag==0){
                                    H=H.plus(
                                    theDomain.getFundamentalSolution().get_u_fund().times(val)
                                    );
                                } else {
                                    TimeFundamentalSolution tFunf;
                                    tFunf=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
                                    H=H.plus(
                                    tFunf.get_v_fund().times(val)
                                    );
                                }
                            }
                        }
                    }
                }             
            }
        }
        return H;
    }
    
    public double IntegrateArea(Element elem_){
        EQuad elem=(EQuad) elem_;
        double Area=0.;
        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            double[] CoordsOnElem=new double[2];
            double[] weight = new double[2];
            double xsi,eta;
            int hierInelem=elem.getHierOfNode(theNode.getID());
            
            int xmesh=this.getXsi_mesh();
            int ymesh=this.getEta_mesh();
            Transformation_002 theTrans2;
            double transJac2;
            int[] subindex2 = new int[4];
            theTrans2 = new Transformation_002(subindex2);
            
            int subs=0;
            switch(hierInelem){
                case 1: subs=2; break;
                case 2: subs=2; break;
                case 3: subs=2; break;
                case 4: subs=2; break;
                
                case 5: subs=3; break;
                case 6: subs=3; break;
                case 7: subs=3; break;
                case 8: subs=3; break;
                
                case 9: subs=4; break;
            }
            Transformation_001 theTrans;
            double transJac;
            
            int[] subindex = new int[3];
            for(int ii=1; ii<=subs; ii++){
                if(ii==1){
                    switch(hierInelem){
                        case 1: subindex[0]=2; subindex[1]=3; subindex[2]=1; break;
                        case 2: subindex[0]=3; subindex[1]=4; subindex[2]=2; break;
                        case 3: subindex[0]=1; subindex[1]=2; subindex[2]=3; break;
                        case 4: subindex[0]=1; subindex[1]=2; subindex[2]=4; break;
                        
                        case 5: subindex[0]=4; subindex[1]=1; subindex[2]=5; break;
                        case 6: subindex[0]=1; subindex[1]=2; subindex[2]=6; break;
                        case 7: subindex[0]=4; subindex[1]=1; subindex[2]=7; break;
                        case 8: subindex[0]=1; subindex[1]=2; subindex[2]=8; break;
                        
                        case 9: subindex[0]=1; subindex[1]=2; subindex[2]=9; break;
                    }
                }else if(ii==2){
                    switch(hierInelem){
                        case 1: subindex[0]=3; subindex[1]=4; subindex[2]=1; break;
                        case 2: subindex[0]=4; subindex[1]=1; subindex[2]=2; break;
                        case 3: subindex[0]=4; subindex[1]=1; subindex[2]=3; break;
                        case 4: subindex[0]=2; subindex[1]=3; subindex[2]=4; break;
                        
                        case 5: subindex[0]=2; subindex[1]=3; subindex[2]=5; break;
                        case 6: subindex[0]=3; subindex[1]=4; subindex[2]=6; break;
                        case 7: subindex[0]=2; subindex[1]=3; subindex[2]=7; break;
                        case 8: subindex[0]=3; subindex[1]=4; subindex[2]=8; break;
                        
                        case 9: subindex[0]=2; subindex[1]=3; subindex[2]=9; break;
                    }                    
                }else if(ii==3){
                    switch(hierInelem){                        
                        case 5: subindex[0]=3; subindex[1]=4; subindex[2]=5; break;
                        case 6: subindex[0]=4; subindex[1]=1; subindex[2]=6; break;
                        case 7: subindex[0]=1; subindex[1]=2; subindex[2]=7; break;
                        case 8: subindex[0]=2; subindex[1]=3; subindex[2]=8; break;
                        
                        case 9: subindex[0]=3; subindex[1]=4; subindex[2]=9; break;
                    }
                }else if(ii==4){
                    switch(hierInelem){
                        case 9: subindex[0]=4; subindex[1]=1; subindex[2]=9; break;
                    }
                }
                theTrans = new Transformation_001(subindex);
                double Area1=0.;
                for(int xi=1; xi<=xmesh; xi++){
                    for(int yi=1; yi<=ymesh; yi++){
                        subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                        theTrans2.setNodes(subindex2);
                        for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                            for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                                double val=1.0;
                                double xsi_,eta_;
                                CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                                CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                                weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                                weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                                xsi_=theTrans.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                                eta_=theTrans.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                                transJac=theTrans.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                                //xsi=CoordsOnElem[0];
                                //eta=CoordsOnElem[1];

                                xsi=theTrans2.getXsi(xsi_, eta_);
                                eta=theTrans2.getEta(xsi_, eta_);
                                transJac2=theTrans2.getJacobian(xsi_, eta_);

                                val*=elem.getShapeFunction(hierInelem, xsi, eta);
                                val*=elem.getJacobian(xsi, eta);
                                val*=weight[0]*weight[1];
                                val*=transJac;
                                val*=transJac2;

                                Area1+=val;
                            }
                        }
                    }
                }
                Area+=Area1;
                //System.out.print(ii);System.out.print(" ");
                //System.out.print(Area1);System.out.print(" ");
                //System.out.println(Area);
            } 
            
            /*
            int xmesh=this.getXsi_mesh();
            int ymesh=this.getEta_mesh();
            Transformation_002 theTrans2;
            double transJac2;
            int[] subindex = new int[4];
            theTrans2 = new Transformation_002(subindex);
            double Area1=0.;
            for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex[0]=xi; subindex[1]=yi; subindex[2]=xmesh; subindex[3]=ymesh;
                    theTrans2.setNodes(subindex);
                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){

                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                            //xsi=CoordsOnElem[0];
                            //eta=CoordsOnElem[1];
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                            val*=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;

                            Area1+=val;
                        }
                    }                    
                    
                }
            }
            Area+=Area1;
            //System.out.print(Area1);System.out.print(" ");System.out.println(Area);
            */
        }
        return Area;
    }
    
    public double IntegrateAreaDeformed(Element elem_, int step, int state){
        // used: http://en.wikipedia.org/wiki/Brahmagupta%27s_formula
        EQuad elem=(EQuad) elem_;
        int numNodes=elem.getNumNodes();
        double Area=0.;
        if(numNodes!=4){
            System.err.println("IntegrateAreaDeformed only for four node element, now the area is given zero value");
        }else{
            double a,b,c,d,s;
            a=0.;
            for(int i=0;i<=2;i++){
                int n1=1; int n2=2;
                a+=(elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state])*
                        (elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state]);
            }
            a=Math.sqrt(a);
            b=0.;
            for(int i=0;i<=2;i++){
                int n1=2; int n2=3;
                b+=(elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state])*
                        (elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state]);
            }
            b=Math.sqrt(b);
            c=0.;
            for(int i=0;i<=2;i++){
                int n1=3; int n2=4;
                c+=(elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state])*
                        (elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state]);
            }
            c=Math.sqrt(c);
            d=0.;
            for(int i=0;i<=2;i++){
                int n1=4; int n2=1;
                d+=(elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state])*
                        (elem.getNodeHier(n1).getCoordinates()[i]+elem.getNodeHier(n1).getu()[i][step][state]-elem.getNodeHier(n2).getCoordinates()[i]-elem.getNodeHier(n2).getu()[i][step][state]);
            }
            d=Math.sqrt(d);
            s=(a+b+c+d)/2.;
            Area=Math.sqrt((s-a)*(s-b)*(s-c)*(s-d));
        }
        
        
        return Area;
    }
    
    @Override
    AbstractMatrix IntegrateDp_res(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        AbstractMatrix H;
        EQuad elem=(EQuad) elem_;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        boolean colonelem=false;
        if(elem.getHierOfNode(colNode.getID())!=0){colonelem=true;}
        boolean coidentity=false;
        if(colNode.getID()==elemNode.getID()){coidentity=true;}
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int xmesh=this.getXsi_mesh();
        int ymesh=this.getEta_mesh();
        Transformation_002 theTrans2;
        double transJac2;
        int[] subindex2 = new int[4];
        theTrans2 = new Transformation_002(subindex2);
        
        if(!colonelem){
            for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                    theTrans2.setNodes(subindex2);
                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){

                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                            //xsi=CoordsOnElem[0];
                            //eta=CoordsOnElem[1];
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);
                            
                            R=elem.getDistance(colNode, xsi, eta);
                            Normal=elem.getNormal(xsi, eta);
                            FundamentalSolution.theFSdata.setR(R);
                            FundamentalSolution.theFSdata.setOutwardNormal(Normal);


                            val=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;
                            TimeFundamentalSolution tFund;
                            tFund=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
                            H=H.plus(
                                    tFund.get_pres_fund().times(val)
                                    //tFund.get_pdif_fund().times(val)
                                    );
                        }
                    }                    
                    
                }
            }
        }else if((!coidentity)&&colonelem){
            int colhierInelem=elem.getHierOfNode(colNode.getID());
            int subs=0;
            switch(colhierInelem){
                case 1: subs=2; break;
                case 2: subs=2; break;
                case 3: subs=2; break;
                case 4: subs=2; break;
                
                case 5: subs=3; break;
                case 6: subs=3; break;
                case 7: subs=3; break;
                case 8: subs=3; break;
                
                case 9: subs=4; break;
            }
            Transformation_001 theTrans;
            double transJac;
            int[] subindex = new int[3];
            for(int ii=1; ii<=subs; ii++){
                if(ii==1){
                    switch(colhierInelem){
                        case 1: subindex[0]=2; subindex[1]=3; subindex[2]=1; break;
                        case 2: subindex[0]=3; subindex[1]=4; subindex[2]=2; break;
                        case 3: subindex[0]=1; subindex[1]=2; subindex[2]=3; break;
                        case 4: subindex[0]=1; subindex[1]=2; subindex[2]=4; break;
                        
                        case 5: subindex[0]=4; subindex[1]=1; subindex[2]=5; break;
                        case 6: subindex[0]=1; subindex[1]=2; subindex[2]=6; break;
                        case 7: subindex[0]=4; subindex[1]=1; subindex[2]=7; break;
                        case 8: subindex[0]=1; subindex[1]=2; subindex[2]=8; break;
                        
                        case 9: subindex[0]=1; subindex[1]=2; subindex[2]=9; break;
                    }
                }else if(ii==2){
                    switch(colhierInelem){
                        case 1: subindex[0]=3; subindex[1]=4; subindex[2]=1; break;
                        case 2: subindex[0]=4; subindex[1]=1; subindex[2]=2; break;
                        case 3: subindex[0]=4; subindex[1]=1; subindex[2]=3; break;
                        case 4: subindex[0]=2; subindex[1]=3; subindex[2]=4; break;
                        
                        case 5: subindex[0]=2; subindex[1]=3; subindex[2]=5; break;
                        case 6: subindex[0]=3; subindex[1]=4; subindex[2]=6; break;
                        case 7: subindex[0]=2; subindex[1]=3; subindex[2]=7; break;
                        case 8: subindex[0]=3; subindex[1]=4; subindex[2]=8; break;
                        
                        case 9: subindex[0]=2; subindex[1]=3; subindex[2]=9; break;
                    }                    
                }else if(ii==3){
                    switch(colhierInelem){                        
                        case 5: subindex[0]=3; subindex[1]=4; subindex[2]=5; break;
                        case 6: subindex[0]=4; subindex[1]=1; subindex[2]=6; break;
                        case 7: subindex[0]=1; subindex[1]=2; subindex[2]=7; break;
                        case 8: subindex[0]=2; subindex[1]=3; subindex[2]=8; break;
                        
                        case 9: subindex[0]=3; subindex[1]=4; subindex[2]=9; break;
                    }
                }else if(ii==4){
                    switch(colhierInelem){
                        case 9: subindex[0]=4; subindex[1]=1; subindex[2]=9; break;
                    }
                }
                theTrans = new Transformation_001(subindex);
                
                for(int xi=1; xi<=xmesh; xi++){
                    for(int yi=1; yi<=ymesh; yi++){
                        subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                        theTrans2.setNodes(subindex2);
                        for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                            for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                                double val=1.0;
                                double xsi_,eta_;
                                CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                                CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                                weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                                weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                                xsi_=theTrans.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                                eta_=theTrans.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                                transJac=theTrans.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                                //xsi=CoordsOnElem[0];
                                //eta=CoordsOnElem[1];

                                xsi=theTrans2.getXsi(xsi_, eta_);
                                eta=theTrans2.getEta(xsi_, eta_);
                                transJac2=theTrans2.getJacobian(xsi_, eta_);

                                R=elem.getDistance(colNode, xsi, eta);
                                Normal=elem.getNormal(xsi, eta);
                                FundamentalSolution.theFSdata.setR(R);
                                FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                                val*=elem.getShapeFunction(hierInelem, xsi, eta);
                                val*=elem.getJacobian(xsi, eta);
                                val*=weight[0]*weight[1];
                                val*=transJac;
                                val*=transJac2;
                                TimeFundamentalSolution tFund;
                                tFund=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
                                H=H.plus(
                                        tFund.get_pres_fund().times(val)
                                        //tFund.get_pdif_fund().times(val)
                                        );
                            }
                        }
                    }
                }           
            }           
        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDp_dif(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        AbstractMatrix H;
        EQuad elem=(EQuad) elem_;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        boolean colonelem=false;
        if(elem.getHierOfNode(colNode.getID())!=0){colonelem=true;}
        boolean coidentity=false;
        if(colNode.getID()==elemNode.getID()){coidentity=true;}
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int xmesh=this.getXsi_mesh();
        int ymesh=this.getEta_mesh();
        Transformation_002 theTrans2;
        double transJac2;
        int[] subindex2 = new int[4];
        theTrans2 = new Transformation_002(subindex2);
        if(coidentity){
            for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                    theTrans2.setNodes(subindex2);
                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){

                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                            //xsi=CoordsOnElem[0];
                            //eta=CoordsOnElem[1];
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);
                            
                            R=elem.getDistance(colNode, xsi, eta);
                            Normal=elem.getNormal(xsi, eta);
                            FundamentalSolution.theFSdata.setR(R);
                            FundamentalSolution.theFSdata.setOutwardNormal(Normal);


                            val=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;
                            TimeFundamentalSolution tFund;
                            tFund=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
                            H=H.plus(
                                    //tFund.get_pres_fund().times(val)
                                    tFund.get_pdif_fund().times(val)
                                    );
                        }
                    }                    
                    
                }
            }
        }
        return H;
    }

    @Override
    double IntegrateM(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        double M=0.;
        EQuad elem=(EQuad) elem_;
        int n=theDomain.theFundSol.get_p_DOFs();

        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];

        double xsi,eta;

        for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
            for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){

                double val=1.0;
                CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

                xsi=CoordsOnElem[0];
                eta=CoordsOnElem[1];

                val=elem.getShapeFunction(hierInelem, xsi, eta);
                val*=elem.getJacobian(xsi, eta);
                val*=weight[0]*weight[1];

                M+=val;
            }
        }

        return M;
    }
    
    public double IntegrateInternalPointDisp(EQuad elem, Node elemNode, Domain theDomain, ResultPoint theInternalPoint, int wdisp, int step, int wstate) {
        double disp=0;
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int xmesh=this.getXsi_mesh();
        int ymesh=this.getEta_mesh();
        Transformation_002 theTrans2;
        double transJac2;
        int[] subindex2 = new int[4];
        theTrans2 = new Transformation_002(subindex2);
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                    theTrans2.setNodes(subindex2);
                    
                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);

//                            xsi=CoordsOnElem[0]; eta=CoordsOnElem[1];
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                            R=elem.getDistance(theInternalPoint, xsi, eta);
                            Normal=elem.getNormal(xsi,eta);
                            FundamentalSolution.theFSdata.setR(R);
                            FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                            val*=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;
                            val*=theDomain.getFundamentalSolution().get_u_fund().get(wdisp, 0)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[0][step]+
                                                                                                ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[0])
                                    +theDomain.getFundamentalSolution().get_u_fund().get(wdisp, 1)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[1][step]+
                                                                                                ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[1])
                                    +theDomain.getFundamentalSolution().get_u_fund().get(wdisp, 2)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[2][step]+
                                                                                                ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[2])
                                    -theDomain.getFundamentalSolution().get_p_fund().get(wdisp, 0)*elem.getNodeHier(hierInelem).getu()[0][step][wstate]
                                    -theDomain.getFundamentalSolution().get_p_fund().get(wdisp, 1)*elem.getNodeHier(hierInelem).getu()[1][step][wstate]
                                    -theDomain.getFundamentalSolution().get_p_fund().get(wdisp, 2)*elem.getNodeHier(hierInelem).getu()[2][step][wstate]
                                    ;
                            disp+=val;
                        }
                    }
                }
        }
        
        return disp;
    }
    
    public double IntegrateInternalPointStress(EQuad elem, Node elemNode, Domain theDomain, ResultPoint theInternalPoint, int wstress, int step, int wstate) {
        double disp=0;
        double[] CoordsOnElem=new double[2];
        double[] weight = new double[2];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int xmesh=this.getXsi_mesh();
        int ymesh=this.getEta_mesh();
        Transformation_002 theTrans2;
        double transJac2;
        int[] subindex2 = new int[4];
        theTrans2 = new Transformation_002(subindex2);
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        for(int xi=1; xi<=xmesh; xi++){
                for(int yi=1; yi<=ymesh; yi++){
                    
                    subindex2[0]=xi; subindex2[1]=yi; subindex2[2]=xmesh; subindex2[3]=ymesh;
                    theTrans2.setNodes(subindex2);
                    
                    for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){
                        for(int j=0; j<SpaceQuadIntegrator.numofGauss; j++){
                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            CoordsOnElem[1]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(j);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                            weight[1]=SpaceQuadIntegrator.theGaussData.getGaussWeight(j);
                            
                            xsi=theTrans2.getXsi(CoordsOnElem[0], CoordsOnElem[1]);
                            eta=theTrans2.getEta(CoordsOnElem[0], CoordsOnElem[1]);
                            transJac2=theTrans2.getJacobian(CoordsOnElem[0], CoordsOnElem[1]);

                            R=elem.getDistance(theInternalPoint, xsi, eta);
                            Normal=elem.getNormal(xsi,eta);
                            FundamentalSolution.theFSdata.setR(R);
                            FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                            val*=elem.getShapeFunction(hierInelem, xsi, eta);
                            val*=elem.getJacobian(xsi, eta);
                            val*=weight[0]*weight[1];
                            val*=transJac2;
                            val*=    theDomain.getFundamentalSolution().get_s_fund().get(0, wstress)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[0][step]+
                                                                                                        ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[0])
                                    +theDomain.getFundamentalSolution().get_s_fund().get(1, wstress)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[1][step]+
                                                                                                        ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[1])
                                    +theDomain.getFundamentalSolution().get_s_fund().get(2, wstress)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[2][step]+
                                                                                                        ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[2])
                                    -theDomain.getFundamentalSolution().get_r_fund().get(0, wstress)*elem.getNodeHier(hierInelem).getu()[0][step][wstate]
                                    -theDomain.getFundamentalSolution().get_r_fund().get(1, wstress)*elem.getNodeHier(hierInelem).getu()[1][step][wstate]
                                    -theDomain.getFundamentalSolution().get_r_fund().get(2, wstress)*elem.getNodeHier(hierInelem).getu()[2][step][wstate]
                                    ;
                            disp+=val;
                        }
                    }
                }
        }
        
        return disp;
    }
}
