/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
import java.util.Iterator;
/**
* This class xxxx.
* @author pchr
* @version xx.xx
*/
public class SpaceLineIntegrator extends SpaceIntegrator{
    private int numGP=5; // gauss point used for the computation of energetics

    // constructor
    public SpaceLineIntegrator(){}

    public SpaceLineIntegrator(int numofGP){
        this.setNumofGauss(numofGP);
    }

    @Override
    AbstractMatrix IntegrateDu(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        return IntegrateDu_Dv(colNode, theDomain, elem_, elemNode, 0);
    }

     @Override
    AbstractMatrix IntegrateDv(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        return IntegrateDu_Dv(colNode, theDomain, elem_, elemNode, 1);
    }

    private AbstractMatrix IntegrateDu_Dv(Node colNode, Domain theDomain, Element elem_, Node elemNode,int flag) {
        // flag 0 for u fundamental solution
        // flag 1 for v fundamental solution
        AbstractMatrix H;
        ELine elem=(ELine) elem_;
        LogFunction Fund= (LogFunction) theDomain.getFundamentalSolution();
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        boolean colonelem=false;
        if(elem.getHierOfNode(colNode.getID())!=0){colonelem=true;}
        boolean coidentity=false;
        if(colNode.getID()==elemNode.getID()){coidentity=true;}
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double[] R;
        double[] Normal;
        double xsi,eta;
        if(!coidentity){
            for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){

                    double val=1.0;
                    CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
                    weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

                    xsi=CoordsOnElem[0];

                    R=elem.getDistance(colNode, xsi);
                    Normal=elem.getNormal(xsi);
                    FundamentalSolution.theFSdata.setR(R);
                    FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                    val*=elem.getShapeFunction(hierInelem, xsi);
                    val*=elem.getJacobian(xsi);
                    val*=weight[0];

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
        }else{
            int nnod=elem.getNumNodes();
            int colhierInelem=elem.getHierOfNode(colNode.getID());
            // nonlog part
            for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){

                    double val=1.0;
                    CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
                    weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

                    xsi=CoordsOnElem[0];

                    R=elem.getDistance(colNode, xsi);
                    FundamentalSolution.theFSdata.setR(R);
                    double f1,f2;
                    f1=0.;f2=0.;

                    switch(nnod){
                        case 2:
                            f1=elem.getNodeHier(2).getCoordinates()[0]-elem.getNodeHier(1).getCoordinates()[0];
                            f2=elem.getNodeHier(2).getCoordinates()[1]-elem.getNodeHier(1).getCoordinates()[1];
                            break;
                        case 3:
                            switch(colhierInelem){
                            case 1:
                                f1=-(2.-xsi)*elem.getNodeHier(1).getCoordinates()[0]+xsi*elem.getNodeHier(2).getCoordinates()[0]
                                        +2.*(1.-xsi)*elem.getNodeHier(3).getCoordinates()[0];
                                f2=-(2.-xsi)*elem.getNodeHier(1).getCoordinates()[1]+xsi*elem.getNodeHier(2).getCoordinates()[1]
                                        +2.*(1.-xsi)*elem.getNodeHier(3).getCoordinates()[1];
                                break;
                            case 2:
                                f1=-(2.+xsi)*elem.getNodeHier(2).getCoordinates()[0]-xsi*elem.getNodeHier(1).getCoordinates()[0]
                                        +2.*(1.+xsi)*elem.getNodeHier(3).getCoordinates()[0];
                                f2=-(2.+xsi)*elem.getNodeHier(2).getCoordinates()[1]-xsi*elem.getNodeHier(1).getCoordinates()[1]
                                        +2.*(1.+xsi)*elem.getNodeHier(3).getCoordinates()[1];
                                break;
                            case 3:
                                f1=0.5*((xsi-1.)*elem.getNodeHier(1).getCoordinates()[0]+(xsi+1.)*elem.getNodeHier(2).getCoordinates()[0])
                                        -xsi*elem.getNodeHier(3).getCoordinates()[0];
                                f2=0.5*((xsi-1.)*elem.getNodeHier(1).getCoordinates()[1]+(xsi+1.)*elem.getNodeHier(2).getCoordinates()[1])
                                        -xsi*elem.getNodeHier(3).getCoordinates()[1];
                                break;
                        };
                        break;
                    }

                    val*=elem.getShapeFunction(hierInelem, xsi);
                    val*=elem.getJacobian(xsi);
                    val*=weight[0];

                    if(flag==0){
                        H=H.plus(
                                Fund.getuNonLogPart(-0.5*Math.log(f2*f2+f1*f1)).times(val)
                            );
                    } else {
                        H=H.plus(
                                Fund.getvNonLogPart().times(val)
                            );
                    }
            }
            // end of nonlog part
            // log part
            int subs=0;
            switch(colhierInelem){
                case 1: subs=1; break;
                case 2: subs=1; break;

                case 3: subs=2; break;
            }
            int tempNumGauss=SpaceLineIntegrator.numofGauss;
            this.setNumofGauss(8);
            SpaceLineIntegrator.theGaussData.setGaussType(2);
            SpaceLineIntegrator.theGaussData.computeGaussPoints();
            for(int ii=1; ii<=subs; ii++){
                for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){
                    double transJac=0.;
                    double val=1.0;
                    CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                    weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);
                    double xp=0.;
                    xsi=0.;
                    R=new double[3];
                    //switch(nnod)
                    switch(nnod){
                        case 2:
                            transJac=2.0;
                            switch(colhierInelem){
                                case 1: xp=-1. ;break;
                                case 2: xp= 1. ;break;
                            }
                            xsi=(2.*CoordsOnElem[0]-1.)*(-xp);
                            R[0]=0.5*(1.+xsi*xp); R[1]=0.; R[2]=0.;
                            //R[0]=0.5*(1.-xsi*xp); R[1]=0.; R[2]=0.;
                            break;
                        case 3:
                            switch(colhierInelem){
                                case 3:
                                    transJac=1.;
                                    switch(ii){
                                        case 1: xp=-1. ;break;
                                        case 2: xp= 1. ;break;
                                    }
                                    xsi=Math.abs(CoordsOnElem[0]);
                                    R[0]=xp*xsi; R[1]=0.; R[2]=0.;
                                    break;
                                default:
                                    transJac=2.;
                                    switch(colhierInelem){
                                        case 1: xp=-1. ;break;
                                        case 2: xp= 1. ;break;
                                    }
                                    xsi=(2.*CoordsOnElem[0]-1.)*(-xp);
                                    R[0]=0.5*(1.+xsi*xp); R[1]=0.; R[2]=0.;
                                    //R[0]=0.5*(1.-xsi*xp); R[1]=0.; R[2]=0.;
                                    break;
                            }
                            break;
                    }

                    Normal=elem.getNormal(xsi);
                    FundamentalSolution.theFSdata.setR(R);
                    FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                    val*=elem.getShapeFunction(hierInelem, xsi);
                    val*=elem.getJacobian(xsi);
                    val*=weight[0];
                    val*=transJac;
                    if(flag==0){
                        H=H.plus(
                                Fund.getuLogPart().times(val)
                                );
                    } else {
                        H=H.plus(
                                Fund.getvLogPart().times(val)
                                );
                    }
                }
            }
            // end of log part
            SpaceLineIntegrator.theGaussData.setGaussType(1);
            this.setNumofGauss(tempNumGauss);

        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDp(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        AbstractMatrix H;
        ELine elem=(ELine) elem_;
        int n=theDomain.theFundSol.get_p_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        boolean colonelem=false;
        if(elem.getHierOfNode(colNode.getID())!=0){colonelem=true;}
        boolean coidentity=false;
        if(colNode.getID()==elemNode.getID()){coidentity=true;}
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double[] R;
        double[] Normal;
        double xsi,eta;
        if(!coidentity){
            for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){

                    double val=1.0;
                    CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
                    weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

                    xsi=CoordsOnElem[0];

                    R=elem.getDistance(colNode, xsi);
                    Normal=elem.getNormal(xsi);
                    FundamentalSolution.theFSdata.setR(R);
                    FundamentalSolution.theFSdata.setOutwardNormal(Normal);

                    val*=elem.getShapeFunction(hierInelem, xsi);
                    val*=elem.getJacobian(xsi);
                    val*=weight[0];

                    H=H.plus(
                            theDomain.getFundamentalSolution().get_p_fund().times(val)
                            );
            }
        }
        return H;
    }

    public double IntegrateLength(ELine elem){
        double Length=0.;
        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            double[] CoordsOnElem=new double[1];
            double[] weight = new double[1];
            double xsi;
            int hierInelem=elem.getHierOfNode(theNode.getID());

            for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){

                    double val=1.0;
                    CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                    weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);

                    xsi=CoordsOnElem[0];


                    val*=elem.getShapeFunction(hierInelem, xsi);
                    val*=elem.getJacobian(xsi);
                    val*=weight[0];

                    Length+=val;
            }
        }
        return Length;
    }

    public double IntegrateLength(IELine elem){
        double Length=0.;
        for(Iterator<InterfaceNode>it=elem.getNodes().values().iterator(); it.hasNext();){
            InterfaceNode theNode = it.next();
            double[] CoordsOnElem=new double[1];
            double[] weight = new double[1];
            double xsi;
            int hierInelem=elem.getHierOfINode(theNode.getID());

            for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){

                    double val=1.0;
                    CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                    weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);

                    xsi=CoordsOnElem[0];


                    val*=elem.getShapeFunction(hierInelem, xsi);
                    val*=elem.getJacobian(xsi);
                    val*=weight[0];

                    Length+=val;
            }
        }
        return Length;
    }


    public double IntegrateSideInequality(ELine elem, Domain theDomain, int step, int bound){
        double Length=0.;
        double temp=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;
        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){

                        double val=1.0;
                        CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                        weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);

                        xsi=CoordsOnElem[0];

                        ux=elem.getNodeHier(hierInelem_i).getuh()[0];
                        uy=elem.getNodeHier(hierInelem_i).getuh()[1];

                        px=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0];
                        py=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1];

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        if(bound!=3){
            temp=Length;
            Length=0.;
            int sidestep=step;
            if(bound==1)sidestep=step-1;
            for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
                theNode_i= it.next();
                for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                    theNode_j= jt.next();

                    double[] CoordsOnElem=new double[1];
                    double[] weight = new double[1];
                    double xsi;

                    int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                    int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                    for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){

                            double val=1.0;
                            CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                            weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);

                            xsi=CoordsOnElem[0];

                            ux=elem.getNodeHier(hierInelem_i).getu()[0][sidestep][0];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][sidestep][0];

                            px=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0];
                            py=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1];

                            val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                            val*=elem.getJacobian(xsi);
                            val*=weight[0];

                            Length+=val;
                    }
                }
            }
            if(bound==0){Length=Length-temp/2.;}else{Length=Length+temp/2.;}
        }
        return Length;
    }

    @Override
    AbstractMatrix IntegrateDp_res(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    AbstractMatrix IntegrateDp_dif(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }


    // **************************************************************************
    // NEW METHODS On 27/07/2010 Panagiotopoulos Christos
    // **************************************************************************
    public double IntegrateWork(ELine elem, Domain theDomain, int step){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][0];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][0];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][0];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][0];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }

    public double IntegrateWork_hmg(ELine elem, Domain theDomain){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux=elem.getNodeHier(hierInelem_i).getuh()[0];
                        uy=elem.getNodeHier(hierInelem_i).getuh()[1];

                        px=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0];
                        py=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1];

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }

    public double IntegrateWorkX(ELine elem, Domain theDomain, int step){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][0];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][0];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][0];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][0];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }
                        
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }

    public double IntegrateWorkY(ELine elem, Domain theDomain, int step){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }


        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][0];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][0];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][0];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][0];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }
                        
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }


    public double IntegrateGeneralisedForce(ELine elem, Domain theDomain, int step, int nodeID, int dof){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux=0.;
                        uy=0.;

                        if(theNode_i.getID() == nodeID){
                            switch(dof){
                                case 1: ux=ux=1.; break;
                                case 2: uy=uy=1.; break;
                            }
                        }
                        
                        if(this.DomainAuxiliaryField){
                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }else{
                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];
                        }
                        

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateGeneralisedForce(ELine elem, Domain theDomain, int DispStep, int TracStep, int DispState, int TracState, int nodeID, int dof){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux=0.;
                        uy=0.;

                        if(theNode_i.getID() == nodeID){
                            switch(dof){
                                case 1: ux=1.; break;
                                case 2: uy=1.; break;
                            }
                        }
                        
                        if(this.DomainAuxiliaryField){
                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[0][TracStep];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[1][TracStep];
                        }else{
                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[0][TracStep];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[1][TracStep];
                        }
                        
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }

    public double IntegrateGeneralisedForce_hmg(ELine elem, Domain theDomain, int nodeID, int dof){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux=0.;
                        uy=0.;

                        if(theNode_i.getID() == nodeID){
                            switch(dof){
                                case 1: ux=elem.getNodeHier(hierInelem_i).getuh()[0]; ux=1.; break;
                                case 2: uy=elem.getNodeHier(hierInelem_i).getuh()[1]; uy=1.; break;
                            }
                        }

                        px=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0];
                        py=elem.getNodeHier(hierInelem_j).getph(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1];

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }



    /////////////////////////////////////////////////////////////////////////////////////////
    // BEGIN OF: methods for integrals on interfaces and interface elements
    // 27 Abril 2011
    /////////////////////////////////////////////////////////////////////////////////////////
    /**
     * It is developed for 2-node interface line elements only.
     * Computes after integration on element, the stored energy of normal spring distribution
     * @author pchr
     * @return The value input as a double.
     * @version dd/mm/yy
     */
    public double IntegrateNormalStoredIEnergy(IELine elem, Interface theInterface, int wstep){
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        double kh = 0.0;
        double k0 = 0.0; double r=1.;
        double ge = 0.0; // gradient plasticity

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsNormalIncluded() && theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
            kh = theInterface.getNormalSpring().getIsotropicHardening();
            k0 = theInterface.getNormalSpring().get_GradientZ();
            ge = theInterface.getNormalSpring().get_GradientPlasticity();
            if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();

            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double s1 = 0,s2 = 0;
            double du1,du2;
            double N1,N2;
            double dN1,dN2;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                s1 = 0;s2 = 0;
                du1 = 0;du2 = 0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                if(theInterface.IsPlastIncluded())s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                if(theInterface.IsPlastIncluded())s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                du1=(ua1-ub1)-s1;
                du2=(ua2-ub2)-s2;

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}

                val*=0.5*kn*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                        2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                        (N1*s1+N2*s2)*(N1*s1+N2*s2)
                        )*(N1*z1+N2*z2)
                    +0.5*kh*(N1*s1+N2*s2)*(N1*s1+N2*s2)+k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r+ge*Math.pow(Math.abs(dN1*s1+dN2*s2), 2.);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }
    
    public double IntegrateDamageDrivingForce(IELine elem, Interface theInterface, int wstep){
        // writen in Crete for paper prmdelam-mix
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        double kt = 0.0;
        double N1,N2;
        
        double z1 ,z2 ;
        double zw1 ,zw2 ;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        double s1 ,s2 ;
        double p1 ,p2 ;
        
        int pstep = 0;
        if(wstep>0){
            
            pstep=wstep-1;
        
            if(theInterface.IsTangentialIncluded())kt = theInterface.getTangentSpring().getStiffness();
            if(theInterface.IsNormalIncluded())kn = theInterface.getNormalSpring().getStiffness();

            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zw1 = 1.; zw2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                va1=0.;va2=0.;
                vb1=0.;vb2=0.;
                s1=0.; s2=0.;
                p1=0.; p2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsSlipIncluded()){
                    s1 = elem.getINodeHier(1).gets(elem.getID())[pstep]; 
                    s2 = elem.getINodeHier(2).gets(elem.getID())[pstep];
                }
                
                if(theInterface.IsPlastIncluded()){
                    p1 = elem.getINodeHier(1).getp(elem.getID())[pstep]; 
                    p2 = elem.getINodeHier(2).getp(elem.getID())[pstep];
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[pstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[pstep];

                if(theInterface.IsDamageIncluded())zw1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())zw2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                if(theInterface.IsNormalIncluded())ua1 = elem.getINodeHier(1).getun()[pstep];
                if(theInterface.IsNormalIncluded())ua2 = elem.getINodeHier(2).getun()[pstep];

                if(theInterface.IsTangentialIncluded())va1 = elem.getINodeHier(1).getut()[pstep];
                if(theInterface.IsTangentialIncluded())va2 = elem.getINodeHier(2).getut()[pstep];

                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        if(theInterface.IsNormalIncluded()){
                            ua1 = elem.getINodeHier(1).getun()[pstep];
                            ub1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                            ua2 = elem.getINodeHier(2).getun()[pstep];
                            ub2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                        }
                        if(theInterface.IsTangentialIncluded()){
                            va1 = elem.getINodeHier(1).getut()[pstep];
                            vb1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                            va2 = elem.getINodeHier(2).getut()[pstep];
                            vb2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                        }
                    }else{
                        if(theInterface.IsNormalIncluded()){
                            ub1 = elem.getINodeHier(1).getun()[pstep];
                            ua1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                            ub2 = elem.getINodeHier(2).getun()[pstep];
                            ua2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                        }
                        if(theInterface.IsTangentialIncluded()){
                            vb1 = elem.getINodeHier(1).getut()[pstep];
                            va1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                            vb2 = elem.getINodeHier(2).getut()[pstep];
                            va2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                        }
                    }
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                val*=0.5*(kn*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*p1+N2*p2)+
                        2.*(N1*ub1+N2*ub2)*(N1*p1+N2*p2)+
                        (N1*p1+N2*p2)*(N1*p1+N2*p2)
                        )
                        +
                        kt*(
                        (N1*va1+N2*va2)*(N1*va1+N2*va2)+
                        (N1*vb1+N2*vb2)*(N1*vb1+N2*vb2)-
                        2.*(N1*va1+N2*va2)*(N1*vb1+N2*vb2)-
                        2.*(N1*va1+N2*va2)*(N1*s1+N2*s2)+
                        2.*(N1*vb1+N2*vb2)*(N1*s1+N2*s2)+
                        (N1*s1+N2*s2)*(N1*s1+N2*s2)
                        )
                        )*(N1*(z1-zw1)+N2*(z2-zw2)); //*(N1*z1+N2*z2)*(N1*(z1-zw1)+N2*(z2-zw2));

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
            
        }
        
        return stored;
    }
    
    public double IntegrateSlipDrivingForce(IELine elem, Interface theInterface, int wstep){
        // writen in Crete for paper prmdelam-mix
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kh = 0.0;
        double kt = 0.0;
        double N1,N2;
        
        double z1 ,z2 ;
        double zw1 ,zw2 ;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        double sw1 ,sw2 ;
        
        int pstep = 0;
        if(wstep>0){
            
            pstep=wstep-1;
        
            kt = theInterface.getTangentSpring().getStiffness();
            kh = theInterface.getTangentSpring().getIsotropicHardening();

            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zw1 = 1.; zw2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                sw1=0.; sw2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsSlipIncluded()){
                    s1 = elem.getINodeHier(1).gets(elem.getID())[pstep]; 
                    s2 = elem.getINodeHier(2).gets(elem.getID())[pstep];
                    
                    sw1 = elem.getINodeHier(1).gets(elem.getID())[wstep]; 
                    sw2 = elem.getINodeHier(2).gets(elem.getID())[wstep];
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[pstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[pstep];

                if(theInterface.IsDamageIncluded())zw1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())zw2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1 = elem.getINodeHier(1).getut()[pstep];
                ua2 = elem.getINodeHier(2).getut()[pstep];

                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[pstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                        ua2 = elem.getINodeHier(2).getut()[pstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[pstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                        ub2 = elem.getINodeHier(2).getut()[pstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                    }
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                val*=(kt*((N1*s1+N2*s2)-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2))*(N1*z1+N2*z2)+
                        kh*(N1*s1+N2*s2)
                        )*(N1*(s1-sw1)+N2*(s2-sw2));

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
            
        }
        return stored;
    }
    
    public double IntegratePlastDrivingForce(IELine elem, Interface theInterface, int wstep){
        // writen in Crete for paper prmdelam-mix
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kh = 0.0;
        double kn = 0.0;
        double N1,N2;
        
        double z1 ,z2 ;
        double zw1 ,zw2 ;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        double sw1 ,sw2 ;
        
        int pstep = 0;
        if(wstep>0){
            
            pstep=wstep-1;
        
            kn = theInterface.getNormalSpring().getStiffness();
            kh = theInterface.getNormalSpring().getIsotropicHardening();

            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zw1 = 1.; zw2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                sw1=0.; sw2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsPlastIncluded()){
                    s1 = elem.getINodeHier(1).getp(elem.getID())[pstep]; 
                    s2 = elem.getINodeHier(2).getp(elem.getID())[pstep];
                    
                    sw1 = elem.getINodeHier(1).getp(elem.getID())[wstep]; 
                    sw2 = elem.getINodeHier(2).getp(elem.getID())[wstep];
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[pstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[pstep];

                if(theInterface.IsDamageIncluded())zw1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())zw2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1 = elem.getINodeHier(1).getun()[pstep];
                ua2 = elem.getINodeHier(2).getun()[pstep];

                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[pstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                        ua2 = elem.getINodeHier(2).getun()[pstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[pstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                        ub2 = elem.getINodeHier(2).getun()[pstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                    }
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                val*=(kn*((N1*s1+N2*s2)-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2))*(N1*z1+N2*z2)+
                        kh*(N1*s1+N2*s2)
                        )*(N1*(s1-sw1)+N2*(s2-sw2));

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
            
        }
        return stored;
    }
    
    public double IntegrateNormalStoredIEnergy(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        double kh = 0.0;
        double k0 = 0.0; double r=1.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsNormalIncluded() && theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
            kh = theInterface.getNormalSpring().getIsotropicHardening();
            k0 = theInterface.getNormalSpring().get_GradientZ();
            if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();

            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double N1,N2;
            double dN1,dN2;
            double va1,va2;
            double vb1 = 0,vb2 = 0; 

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }
                }
                pstep=wstep-1;if(wstep==0)pstep=0;
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getun_(vMaterial_main,tau)[pstep];
                va2 = elem.getINodeHier(2).getun_(vMaterial_main,tau)[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}
                double nu = vMaterial_main.getDispRelaxationTime_1();
                val*=0.5*kn*(tau*tau/((tau+nu)*(tau+nu)))*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)
                        )*(N1*z1+N2*z2)
                        +0.5*kn*(nu*nu/((tau+nu)*(tau+nu)))*(
                        (N1*va1+N2*va2)*(N1*va1+N2*va2)+
                        (N1*vb1+N2*vb2)*(N1*vb1+N2*vb2)-
                        2.*(N1*va1+N2*va2)*(N1*vb1+N2*vb2)
                        )*(N1*z1+N2*z2)
                        +kn*(tau*nu/((tau+nu)*(tau+nu)))*(
                        (N1*ua1+N2*ua2-N1*ub1-N2*ub2)*(N1*va1+N2*va2-N1*vb1-N2*vb2))*(N1*z1+N2*z2)
                        +k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r;
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val*((tau+nu)/tau);
            }
        }
        return stored;
    }
    
    public double IntegrateNormalStoredIEnergyQuad(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsNormalIncluded() && theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();

            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double N1,N2;
            double dN1,dN2;
            double va1,va2;
            double vb1 = 0,vb2 = 0; 

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }
                }
                
                pstep=wstep-1; if(wstep==0)pstep=0;
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getun_(vMaterial_main,tau)[pstep];
                va2 = elem.getINodeHier(2).getun_(vMaterial_main,tau)[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                val*=0.5*kn*(tau*tau/((tau+vMaterial_main.getDispRelaxationTime_1())*(tau+vMaterial_main.getDispRelaxationTime_1())))*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)
                        )*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val*((tau+vMaterial_main.getDispRelaxationTime_1())/tau);
            }
        }
        return stored;
    }
    
    public double IntegrateNormalStoredIEnergyLinear(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsNormalIncluded() && theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();

            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double N1,N2;
            double dN1,dN2;
            double va1,va2;
            double vb1 = 0,vb2 = 0; 

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }
                }
                
                pstep=wstep-1; if(wstep==0)pstep=0;
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getun_(vMaterial_main,tau)[pstep];
                va2 = elem.getINodeHier(2).getun_(vMaterial_main,tau)[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}

                val*=kn*(tau*vMaterial_main.getDispRelaxationTime_1()/((tau+vMaterial_main.getDispRelaxationTime_1())*(tau+vMaterial_main.getDispRelaxationTime_1())))*(
                        (N1*ua1+N2*ua2-N1*ub1-N2*ub2)*(N1*va1+N2*va2-N1*vb1-N2*vb2))*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val*((tau+vMaterial_main.getDispRelaxationTime_1())/tau);
            }
        }
        return stored;
    }
    
    public double IntegrateNormalStoredIEnergy_split(IELine elem, Interface theInterface, int wstep){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        double kh = 0.0;
        double k0 = 0.0; double r=1.;
        double ge = 0.0; // gradient plasticity

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsNormalIncluded() && theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
            kh = theInterface.getNormalSpring().getIsotropicHardening();
            k0 = theInterface.getNormalSpring().get_GradientZ();
            ge = theInterface.getNormalSpring().get_GradientPlasticity();
            if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();

            double z1 = 1.,z2 = 1.;
            double ua1pos,ua2pos;
            double ub1pos,ub2pos;
            double ua1neg,ua2neg;
            double ub1neg,ub2neg;
            double s1 = 0,s2 = 0;
            double N1,N2;
            double dN1,dN2;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1pos=0.;ub2pos=0.;
                ub1neg=0.;ub2neg=0.;
                z1 = 1.;z2 = 1.;
                s1 = 0;s2 = 0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1pos = elem.getINodeHier(1).getun_pos()[wstep];
                ua2pos = elem.getINodeHier(2).getun_pos()[wstep];
                ua1neg = elem.getINodeHier(1).getun_neg()[wstep];
                ua2neg = elem.getINodeHier(2).getun_neg()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1pos = elem.getINodeHier(1).getun_pos()[wstep];
                        ub1pos = elem.getINodeHier(1).getTwin().getun_pos()[wstep];
                        ua2pos = elem.getINodeHier(2).getun_pos()[wstep];
                        ub2pos = elem.getINodeHier(2).getTwin().getun_pos()[wstep];
                        
                        ua1neg = elem.getINodeHier(1).getun_neg()[wstep];
                        ub1neg = elem.getINodeHier(1).getTwin().getun_neg()[wstep];
                        ua2neg = elem.getINodeHier(2).getun_neg()[wstep];
                        ub2neg = elem.getINodeHier(2).getTwin().getun_neg()[wstep];
                    }else{
                        ub1pos = elem.getINodeHier(1).getun_pos()[wstep];
                        ua1pos = elem.getINodeHier(1).getTwin().getun_pos()[wstep];
                        ub2pos = elem.getINodeHier(2).getun_pos()[wstep];
                        ua2pos = elem.getINodeHier(2).getTwin().getun_pos()[wstep];
                        
                        ub1neg = elem.getINodeHier(1).getun_neg()[wstep];
                        ua1neg = elem.getINodeHier(1).getTwin().getun_neg()[wstep];
                        ub2neg = elem.getINodeHier(2).getun_neg()[wstep];
                        ua2neg = elem.getINodeHier(2).getTwin().getun_neg()[wstep];
                    }
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                if(theInterface.IsPlastIncluded())s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                if(theInterface.IsPlastIncluded())s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];


                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}

                val*=(kn*(Math.pow(N1,2)*Math.pow(s1,2) + 2*N1*N2*s1*s2 + Math.pow(N2,2)*Math.pow(s2,2) - 2*Math.pow(N1,2)*s1*ua1pos - 2*N1*N2*s2*ua1pos + Math.pow(N1,2)*Math.pow(ua1pos,2) - 2*N1*N2*s1*ua2pos -
       2*Math.pow(N2,2)*s2*ua2pos + 2*N1*N2*ua1pos*ua2pos + Math.pow(N2,2)*Math.pow(ua2pos,2) + 2*Math.pow(N1,2)*s1*ub1pos + 2*N1*N2*s2*ub1pos - 2*Math.pow(N1,2)*ua1pos*ub1pos - 2*N1*N2*ua2pos*ub1pos +
       Math.pow(N1,2)*Math.pow(ub1pos,2) + 2*N1*N2*s1*ub2pos + 2*Math.pow(N2,2)*s2*ub2pos - 2*N1*N2*ua1pos*ub2pos - 2*Math.pow(N2,2)*ua2pos*ub2pos + 2*N1*N2*ub1pos*ub2pos + Math.pow(N2,2)*Math.pow(ub2pos,2) -
       2*Math.pow(N1,2)*s1*ua1neg*(N1*z1 + N2*z2) - 2*N1*N2*s2*ua1neg*(N1*z1 + N2*z2) + Math.pow(N1,2)*Math.pow(ua1neg,2)*(N1*z1 + N2*z2) + 2*Math.pow(N1,2)*ua1neg*ua1pos*(N1*z1 + N2*z2) -
       2*N1*N2*s1*ua2neg*(N1*z1 + N2*z2) - 2*Math.pow(N2,2)*s2*ua2neg*(N1*z1 + N2*z2) + 2*N1*N2*ua1neg*ua2neg*(N1*z1 + N2*z2) + 2*N1*N2*ua1pos*ua2neg*(N1*z1 + N2*z2) +
       Math.pow(N2,2)*Math.pow(ua2neg,2)*(N1*z1 + N2*z2) + 2*N1*N2*ua1neg*ua2pos*(N1*z1 + N2*z2) + 2*Math.pow(N2,2)*ua2neg*ua2pos*(N1*z1 + N2*z2) + 2*Math.pow(N1,2)*s1*ub1neg*(N1*z1 + N2*z2) +
       2*N1*N2*s2*ub1neg*(N1*z1 + N2*z2) - 2*Math.pow(N1,2)*ua1neg*ub1neg*(N1*z1 + N2*z2) - 2*Math.pow(N1,2)*ua1pos*ub1neg*(N1*z1 + N2*z2) - 2*N1*N2*ua2neg*ub1neg*(N1*z1 + N2*z2) -
       2*N1*N2*ua2pos*ub1neg*(N1*z1 + N2*z2) + Math.pow(N1,2)*Math.pow(ub1neg,2)*(N1*z1 + N2*z2) - 2*Math.pow(N1,2)*ua1neg*ub1pos*(N1*z1 + N2*z2) - 2*N1*N2*ua2neg*ub1pos*(N1*z1 + N2*z2) +
       2*Math.pow(N1,2)*ub1neg*ub1pos*(N1*z1 + N2*z2) + 2*N1*N2*s1*ub2neg*(N1*z1 + N2*z2) + 2*Math.pow(N2,2)*s2*ub2neg*(N1*z1 + N2*z2) - 2*N1*N2*ua1neg*ub2neg*(N1*z1 + N2*z2) -
       2*N1*N2*ua1pos*ub2neg*(N1*z1 + N2*z2) - 2*Math.pow(N2,2)*ua2neg*ub2neg*(N1*z1 + N2*z2) - 2*Math.pow(N2,2)*ua2pos*ub2neg*(N1*z1 + N2*z2) + 2*N1*N2*ub1neg*ub2neg*(N1*z1 + N2*z2) +
       2*N1*N2*ub1pos*ub2neg*(N1*z1 + N2*z2) + Math.pow(N2,2)*Math.pow(ub2neg,2)*(N1*z1 + N2*z2) - 2*N1*N2*ua1neg*ub2pos*(N1*z1 + N2*z2) - 2*Math.pow(N2,2)*ua2neg*ub2pos*(N1*z1 + N2*z2) +
       2*N1*N2*ub1neg*ub2pos*(N1*z1 + N2*z2) + 2*Math.pow(N2,2)*ub2neg*ub2pos*(N1*z1 + N2*z2)))/2.
                        +0.5*kh*(N1*s1+N2*s2)*(N1*s1+N2*s2)+k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r+ge*Math.pow(Math.abs(dN1*s1+dN2*s2), 2.);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }

    public double IntegrateTangentialStoredIEnergy(IELine elem, Interface theInterface, int wstep){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kt = 0.0;
        double kt_0 = 0.0;
        double kh = 0.0;
        double k0 = 0.0; double r=1.;
        double ge = 0.0; // gradient plasticity
        double a0 = 0.0;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsTangentialIncluded() && theInterface.getTangentSpring()!=null){
            kt = theInterface.getTangentSpring().getStiffness();
            kt_0 = theInterface.getTangentSpring().getStiffness_0();
            kh = theInterface.getTangentSpring().getIsotropicHardening();
            k0 = theInterface.getTangentSpring().get_GradientZ();
            ge = theInterface.getTangentSpring().get_GradientPlasticity();
            if(theInterface.getEnergyDissipator()!=null)a0 = theInterface.getEnergyDissipator().getAlpha0();
            if(Math.abs(k0)>1.e-10)r = theInterface.getTangentSpring().get_r();


            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double s1,s2;
            double N1,N2;
            double dN1,dN2;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                s1=0.; s2=0.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                if(theInterface.IsSlipIncluded()){
                    s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
                    s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}

                val*=0.5*kt*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                        2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                        (N1*s1+N2*s2)*(N1*s1+N2*s2)
                        )*(N1*z1+N2*z2)
                        +0.5*kh*(N1*s1+N2*s2 )*(N1*s1+N2*s2 )+k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r
                        +ge*Math.pow(Math.abs(dN1*s1+dN2*s2), 2.0)
                        +0.5*kt_0*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                        2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                        (N1*s1+N2*s2)*(N1*s1+N2*s2)
                        )
                        -a0*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }
    
    public double IntegrateTangentialStoredIEnergy(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kt = 0.0;
        double k0 = 0.0; double r=1.;
        double a0 = 0.0;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsTangentialIncluded() && theInterface.getTangentSpring()!=null){
            kt = theInterface.getTangentSpring().getStiffness();
            k0 = theInterface.getTangentSpring().get_GradientZ();
            if(theInterface.getEnergyDissipator()!=null)a0 = theInterface.getEnergyDissipator().getAlpha0();
            if(Math.abs(k0)>1.e-10)r = theInterface.getTangentSpring().get_r();


            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double N1,N2;
            double dN1,dN2;
            double va1,va2;
            double vb1 = 0,vb2 = 0;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }
                }
                pstep=wstep-1; if(wstep==0)pstep=0;
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getut_(vMaterial_main,tau)[pstep];
                va2 = elem.getINodeHier(2).getut_(vMaterial_main,tau)[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];


                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}
                double nu = vMaterial_main.getDispRelaxationTime_1();
                val*=0.5*kt*(tau*tau/((tau+nu)*(tau+nu)))*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)
                        )*(N1*z1+N2*z2)
                        +0.5*kt*(nu*nu/((tau+nu)*(tau+nu)))*(
                        (N1*va1+N2*va2)*(N1*va1+N2*va2)+
                        (N1*vb1+N2*vb2)*(N1*vb1+N2*vb2)-
                        2.*(N1*va1+N2*va2)*(N1*vb1+N2*vb2)
                        )*(N1*z1+N2*z2)
                        +kt*(tau*nu/((tau+nu)*(tau+nu)))*(
                        (N1*ua1+N2*ua2-N1*ub1-N2*ub2)*(N1*va1+N2*va2-N1*vb1-N2*vb2))*(N1*z1+N2*z2)
                        +k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r
                        -a0*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val*((tau+nu)/tau);
            }
        }
        return stored;
    }
    
    public double IntegrateTangentialStoredIEnergyQuad(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kt = 0.0;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsTangentialIncluded() && theInterface.getTangentSpring()!=null){
            kt = theInterface.getTangentSpring().getStiffness();
            
            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double N1,N2;
            double dN1,dN2;
            double va1,va2;
            double vb1 = 0,vb2 = 0;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }
                }
                
                pstep=wstep-1; if(wstep==0)pstep=0;
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getut_(vMaterial_main,tau)[pstep];
                va2 = elem.getINodeHier(2).getut_(vMaterial_main,tau)[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}
                double nu = vMaterial_main.getDispRelaxationTime_1();
                val*=0.5*kt*(tau*tau/((tau+nu)*(tau+nu)))*(
                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)
                        )*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val*((tau+nu)/tau);
            }
        }
        return stored;
    }
    
    public double IntegrateTangentialStoredIEnergyLinear(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kt = 0.0;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsTangentialIncluded() && theInterface.getTangentSpring()!=null){
            kt = theInterface.getTangentSpring().getStiffness();
            
            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double N1,N2;
            double dN1,dN2;
            double va1,va2;
            double vb1 = 0,vb2 = 0;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }
                }
                
                pstep=wstep-1; if(wstep==0)pstep=0;
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getut_(vMaterial_main,tau)[pstep];
                va2 = elem.getINodeHier(2).getut_(vMaterial_main,tau)[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}
                double nu = vMaterial_main.getDispRelaxationTime_1();
                val*=kt*(tau*nu/((tau+nu)*(tau+nu)))*(
                        (N1*ua1+N2*ua2-N1*ub1-N2*ub2)*(N1*va1+N2*va2-N1*vb1-N2*vb2))*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val*((tau+nu)/tau);
            }
        }
        return stored;
    }

    public double IntegrateDissipatedDamageIEnergy(IELine elem, Interface theInterface, int wstep){
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double aI=theInterface.getEnergyDissipator().getAlphaI();
        double e=theInterface.getEnergyDissipator().getEpsilon();

        double z1 = 1.,z2 = 1.;
        double zp1 = 1.,zp2 = 1.;
        double dz1,dz2;
        
        double uta1,uta2;
        double utb1,utb2;
        
        double una1,una2;
        double unb1,unb2;
        
        double kt=0.,kn;
        double kh;
        
        double lam = theInterface.getEnergyDissipator().getLambda();
        
        int model = theInterface.getEnergyDissipator().getModel();

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            utb1=0.;utb2=0.;
            uta1=0.;uta2=0.;
            unb1=0.;unb2=0.;
            una1=0.;una2=0.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];

            if(theInterface.IsDamageIncluded()){
                z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                if(wstep>0){
                    zp1 = elem.getINodeHier(1).getz(elem.getID())[wstep-1];
                    zp2 = elem.getINodeHier(2).getz(elem.getID())[wstep-1];
                }else{
                    zp1 = elem.getINodeHier(1).getz(elem.getID())[0];
                    zp2 = elem.getINodeHier(2).getz(elem.getID())[0];
                }
            }
            
            int pstep=wstep-1; if(wstep==0)pstep=0;
            if(theInterface.IsTangentialIncluded()){
                kt = theInterface.getTangentSpring().getStiffness();
                uta1 = elem.getINodeHier(1).getut()[pstep];
                uta2 = elem.getINodeHier(2).getut()[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        uta1 = elem.getINodeHier(1).getut()[pstep];
                        utb1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                        uta2 = elem.getINodeHier(2).getut()[pstep];
                        utb2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                    }else{
                        utb1 = elem.getINodeHier(1).getut()[pstep];
                        uta1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                        utb2 = elem.getINodeHier(2).getut()[pstep];
                        uta2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                    }
                }
                //||ut||=(N1*uta1+N2*uta2-N1*utb1-N2*utb2)
            }
            kn = theInterface.getNormalSpring().getStiffness();
            una1 = elem.getINodeHier(1).getun()[pstep];
            una2 = elem.getINodeHier(2).getun()[pstep];
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                    una1 = elem.getINodeHier(1).getun()[pstep];
                    unb1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                    una2 = elem.getINodeHier(2).getun()[pstep];
                    unb2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                }else{
                    unb1 = elem.getINodeHier(1).getun()[pstep];
                    una1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                    unb2 = elem.getINodeHier(2).getun()[pstep];
                    una2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                }
            }
            //||un||=(N1*una1+N2*una2-N1*unb1-N2*unb2)
            

            dz1=Math.abs(zp1-z1);
            dz2=Math.abs(zp2-z2);
            
            double N1=elem.getShapeFunction(1, xsi);
            double N2=elem.getShapeFunction(2, xsi);
            
            double psi;
            switch(theInterface.getEnergyDissipator().getFractureAngle()){
                case 1: 
                    psi = Math.atan2((N1*uta1+N2*uta2-N1*utb1-N2*utb2),(N1*una1+N2*una2-N1*unb1-N2*unb2));
                    break;
                case 2:
                    psi = (kt/kn)*Math.atan2((N1*uta1+N2*uta2-N1*utb1-N2*utb2),(N1*una1+N2*una2-N1*unb1-N2*unb2));
                    break;
                default:
                    psi = Math.atan2(Math.sqrt(kt*(N1*uta1+N2*uta2-N1*utb1-N2*utb2)*(N1*uta1+N2*uta2-N1*utb1-N2*utb2)),Math.sqrt(kn*(N1*una1+N2*una2-N1*unb1-N2*unb2)*(N1*una1+N2*una2-N1*unb1-N2*unb2)));
                    break;
            }

            switch(model){
                case 1: val*=aI*(1.+Math.tan((1.-lam)*psi)*Math.tan((1.-lam)*psi))*(N1*dz1+N2*dz2)
                    +e*(N1*dz1+N2*dz2)*(N1*dz1+N2*dz2);
                    break;
                case 2: 
                    throw new UnsupportedOperationException("Not supported yet.");
                case 3:
                    val*=theInterface.getEnergyDissipator().getAlpha(theInterface.getTangentSpring(), theInterface.getNormalSpring(), model, psi)*(N1*dz1+N2*dz2)
                    +e*(N1*dz1+N2*dz2)*(N1*dz1+N2*dz2);
                    break;
                default : val*=aI*(N1*dz1+N2*dz2)
                    +e*(N1*dz1+N2*dz2)*(N1*dz1+N2*dz2);
                    break;
            }
            
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }

    public double IntegrateDissipatedSlipIEnergy(IELine elem, Interface theInterface, int wstep){
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double a2=theInterface.getEnergyDissipator().getAlpha2();
        double a2_1=theInterface.getEnergyDissipator().getAlpha2_1();

        double s1 = 0.,s2 = 0.;
        double sp1 = 0.,sp2 = 0.;


        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }


        for(int i=0; i<numGP; i++){
            val=1.0;
            s1 = 0. ; s2 = 0.;
            sp1 = 0. ; sp2 = 0.;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];

            if(theInterface.IsSlipIncluded()){
                s1 = Math.abs(elem.getINodeHier(1).gets(elem.getID())[wstep]);
                s2 = Math.abs(elem.getINodeHier(2).gets(elem.getID())[wstep]);
                
                if(wstep>0){
                    sp1=Math.abs(elem.getINodeHier(1).gets(elem.getID())[wstep-1]);
                    sp2=Math.abs(elem.getINodeHier(2).gets(elem.getID())[wstep-1]);
                }
            }
            

            val*=a2*(elem.getShapeFunction(1, xsi)*Math.abs(s1-sp1)+elem.getShapeFunction(2, xsi)*Math.abs(s2-sp2));
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }
    
    public double IntegrateDissipatedSlipIEnergy_min(IELine elem, Interface theInterface, int wstep){
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double a2=theInterface.getEnergyDissipator().getAlpha2();
        double a2_1=theInterface.getEnergyDissipator().getAlpha2_1();
        
        double a1=theInterface.getEnergyDissipator().getAlphaI();


        double s_pos1 = 0.,s_pos2 = 0.;
        double s_neg1 = 0.,s_neg2 = 0.;
        double z1 = 1.,z2 = 1.;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }


        for(int i=0; i<numGP; i++){
            val=1.0;
            s_pos1 = 0. ; s_pos2 = 0.;
            s_neg1 = 0. ; s_neg2 = 0.;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            z1 = 1.; z2 = 1.;
            
            if(theInterface.IsDamageIncluded()){
                z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];
            }
            if(z1<=1.e-16){z1=0.0;}else{z1=1.0;}
            if(z2<=1.e-16){z2=0.0;}else{z2=1.0;}
//            z1 = 1.; z2 = 1.;

            if(theInterface.IsSlipIncluded()){
                s_pos1 = elem.getINodeHier(1).gets_pos(elem.getID())[wstep];
                s_pos2 = elem.getINodeHier(2).gets_pos(elem.getID())[wstep];
                
                s_neg1 = elem.getINodeHier(1).gets_neg(elem.getID())[wstep];
                s_neg2 = elem.getINodeHier(2).gets_neg(elem.getID())[wstep];
//                
//                if(wstep>0){
//                    s_pos1-=elem.getINodeHier(1).gets_pos(elem.getID())[wstep-1];
//                    s_pos2-=elem.getINodeHier(2).gets_pos(elem.getID())[wstep-1];
//                
//                    s_neg1-=-elem.getINodeHier(1).gets_neg(elem.getID())[wstep-1];
//                    s_neg2-=-elem.getINodeHier(2).gets_neg(elem.getID())[wstep-1];
//                    
////                    s_pos1=Math.abs(s_pos1);
////                    s_pos2=Math.abs(s_pos2);
////                    s_neg1=Math.abs(s_neg1);
////                    s_neg2=Math.abs(s_neg2);
//                }
            }

//            val*=a2*(elem.getShapeFunction(1, xsi)*s_pos1+elem.getShapeFunction(2, xsi)*s_pos2)*(elem.getShapeFunction(1, xsi)*z1+elem.getShapeFunction(2, xsi)*z2)+
//                    a2*(elem.getShapeFunction(1, xsi)*s_neg1+elem.getShapeFunction(2, xsi)*s_neg2)*(elem.getShapeFunction(1, xsi)*z1+elem.getShapeFunction(2, xsi)*z2);
            val*=a2*(elem.getShapeFunction(1, xsi)*s_pos1+elem.getShapeFunction(2, xsi)*s_pos2)+
                    a2*(elem.getShapeFunction(1, xsi)*s_neg1+elem.getShapeFunction(2, xsi)*s_neg2);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }

    public double IntegrateDissipatedPlastIEnergy(IELine elem, Interface theInterface, int wstep){
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double a2=theInterface.getEnergyDissipator().getAlpha3();

        double s1 = 0.,s2 = 0.;
        double sp1 = 0.,sp2 = 0.;
        double ds1,ds2;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            val=1.0;
            s1 = 0. ; s2 = 0.;
            sp1= 0. ; sp2=0.;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];

            if(theInterface.IsPlastIncluded()){
                s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                if(wstep>0){
                    sp1 = elem.getINodeHier(1).getp(elem.getID())[wstep-1];
                    sp2 = elem.getINodeHier(2).getp(elem.getID())[wstep-1];
                }
            }

            ds1=Math.abs(s1-sp1); //Math.abs
            ds2=Math.abs(s2-sp2); //Math.abs

            val*=a2*(elem.getShapeFunction(1, xsi)*ds1+elem.getShapeFunction(2, xsi)*ds2);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }
    
    public double IntegrateDissipatedPlastIEnergy_min(IELine elem, Interface theInterface, int wstep){
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double a2=theInterface.getEnergyDissipator().getAlpha3();

        double p_pos1 = 0.,p_pos2 = 0.;
        double p_neg1 = 0.,p_neg2 = 0.;
        double z1 = 1.,z2 = 1.;


        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }


        for(int i=0; i<numGP; i++){
            val=1.0;
            p_pos1 = 0. ; p_pos2 = 0.;
            p_neg1 = 0. ; p_neg2 = 0.;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            z1 = 1.; z2 = 1.;
            
            if(theInterface.IsDamageIncluded()){
                z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];
            }
            if(z1<=1.e-16){z1=0.0;}else{z1=1.0;}
            if(z2<=1.e-16){z2=0.0;}else{z2=1.0;}

            if(theInterface.IsSlipIncluded()){
                p_pos1 = elem.getINodeHier(1).getp_pos(elem.getID())[wstep];
                p_pos2 = elem.getINodeHier(2).getp_pos(elem.getID())[wstep];
                
                p_neg1 = elem.getINodeHier(1).getp_neg(elem.getID())[wstep];
                p_neg2 = elem.getINodeHier(2).getp_neg(elem.getID())[wstep];
                

            }

            val*=a2*(elem.getShapeFunction(1, xsi)*p_pos1+elem.getShapeFunction(2, xsi)*p_pos2)*(elem.getShapeFunction(1, xsi)*z1+elem.getShapeFunction(2, xsi)*z2)+
                    a2*(elem.getShapeFunction(1, xsi)*p_neg1+elem.getShapeFunction(2, xsi)*p_neg2)*(elem.getShapeFunction(1, xsi)*z1+elem.getShapeFunction(2, xsi)*z2);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }


    public double IntegrateIntGeneralisedNormalForce(IELine elem, Interface theInterface, int wstep, int nodeID){
        // the partial derivative of spring energy with respect to normal displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        // nodeID refers to the ID of the InterfaceNode
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kn = 0.;
        if(theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
        }

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1=0., s2=0;
        double du1,du2;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            s1=0.; s2=0;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            if(theInterface.IsPlastIncluded())s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
            if(theInterface.IsPlastIncluded())s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

            if(theInterface.getNumOfConnectedDomains()==2){                
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    s1=-s1; s2=-s2;
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getun()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ub2 = elem.getINodeHier(2).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                }
            }

            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);

            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=kn*N1*(N1*(ua1-ub1-s1)+N2*(ua2-ub2-s2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=kn*N2*(N1*(ua1-ub1-s1)+N2*(ua2-ub2-s2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=kn*N1*(N1*(ua1-s1)+N2*(ua2-s2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=kn*N2*(N1*(ua1-s1)+N2*(ua2-s2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val;
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedNormalForce(IELine elem, Interface theInterface, int wstep, int nodeID, double tau){
        // the partial derivative of spring energy with respect to normal displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        // nodeID refers to the ID of the InterfaceNode
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kn = 0.;
        if(theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
        }

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            if(theInterface.getNumOfConnectedDomains()==2){                
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getun()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ub2 = elem.getINodeHier(2).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                }
            }
            
            pstep=wstep-1;if(wstep==0)pstep=0;
            vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
            if(theInterface.getNumOfConnectedDomains()==2){                
                throw new UnsupportedOperationException("Not supported yet.");
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    va1 = elem.getINodeHier(1).getun_(vMaterial_main, tau)[pstep];
                    va2 = elem.getINodeHier(2).getun_(vMaterial_main, tau)[pstep];
                }
            }
            
            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            double nu = vMaterial_main.getDispRelaxationTime_1();
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N1*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2)
                            +(tau*nu/((tau+nu)*(tau+nu)))*kn*N1*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N2*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2)
                            +(tau*nu/((tau+nu)*(tau+nu)))*kn*N2*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N1*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2)
                                                            +(tau*nu/((tau+nu)*(tau+nu)))*kn*N1*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N2*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2)
                                                            +(tau*nu/((tau+nu)*(tau+nu)))*kn*N2*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val*((tau+nu)/tau);
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedNormalForceQuad(IELine elem, Interface theInterface, int wstep, int nodeID, double tau){
        // the partial derivative of spring energy with respect to normal displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        // nodeID refers to the ID of the InterfaceNode
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kn = 0.;
        if(theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
        }

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            if(theInterface.getNumOfConnectedDomains()==2){                
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getun()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ub2 = elem.getINodeHier(2).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                }
            }
            
            pstep=wstep-1;if(wstep==0)pstep=0;
            vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
            if(theInterface.getNumOfConnectedDomains()==2){                
                throw new UnsupportedOperationException("Not supported yet.");
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    va1 = elem.getINodeHier(1).getun_(vMaterial_main, tau)[pstep];
                    va2 = elem.getINodeHier(2).getun_(vMaterial_main, tau)[pstep];
                }
            }
            
            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            double nu = vMaterial_main.getDispRelaxationTime_1();
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N1*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N2*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N1*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kn*N2*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val*((tau+nu)/tau);
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedNormalForceLinear(IELine elem, Interface theInterface, int wstep, int nodeID, double tau){
        // the partial derivative of spring energy with respect to normal displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        // nodeID refers to the ID of the InterfaceNode
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kn = 0.;
        if(theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
        }

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            if(theInterface.getNumOfConnectedDomains()==2){                
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getun()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                    ub2 = elem.getINodeHier(2).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                }
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getun()[wstep];
                    ua2 = elem.getINodeHier(2).getun()[wstep];
                }
            }
            
            pstep=wstep-1;if(wstep==0)pstep=0;
            vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
            if(theInterface.getNumOfConnectedDomains()==2){                
                throw new UnsupportedOperationException("Not supported yet.");
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    va1 = elem.getINodeHier(1).getun_(vMaterial_main, tau)[pstep];
                    va2 = elem.getINodeHier(2).getun_(vMaterial_main, tau)[pstep];
                }
            }
            
            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            double nu = vMaterial_main.getDispRelaxationTime_1();
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=(tau*nu/((tau+nu)*(tau+nu)))*kn*N1*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=(tau*nu/((tau+nu)*(tau+nu)))*kn*N2*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=(tau*nu/((tau+nu)*(tau+nu)))*kn*N1*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=(tau*nu/((tau+nu)*(tau+nu)))*kn*N2*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val*((tau+nu)/tau);
        }
        return stored;
    }

    public double IntegrateIntGeneralisedTangentForce(IELine elem, Interface theInterface, int wstep, int nodeID){
        // the partial derivative of spring energy with respect to tangential displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kt = 0.;
        if(theInterface.getTangentSpring()!=null)kt = theInterface.getTangentSpring().getStiffness();
        double kt_0 = 0.;
        if(theInterface.getTangentSpring()!=null)kt_0 = theInterface.getTangentSpring().getStiffness_0();

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 = 0.,s2 = 0.;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            s1=0.; s2=0.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];


            if(theInterface.IsSlipIncluded())s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
            if(theInterface.IsSlipIncluded())s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];

            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    s1=-s1; s2=-s2;
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getut()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ub2 = elem.getINodeHier(2).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }
            }

            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                }
            }

            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);

            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=kt*N1*(N1*(ua1-ub1-s1)+N2*(ua2-ub2-s2))*(N1*z1+N2*z2)+
                            kt_0*N1*(N1*(ua1-ub1-s1)+N2*(ua2-ub2-s2));
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=kt*N2*(N1*(ua1-ub1-s1)+N2*(ua2-ub2-s2))*(N1*z1+N2*z2)
                            +kt_0*N2*(N1*(ua1-ub1-s1)+N2*(ua2-ub2-s2));
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=kt*N1*(N1*(ua1-s1)+N2*(ua2-s2))*(N1*z1+N2*z2)+kt_0*N1*(N1*(ua1-s1)+N2*(ua2-s2));
               if(elem.getINodeHier(2).getID() == nodeID) val*=kt*N2*(N1*(ua1-s1)+N2*(ua2-s2))*(N1*z1+N2*z2)+kt_0*N2*(N1*(ua1-s1)+N2*(ua2-s2));

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val;
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedTangentForce(IELine elem, Interface theInterface, int wstep, int nodeID, double tau){
        // the partial derivative of spring energy with respect to tangential displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kt = 0.;
        if(theInterface.getTangentSpring()!=null)kt = theInterface.getTangentSpring().getStiffness();

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];

            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getut()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ub2 = elem.getINodeHier(2).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }
            }

            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                }
            }
            
            pstep=wstep-1; if(wstep==0)pstep=0;
            vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
            if(theInterface.getNumOfConnectedDomains()==2){                
                throw new UnsupportedOperationException("Not supported yet.");
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    va1 = elem.getINodeHier(1).getut_(vMaterial_main, tau)[pstep];
                    va2 = elem.getINodeHier(2).getut_(vMaterial_main, tau)[pstep];
                }
            }

            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            double nu = vMaterial_main.getDispRelaxationTime_1();
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N1*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2)
                            +(tau*nu/((tau+nu)*(tau+nu)))*kt*N1*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N2*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2)
                            +(tau*nu/((tau+nu)*(tau+nu)))*kt*N2*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N1*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2)
                                                            +(tau*nu/((tau+nu)*(tau+nu)))*kt*N1*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N2*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2)
                                                            +(tau*nu/((tau+nu)*(tau+nu)))*kt*N2*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val*((tau+nu)/tau);
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedTangentForceQuad(IELine elem, Interface theInterface, int wstep, int nodeID, double tau){
        // the partial derivative of spring energy with respect to tangential displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kt = 0.;
        if(theInterface.getTangentSpring()!=null)kt = theInterface.getTangentSpring().getStiffness();

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];

            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getut()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ub2 = elem.getINodeHier(2).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }
            }

            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                }
            }
            
            pstep=wstep-1; if(wstep==0)pstep=0;
            vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
            if(theInterface.getNumOfConnectedDomains()==2){                
                throw new UnsupportedOperationException("Not supported yet.");
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    va1 = elem.getINodeHier(1).getut_(vMaterial_main, tau)[pstep];
                    va2 = elem.getINodeHier(2).getut_(vMaterial_main, tau)[pstep];
                }
            }

            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            double nu = vMaterial_main.getDispRelaxationTime_1();
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N1*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N2*(N1*(ua1-ub1)+N2*(ua2-ub2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N1*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=(tau*tau/((tau+nu)*(tau+nu)))*kt*N2*(N1*(ua1)+N2*(ua2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val*((tau+nu)/tau);
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedTangentForceLinear(IELine elem, Interface theInterface, int wstep, int nodeID, double tau){
        // the partial derivative of spring energy with respect to tangential displacement
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double kt = 0.;
        if(theInterface.getTangentSpring()!=null)kt = theInterface.getTangentSpring().getStiffness();

        double z1 = 1.,z2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double va1=0.,va2=0.;
        double vb1=0.,vb2=0.;
        int pstep;
        ViscousMaterial vMaterial_main;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(int i=0; i<numGP; i++){
            ua1=0.;ua2=0.;
            ub1=0.;ub2=0.;
            z1=1.; z2=1.;
            val=1.0;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];

            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                    ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }else if(elem.getINodeHier(1).getTwinID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    ub1 = elem.getINodeHier(1).getut()[wstep];
                    ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                    ub2 = elem.getINodeHier(2).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                }
            }

            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    ua1 = elem.getINodeHier(1).getut()[wstep];
                    ua2 = elem.getINodeHier(2).getut()[wstep];
                }
            }
            
            pstep=wstep-1; if(wstep==0)pstep=0;
            vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
            if(theInterface.getNumOfConnectedDomains()==2){                
                throw new UnsupportedOperationException("Not supported yet.");
            }
            
            if(theInterface.getNumOfConnectedDomains()==1){
                if(elem.getINodeHier(1).getID()  == nodeID || elem.getINodeHier(2).getID() == nodeID){
                    va1 = elem.getINodeHier(1).getut_(vMaterial_main, tau)[pstep];
                    va2 = elem.getINodeHier(2).getut_(vMaterial_main, tau)[pstep];
                }
            }

            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            double nu = vMaterial_main.getDispRelaxationTime_1();
            if(theInterface.getNumOfConnectedDomains()==2){
                if(elem.getINodeHier(1).getID() == nodeID || elem.getINodeHier(1).getTwinID() == nodeID){
                    val*=(tau*nu/((tau+nu)*(tau+nu)))*kt*N1*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID || elem.getINodeHier(2).getTwinID() == nodeID){
                    val*=(tau*nu/((tau+nu)*(tau+nu)))*kt*N2*(N1*(va1-vb1)+N2*(va2-vb2))*(N1*z1+N2*z2);
                }
            }
            if(theInterface.getNumOfConnectedDomains()==1){
               if(elem.getINodeHier(1).getID() == nodeID) val*=(tau*nu/((tau+nu)*(tau+nu)))*kt*N1*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);
               if(elem.getINodeHier(2).getID() == nodeID) val*=(tau*nu/((tau+nu)*(tau+nu)))*kt*N2*(N1*(va1)+N2*(va2))*(N1*z1+N2*z2);

            }
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            stored+=val*((tau+nu)/tau);
        }
        return stored;
    }

    public double IntegrateDissipatedDamageIForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of dissipator energy with respect to damage variable
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double aI=theInterface.getEnergyDissipator().getAlphaI();
        double e=theInterface.getEnergyDissipator().getEpsilon();

        double z1 = 1.,z2 = 1.;
        double zp1 = 1.,zp2 = 1.;
        double N1,N2;
        //double sign1,sign2;
        
        double uta1,uta2;
        double utb1,utb2;
        
        double una1,una2;
        double unb1,unb2;
        
        double kt=0.,kn;
        
        double lam = theInterface.getEnergyDissipator().getLambda();
        
        int model = theInterface.getEnergyDissipator().getModel();
        
        double psi = 0;
        if((elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getzEFTable(elem.getID()) == dof)||
                (elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getzEFTable(elem.getID()) == dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                utb1=0.;utb2=0.;
                uta1=0.;uta2=0.;
                unb1=0.;unb2=0.;
                una1=0.;una2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                int pstep=wstep-1; if(wstep==0)pstep=0;
                if(theInterface.IsTangentialIncluded()){
                    kt = theInterface.getTangentSpring().getStiffness();
                    uta1 = elem.getINodeHier(1).getut()[pstep];
                    uta2 = elem.getINodeHier(2).getut()[pstep];
                    if(theInterface.getNumOfConnectedDomains()==2){
                        if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                            uta1 = elem.getINodeHier(1).getut()[pstep];
                            utb1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                            uta2 = elem.getINodeHier(2).getut()[pstep];
                            utb2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                        }else{
                            utb1 = elem.getINodeHier(1).getut()[pstep];
                            uta1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                            utb2 = elem.getINodeHier(2).getut()[pstep];
                            uta2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                        }
                    }
                    //||ut||=(N1*uta1+N2*uta2-N1*utb1-N2*utb2)
                }
                kn = theInterface.getNormalSpring().getStiffness();
                una1 = elem.getINodeHier(1).getun()[pstep];
                una2 = elem.getINodeHier(2).getun()[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        una1 = elem.getINodeHier(1).getun()[pstep];
                        unb1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                        una2 = elem.getINodeHier(2).getun()[pstep];
                        unb2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                    }else{
                        unb1 = elem.getINodeHier(1).getun()[pstep];
                        una1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                        unb2 = elem.getINodeHier(2).getun()[pstep];
                        una2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                    }
                }
                //||un||=(N1*una1+N2*una2-N1*unb1-N2*unb2)

                zp1 = 1.; zp2 = 1.;
                N1=elem.getShapeFunction(1, xsi); N2=elem.getShapeFunction(2, xsi);
                if(theInterface.IsDamageIncluded()){

                    z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                    z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                    if(wstep>0){
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[wstep-1];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[wstep-1];
                    }else{
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[0];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[0];
                    }

                    //if(z1-zp1>=0){sign1=1.;}else{sign1=-1.;}
                    //if(z2-zp2>=0){sign2=1.;}else{sign2=-1.;}


                    switch(theInterface.getEnergyDissipator().getFractureAngle()){
                        case 1: 
                            psi = Math.atan2((N1*uta1+N2*uta2-N1*utb1-N2*utb2),(N1*una1+N2*una2-N1*unb1-N2*unb2));
                            break;
                        case 2:
                            psi = (kt/kn)*Math.atan2((N1*uta1+N2*uta2-N1*utb1-N2*utb2),(N1*una1+N2*una2-N1*unb1-N2*unb2));
                            break;
                        default:
                            psi = Math.atan2(Math.sqrt(kt*(N1*uta1+N2*uta2-N1*utb1-N2*utb2)*(N1*uta1+N2*uta2-N1*utb1-N2*utb2)),Math.sqrt(kn*(N1*una1+N2*una2-N1*unb1-N2*unb2)*(N1*una1+N2*una2-N1*unb1-N2*unb2)));
                            break;
                    }

                    if(model==2||model==3){
                        double aaI=theInterface.getEnergyDissipator().getAlpha(theInterface.getTangentSpring(), theInterface.getNormalSpring(), model, psi);
                        aI=aaI;
                    }
                    if(model==1){
                        double aaI=theInterface.getEnergyDissipator().getAlphaI()*(1.+Math.tan((1.-lam)*psi)*Math.tan((1.-lam)*psi));
                        aI=aaI;
                    }
//                    if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getzEFTable(elem.getID()) == dof)val*=-aI*N1*sign1+2*e*N1*(N1*(zp1-z1)+N2*(zp2-z2));
//                    if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getzEFTable(elem.getID()) == dof)val*=-aI*N2*sign2+2*e*N2*(N1*(zp1-z1)+N2*(zp2-z2));
                    if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getzEFTable(elem.getID()) == dof)val*=-aI*N1+2*e*N1*(N1*(zp1-z1)+N2*(zp2-z2));
                    if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getzEFTable(elem.getID()) == dof)val*=aI*N2+2*e*N2*(N1*(zp1-z1)+N2*(zp2-z2));

                }
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                dis+=val;
            }
        }

        return dis;
    }
    
    public double IntegrateDissipatedDamageIConstant(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of dissipator energy with respect to damage variable
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double aI=theInterface.getEnergyDissipator().getAlphaI();
        double e=theInterface.getEnergyDissipator().getEpsilon();

        double z1 = 1.,z2 = 1.;
        double zp1 = 1.,zp2 = 1.;
        double N1,N2;
        double sign1,sign2;
        
        double uta1,uta2;
        double utb1,utb2;
        
        double una1,una2;
        double unb1,unb2;
        
        double kt=0.,kn;
        
        double lam = theInterface.getEnergyDissipator().getLambda();
        
        int model = theInterface.getEnergyDissipator().getModel();
        
        double psi;
        
        if((elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getzEFTable(elem.getID()) == dof)||
                (elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getzEFTable(elem.getID()) == dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                utb1=0.;utb2=0.;
                uta1=0.;uta2=0.;
                unb1=0.;unb2=0.;
                una1=0.;una2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                int pstep=wstep-1; if(wstep==0)pstep=0;
                if(theInterface.IsTangentialIncluded()){
                    kt = theInterface.getTangentSpring().getStiffness();
                    uta1 = elem.getINodeHier(1).getut()[pstep];
                    uta2 = elem.getINodeHier(2).getut()[pstep];
                    if(theInterface.getNumOfConnectedDomains()==2){
                        if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                            uta1 = elem.getINodeHier(1).getut()[pstep];
                            utb1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                            uta2 = elem.getINodeHier(2).getut()[pstep];
                            utb2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                        }else{
                            utb1 = elem.getINodeHier(1).getut()[pstep];
                            uta1 = elem.getINodeHier(1).getTwin().getut()[pstep];
                            utb2 = elem.getINodeHier(2).getut()[pstep];
                            uta2 = elem.getINodeHier(2).getTwin().getut()[pstep];
                        }
                    }
                    //||ut||=(N1*uta1+N2*uta2-N1*utb1-N2*utb2)
                }
                kn = theInterface.getNormalSpring().getStiffness();
                una1 = elem.getINodeHier(1).getun()[pstep];
                una2 = elem.getINodeHier(2).getun()[pstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        una1 = elem.getINodeHier(1).getun()[pstep];
                        unb1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                        una2 = elem.getINodeHier(2).getun()[pstep];
                        unb2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                    }else{
                        unb1 = elem.getINodeHier(1).getun()[pstep];
                        una1 = elem.getINodeHier(1).getTwin().getun()[pstep];
                        unb2 = elem.getINodeHier(2).getun()[pstep];
                        una2 = elem.getINodeHier(2).getTwin().getun()[pstep];
                    }
                }
                //||un||=(N1*una1+N2*una2-N1*unb1-N2*unb2)

                zp1 = 1.; zp2 = 1.;
                N1=elem.getShapeFunction(1, xsi); N2=elem.getShapeFunction(2, xsi);
                if(theInterface.IsDamageIncluded()){

                    z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                    z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                    if(wstep>0){
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[wstep-1];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[wstep-1];
                    }else{
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[0];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[0];
                    }

                    if(z1-zp1>=0){sign1=1.;}else{sign1=-1.;}
                    if(z2-zp2>=0){sign2=1.;}else{sign2=-1.;}

                    switch(theInterface.getEnergyDissipator().getFractureAngle()){
                        case 1: 
                            psi = Math.atan2((N1*uta1+N2*uta2-N1*utb1-N2*utb2),(N1*una1+N2*una2-N1*unb1-N2*unb2));
                            break;
                        case 2:
                            psi = (kt/kn)*Math.atan2((N1*uta1+N2*uta2-N1*utb1-N2*utb2),(N1*una1+N2*una2-N1*unb1-N2*unb2));
                            break;
                        default:
                            psi = Math.atan2(Math.sqrt(kt*(N1*uta1+N2*uta2-N1*utb1-N2*utb2)*(N1*uta1+N2*uta2-N1*utb1-N2*utb2)),Math.sqrt(kn*(N1*una1+N2*una2-N1*unb1-N2*unb2)*(N1*una1+N2*una2-N1*unb1-N2*unb2)));
                            break;
                    }

                    if(model==2||model==3){
                        double aaI=theInterface.getEnergyDissipator().getAlpha(theInterface.getTangentSpring(), theInterface.getNormalSpring(), model, psi);
                        aI=aaI;
                    }
                    if(model==1){
                        double aaI=theInterface.getEnergyDissipator().getAlphaI()*(1.+Math.tan((1.-lam)*psi)*Math.tan((1.-lam)*psi));
                        aI=aaI;
                    }
                    if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getzEFTable(elem.getID()) == dof)val*=aI*N1;
                    if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getzEFTable(elem.getID()) == dof)val*=aI*N2;

                }

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                dis+=val;
            }
        }

        return dis;
    }

    public double IntegrateIntGeneralisedNormalDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of spring energy with respect to damage variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;
        double sign=1.;

        double kn = theInterface.getNormalSpring().getStiffness();
        double k0 = theInterface.getNormalSpring().get_GradientZ();
        double r = 2.;
        if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();

        double z1 ,z2 ;
        double zp1 = 1.,zp2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        if((elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getzEFTable(elem.getID()) == dof)||
                (elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getzEFTable(elem.getID()) == dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zp1 = 1.; zp2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsPlastIncluded())s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                if(theInterface.IsPlastIncluded())s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];

                if(theInterface.getNumOfConnectedDomains()==2){

                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                        //s1=-s1;s2=-s2;
                    }

                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                if( (dN1*z1+dN2*z2)>0 ){sign=1.;}else if( (dN1*z1+dN2*z2)<0 ){sign=-1.;}else{sign=0.;}

                if(elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof){
                    val*=0.5*kn*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N1+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N1*sign;
                }
                if(elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof){
                    val*=0.5*kn*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N2+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N2*sign;
                }
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }
    
//    public double IntegrateIntGeneralisedNormalDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof, double nu, double tau){
//        // the partial derivative of spring energy with respect to damage variable
//        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
//        double stored=0.;
//        double[] CoordsOnElem=new double[1];
//        double[] weight = new double[1];
//        double xsi;
//        double val;
//        double N1,N2;
//        double dN1,dN2;
//        double sign=1.;
//
//        double kn = theInterface.getNormalSpring().getStiffness();
//        double k0 = theInterface.getNormalSpring().get_GradientZ();
//        double r = 2.;
//        if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();
//
//        double z1 ,z2 ;
//        double zp1 = 1.,zp2 = 1.;
//        double ua1=0.,ua2=0.;
//        double ub1=0.,ub2=0.;
//        double va1=0.,va2=0.;
//        double vb1=0.,vb2=0.;
//
//        double[] xsis = new double[numGP];
//        double[] weights = new double[numGP];
//        switch(numGP){
//            case 3:
//                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
//                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
//                break;
//            case 4:
//                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
//                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
//                break;
//            case 5:
//                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
//                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
//                break;
//            default: System.exit(345); break;
//        }
//
//        for(int i=0; i<numGP; i++){
//            z1=0.; z2=0.;
//            zp1 = 1.; zp2 = 1.;
//            ua1=0.;ua2=0.;
//            ub1=0.;ub2=0.;
//            val=1.0;
//            CoordsOnElem[0]=xsis[i];
//            weight[0]=weights[i];
//            xsi=CoordsOnElem[0];
//
//            if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
//            if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];
//
//            ua1 = elem.getINodeHier(1).getun()[wstep];
//            ua2 = elem.getINodeHier(2).getun()[wstep];
//
//            if(theInterface.getNumOfConnectedDomains()==2){
//
//                if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
//                    ua1 = elem.getINodeHier(1).getun()[wstep];
//                    ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
//                    ua2 = elem.getINodeHier(2).getun()[wstep];
//                    ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
//                }else{
//                    ub1 = elem.getINodeHier(1).getun()[wstep];
//                    ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
//                    ub2 = elem.getINodeHier(2).getun()[wstep];
//                    ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
//                    //s1=-s1;s2=-s2;
//                }
//
//            }
//            va1 = elem.getINodeHier(1).getun_(nu,tau)[wstep-1];
//            va2 = elem.getINodeHier(2).getun_(nu,tau)[wstep-1];
//            if(theInterface.getNumOfConnectedDomains()==2){
//                if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
//                    va1 = elem.getINodeHier(1).getun_(nu,tau)[wstep-1];
//                    vb1 = elem.getINodeHier(1).getTwin().getun_(nu,tau)[wstep-1];
//                    va2 = elem.getINodeHier(2).getun_(nu,tau)[wstep-1];
//                    vb2 = elem.getINodeHier(2).getTwin().getun_(nu,tau)[wstep-1];
//                }else{
//                    vb1 = elem.getINodeHier(1).getun_(nu,tau)[wstep-1];
//                    va1 = elem.getINodeHier(1).getTwin().getun_(nu,tau)[wstep-1];
//                    vb2 = elem.getINodeHier(2).getun_(nu,tau)[wstep-1];
//                    va2 = elem.getINodeHier(2).getTwin().getun_(nu,tau)[wstep-1];
//                }
//            }
//
//            N1=elem.getShapeFunction(1, xsi);
//            N2=elem.getShapeFunction(2, xsi);
//
//            dN1=elem.getShapeFunction_xsi(1, xsi);
//            dN2=elem.getShapeFunction_xsi(2, xsi);
//
//            if( (dN1*z1+dN2*z2)>0 ){sign=1.;}else if( (dN1*z1+dN2*z2)<0 ){sign=-1.;}else{sign=0.;}
//
//            if(elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof){
//                val*=0.5*kn*(tau*tau/((tau+nu)*(tau+nu)))*(
//                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
//                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
//                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)
//                        )*N1
//                        +0.5*kn*(nu*nu/((tau+nu)*(tau+nu)))*(
//                        (N1*va1+N2*va2)*(N1*va1+N2*va2)+
//                        (N1*vb1+N2*vb2)*(N1*vb1+N2*vb2)-
//                        2.*(N1*va1+N2*va2)*(N1*vb1+N2*vb2)
//                        )*N1
//                        +kn*(tau*nu/((tau+nu)*(tau+nu)))*(
//                        (N1*ua1+N2*ua2-N1*ub1-N2*ub2)*(N1*va1+N2*va2-N1*vb1-N2*vb2))*N1
//                        +k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N1*sign;
//            }
//            if(elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof){
//                val*=0.5*kn*(tau*tau/((tau+nu)*(tau+nu)))*(
//                        (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
//                        (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
//                        2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)
//                        )*N2
//                        +0.5*kn*(nu*nu/((tau+nu)*(tau+nu)))*(
//                        (N1*va1+N2*va2)*(N1*va1+N2*va2)+
//                        (N1*vb1+N2*vb2)*(N1*vb1+N2*vb2)-
//                        2.*(N1*va1+N2*va2)*(N1*vb1+N2*vb2)
//                        )*N2
//                        +kn*(tau*nu/((tau+nu)*(tau+nu)))*(
//                        (N1*ua1+N2*ua2-N1*ub1-N2*ub2)*(N1*va1+N2*va2-N1*vb1-N2*vb2))*N2
//                        +k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N2*sign;
//            }
//            val*=elem.getJacobian(xsi);
//            val*=weight[0];
//            stored+=val*((tau+nu)/tau);
//        }
//        return stored;
//    }
    
    public double IntegrateIntGeneralisedNormalDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof, double tau){
        // the partial derivative of spring energy with respect to damage variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;
        double sign=1.;
        ViscousMaterial vMaterial_main;

        double kn = theInterface.getNormalSpring().getStiffness();
        double a0 = theInterface.getEnergyDissipator().getAlpha0();
        double k0 = theInterface.getTangentSpring().get_GradientZ();
        double r = 2.;
        if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();

        double z1 ,z2 ;
        double zp1 = 1.,zp2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        if((elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof)||
                (elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zp1 = 1.; zp2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsSlipIncluded())s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
                if(theInterface.IsSlipIncluded())s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];


                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                ua1 = elem.getINodeHier(1).getun_(vMaterial_main,tau)[wstep];
                ua2 = elem.getINodeHier(2).getun_(vMaterial_main,tau)[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                if( (dN1*z1+dN2*z2)>0 ){sign=1.;}else if( (dN1*z1+dN2*z2)<0 ){sign=-1.;}else{sign=0.;}

                if(elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof){
                    val*=0.5*(kn)*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N1-a0*N1+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N1*sign;
                }
                if(elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof){
                    val*=0.5*(kn)*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N2-a0*N2+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N2*sign;
                }
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }
    
    public double IntegrateIntGeneralisedNormalDrivingForce_split(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of spring energy with respect to damage variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;
        double sign=1.;

        double kn = theInterface.getNormalSpring().getStiffness();
        double k0 = theInterface.getNormalSpring().get_GradientZ();
        double r = 2.;
        if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();

        double z1 ,z2 ;
        double zp1 = 1.,zp2 = 1.;
        double ua1pos=0.,ua2pos=0.;
        double ub1pos=0.,ub2pos=0.;
        
        double ua1neg=0.,ua2neg=0.;
        double ub1neg=0.,ub2neg=0.;
        double s1 ,s2 ;
        if((elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof)||
                (elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zp1 = 1.; zp2 = 1.;
                ua1pos=0.;ua2pos=0.;
                ub1pos=0.;ub2pos=0.;

                ua1neg=0.;ua2neg=0.;
                ub1neg=0.;ub2neg=0.;
                s1=0.; s2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsPlastIncluded())s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                if(theInterface.IsPlastIncluded())s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1pos = elem.getINodeHier(1).getun_pos()[wstep];
                ua2pos = elem.getINodeHier(2).getun_pos()[wstep];

                ua1neg = elem.getINodeHier(1).getun_neg()[wstep];
                ua2neg = elem.getINodeHier(2).getun_neg()[wstep];

                if(theInterface.getNumOfConnectedDomains()==2){

                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1pos = elem.getINodeHier(1).getun_pos()[wstep];
                        ub1pos = elem.getINodeHier(1).getTwin().getun_pos()[wstep];
                        ua2pos = elem.getINodeHier(2).getun_pos()[wstep];
                        ub2pos = elem.getINodeHier(2).getTwin().getun_pos()[wstep];

                        ua1neg = elem.getINodeHier(1).getun_neg()[wstep];
                        ub1neg = elem.getINodeHier(1).getTwin().getun_neg()[wstep];
                        ua2neg = elem.getINodeHier(2).getun_neg()[wstep];
                        ub2neg = elem.getINodeHier(2).getTwin().getun_neg()[wstep];
                    }else{
                        ub1pos = elem.getINodeHier(1).getun_pos()[wstep];
                        ua1pos = elem.getINodeHier(1).getTwin().getun_pos()[wstep];
                        ub2pos = elem.getINodeHier(2).getun_pos()[wstep];
                        ua2pos = elem.getINodeHier(2).getTwin().getun_pos()[wstep];

                        ub1neg = elem.getINodeHier(1).getun_neg()[wstep];
                        ua1neg = elem.getINodeHier(1).getTwin().getun_neg()[wstep];
                        ub2neg = elem.getINodeHier(2).getun_neg()[wstep];
                        ua2neg = elem.getINodeHier(2).getTwin().getun_neg()[wstep];
                        //s1=-s1; s2=-s2;
                    }

                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                if( (dN1*z1+dN2*z2)>0 ){sign=1.;}else if( (dN1*z1+dN2*z2)<0 ){sign=-1.;}else{sign=0.;}

                if(elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof){
                    val*=(kn*(-2*Math.pow(N1,3)*s1*ua1neg - 2*Math.pow(N1,2)*N2*s2*ua1neg + Math.pow(N1,3)*Math.pow(ua1neg,2) + 2*Math.pow(N1,3)*ua1neg*ua1pos - 2*Math.pow(N1,2)*N2*s1*ua2neg - 2*N1*Math.pow(N2,2)*s2*ua2neg +
                        2*Math.pow(N1,2)*N2*ua1neg*ua2neg + 2*Math.pow(N1,2)*N2*ua1pos*ua2neg + N1*Math.pow(N2,2)*Math.pow(ua2neg,2) + 2*Math.pow(N1,2)*N2*ua1neg*ua2pos + 2*N1*Math.pow(N2,2)*ua2neg*ua2pos +
                        2*Math.pow(N1,3)*s1*ub1neg + 2*Math.pow(N1,2)*N2*s2*ub1neg - 2*Math.pow(N1,3)*ua1neg*ub1neg - 2*Math.pow(N1,3)*ua1pos*ub1neg - 2*Math.pow(N1,2)*N2*ua2neg*ub1neg - 2*Math.pow(N1,2)*N2*ua2pos*ub1neg +
                        Math.pow(N1,3)*Math.pow(ub1neg,2) - 2*Math.pow(N1,3)*ua1neg*ub1pos - 2*Math.pow(N1,2)*N2*ua2neg*ub1pos + 2*Math.pow(N1,3)*ub1neg*ub1pos + 2*Math.pow(N1,2)*N2*s1*ub2neg + 2*N1*Math.pow(N2,2)*s2*ub2neg -
                        2*Math.pow(N1,2)*N2*ua1neg*ub2neg - 2*Math.pow(N1,2)*N2*ua1pos*ub2neg - 2*N1*Math.pow(N2,2)*ua2neg*ub2neg - 2*N1*Math.pow(N2,2)*ua2pos*ub2neg + 2*Math.pow(N1,2)*N2*ub1neg*ub2neg +
                        2*Math.pow(N1,2)*N2*ub1pos*ub2neg + N1*Math.pow(N2,2)*Math.pow(ub2neg,2) - 2*Math.pow(N1,2)*N2*ua1neg*ub2pos - 2*N1*Math.pow(N2,2)*ua2neg*ub2pos + 2*Math.pow(N1,2)*N2*ub1neg*ub2pos +
                        2*N1*Math.pow(N2,2)*ub2neg*ub2pos))/2.
                            +k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N1*sign;
                }
                if(elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof){
                    val*=(kn*(-2*Math.pow(N1,2)*N2*s1*ua1neg - 2*N1*Math.pow(N2,2)*s2*ua1neg + Math.pow(N1,2)*N2*Math.pow(ua1neg,2) + 2*Math.pow(N1,2)*N2*ua1neg*ua1pos - 2*N1*Math.pow(N2,2)*s1*ua2neg - 2*Math.pow(N2,3)*s2*ua2neg +
                        2*N1*Math.pow(N2,2)*ua1neg*ua2neg + 2*N1*Math.pow(N2,2)*ua1pos*ua2neg + Math.pow(N2,3)*Math.pow(ua2neg,2) + 2*N1*Math.pow(N2,2)*ua1neg*ua2pos + 2*Math.pow(N2,3)*ua2neg*ua2pos +
                        2*Math.pow(N1,2)*N2*s1*ub1neg + 2*N1*Math.pow(N2,2)*s2*ub1neg - 2*Math.pow(N1,2)*N2*ua1neg*ub1neg - 2*Math.pow(N1,2)*N2*ua1pos*ub1neg - 2*N1*Math.pow(N2,2)*ua2neg*ub1neg -
                        2*N1*Math.pow(N2,2)*ua2pos*ub1neg + Math.pow(N1,2)*N2*Math.pow(ub1neg,2) - 2*Math.pow(N1,2)*N2*ua1neg*ub1pos - 2*N1*Math.pow(N2,2)*ua2neg*ub1pos + 2*Math.pow(N1,2)*N2*ub1neg*ub1pos +
                        2*N1*Math.pow(N2,2)*s1*ub2neg + 2*Math.pow(N2,3)*s2*ub2neg - 2*N1*Math.pow(N2,2)*ua1neg*ub2neg - 2*N1*Math.pow(N2,2)*ua1pos*ub2neg - 2*Math.pow(N2,3)*ua2neg*ub2neg - 2*Math.pow(N2,3)*ua2pos*ub2neg +
                        2*N1*Math.pow(N2,2)*ub1neg*ub2neg + 2*N1*Math.pow(N2,2)*ub1pos*ub2neg + Math.pow(N2,3)*Math.pow(ub2neg,2) - 2*N1*Math.pow(N2,2)*ua1neg*ub2pos - 2*Math.pow(N2,3)*ua2neg*ub2pos +
                        2*N1*Math.pow(N2,2)*ub1neg*ub2pos + 2*Math.pow(N2,3)*ub2neg*ub2pos))/2.
                            +k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N2*sign;
                }
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }

    public double IntegrateIntGeneralisedTangentDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of spring energy with respect to damage variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;
        double sign=1.;

        double kt = theInterface.getTangentSpring().getStiffness();
        double kt_0 = theInterface.getTangentSpring().getStiffness_0();
        double a0 = theInterface.getEnergyDissipator().getAlpha0();
        double k0 = theInterface.getTangentSpring().get_GradientZ();
        double r = 2.;
        if(Math.abs(k0)>1.e-10)r = theInterface.getTangentSpring().get_r();

        double z1 ,z2 ;
        double zp1 = 1.,zp2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        if((elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof)||
                (elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zp1 = 1.; zp2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];


                if(theInterface.IsSlipIncluded()){
                    s1 = elem.getINodeHier(1).gets(elem.getID())[wstep]; 
                    s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];
                }


                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                        //s1=-s1; s2=-s2;
                    }
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                if( (dN1*z1+dN2*z2)>0 ){sign=1.;}else if( (dN1*z1+dN2*z2)<0 ){sign=-1.;}else{sign=0.;}

                if(elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof){
                    val*=0.5*(kt)*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N1-a0*N1+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N1*sign;
                }
                if(elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof){
                    val*=0.5*(kt)*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N2-a0*N2+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N2*sign;
                }
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }
    
    
    public double IntegrateIntGeneralisedTangentDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof, double tau){
        // the partial derivative of spring energy with respect to damage variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;
        double sign=1.;
        ViscousMaterial vMaterial_main;

        double kt = theInterface.getTangentSpring().getStiffness();
        double a0 = theInterface.getEnergyDissipator().getAlpha0();
        double k0 = theInterface.getTangentSpring().get_GradientZ();
        double r = 2.;
        if(Math.abs(k0)>1.e-10)r = theInterface.getTangentSpring().get_r();

        double z1 ,z2 ;
        double zp1 = 1.,zp2 = 1.;
        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        if((elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof)||
                (elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                z1=0.; z2=0.;
                zp1 = 1.; zp2 = 1.;
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsSlipIncluded())s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
                if(theInterface.IsSlipIncluded())s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];


                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                ua1 = elem.getINodeHier(1).getut_(vMaterial_main,tau)[wstep];
                ua2 = elem.getINodeHier(2).getut_(vMaterial_main,tau)[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                if( (dN1*z1+dN2*z2)>0 ){sign=1.;}else if( (dN1*z1+dN2*z2)<0 ){sign=-1.;}else{sign=0.;}

                if(elem.getINodeHier(1).getID() == nodeID && elem.getINodeHier(1).getzEFTable(elem.getID())==dof){
                    val*=0.5*(kt)*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N1-a0*N1+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N1*sign;
                }
                if(elem.getINodeHier(2).getID() == nodeID && elem.getINodeHier(2).getzEFTable(elem.getID())==dof){
                    val*=0.5*(kt)*(
                            (N1*ua1+N2*ua2)*(N1*ua1+N2*ua2)+
                            (N1*ub1+N2*ub2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*ub1+N2*ub2)-
                            2.*(N1*ua1+N2*ua2)*(N1*s1+N2*s2)+
                            2.*(N1*ub1+N2*ub2)*(N1*s1+N2*s2)+
                            (N1*s1+N2*s2)*(N1*s1+N2*s2)
                            )*N2-a0*N2+k0*(Math.pow((dN1*z1+dN2*z2), r-1.))*N2*sign;
                }
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }

    public double IntegrateDissipatedSlipIForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of dissipator energy with respect to slip variable
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double a2=theInterface.getEnergyDissipator().getAlpha2();
        double a2_1=theInterface.getEnergyDissipator().getAlpha2_1();

        double s1 = 0.,s2 = 0.;
        double sp1 = 0.,sp2 = 0.;
        double zp1 = 1.,zp2 = 1.;
        double ds1,ds2;
        double sign1 = 0,sign2 = 0;
        double N1,N2;
        if( (elem.getINodeHier(1).getID() == nodeID  && ((elem.getINodeHier(1).gets_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(1).gets_negEFTable(elem.getID()) == dof) ) )||
                (elem.getINodeHier(2).getID() == nodeID  && ((elem.getINodeHier(2).gets_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(2).gets_negEFTable(elem.getID()) == dof) ))){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ds1=0.; ds2=0.;
                sp1 = 0.; sp2 = 0.;
                zp1 = 0.; zp1 = 0.;

                if(theInterface.IsDamageIncluded()){
                    if(wstep>0){
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[wstep-1];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[wstep-1];
                    }else{
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[0];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[0];
                    }
                }

                if(theInterface.IsSlipIncluded()){
                    s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
                    s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];

                    if(wstep>0){
                        sp1 = elem.getINodeHier(1).gets(elem.getID())[wstep-1];
                        sp2 = elem.getINodeHier(2).gets(elem.getID())[wstep-1];
                    }
                }

                //pbug
                sign1=Math.signum(s1-sp1);
                sign2=Math.signum(s2-sp2);

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                sign1=Math.signum(N1*(s1-sp1)+N2*(s2-sp2)); sign2=sign1;

                if(elem.getINodeHier(1).getID() == nodeID  && ((elem.getINodeHier(1).gets_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(1).gets_negEFTable(elem.getID()) == dof) ))val*=a2*N1;
                if(elem.getINodeHier(2).getID() == nodeID  && ((elem.getINodeHier(2).gets_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(2).gets_negEFTable(elem.getID()) == dof) ))val*=a2*N2;

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                dis+=val;
            }
        }
        return dis;
    }
    
    public double IntegrateDissipatedSlipIConstant(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of dissipator energy with respect to slip variable
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double a2=theInterface.getEnergyDissipator().getAlpha2();
        double a2_1=theInterface.getEnergyDissipator().getAlpha2_1();

        double s1 = 0.,s2 = 0.;
        double sp1 = 0.,sp2 = 0.;
        double zp1 = 1.,zp2 = 1.;
        double ds1,ds2;
        double sign1 = 0,sign2 = 0;
        double N1,N2;
        if( (elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getsEFTable(elem.getID()) == dof)||
                (elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getsEFTable(elem.getID()) == dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ds1=0.; ds2=0.;
                sp1 = 0.; sp2 = 0.;
                zp1 = 0.; zp1 = 0.;

                if(theInterface.IsDamageIncluded()){
                    if(wstep>0){
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[wstep-1];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[wstep-1];
                    }else{
                        zp1 = elem.getINodeHier(1).getz(elem.getID())[0];
                        zp2 = elem.getINodeHier(2).getz(elem.getID())[0];
                    }
                }

                if(theInterface.IsSlipIncluded()){
                    s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
                    s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];

                    if(wstep>0){
                        sp1 = elem.getINodeHier(1).gets(elem.getID())[wstep-1];
                        sp2 = elem.getINodeHier(2).gets(elem.getID())[wstep-1];
                    }

                    if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getsEFTable(elem.getID()) == dof)ds1=1.;
                    if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getsEFTable(elem.getID()) == dof)ds2=1.;
                }

                sign1=Math.signum(s1-sp1);
                sign2=Math.signum(s2-sp2);

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                sign1=1.; sign2=sign1;

                if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getsEFTable(elem.getID()) == dof)val*=(a2+a2_1*(N1*zp1+N2*zp2))*N1*sign1;
                if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getsEFTable(elem.getID()) == dof)val*=(a2+a2_1*(N1*zp1+N2*zp2))*N2*sign2;

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                dis+=val;
            }
        }
        return dis;
    }

    public double IntegrateDissipatedPlastIForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of dissipator energy with respect to plast variable
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double a2=theInterface.getEnergyDissipator().getAlpha3();

        double s1 = 1.,s2 = 1.;
        double sp1 = 1.,sp2 = 1.;
        double ds1,ds2;
        double sign1 = 0,sign2 = 0;
        if((elem.getINodeHier(1).getID() == nodeID  && ((elem.getINodeHier(1).getp_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(1).getp_negEFTable(elem.getID()) == dof) ))||
                (elem.getINodeHier(2).getID() == nodeID  && ((elem.getINodeHier(2).getp_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(2).getp_negEFTable(elem.getID()) == dof) )) ){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ds1=0.; ds2=0.;
                sp1 = 0.; sp2 = 0.;

                if(theInterface.IsPlastIncluded()){
                    s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                    s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                    if(wstep>0){
                        sp1 = elem.getINodeHier(1).getp(elem.getID())[wstep-1];
                        sp2 = elem.getINodeHier(2).getp(elem.getID())[wstep-1];
                    }

                }

                //pbug
                sign1=Math.signum(s1-sp1);
                sign2=Math.signum(s2-sp2);

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                sign1=Math.signum(N1*(s1-sp1)+N2*(s2-sp2)); sign2=sign1;

                if(elem.getINodeHier(1).getID() == nodeID  && ((elem.getINodeHier(1).getp_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(1).getp_negEFTable(elem.getID()) == dof) ))val*=a2*N1;
                if(elem.getINodeHier(2).getID() == nodeID  && ((elem.getINodeHier(2).getp_posEFTable(elem.getID()) == dof)||(elem.getINodeHier(2).getp_negEFTable(elem.getID()) == dof) ))val*=a2*N2;

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                dis+=val;
            }
        }
        return dis;
    }
    
    public double IntegrateDissipatedPlastIConstant(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of dissipator energy with respect to plast variable
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;

        double a2=theInterface.getEnergyDissipator().getAlpha3();

        double s1 = 1.,s2 = 1.;
        double sp1 = 1.,sp2 = 1.;
        double ds1,ds2;
        double sign1 = 0,sign2 = 0;
        
        if((elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getpEFTable(elem.getID()) == dof)||
                (elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getpEFTable(elem.getID()) == dof)){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ds1=0.; ds2=0.;
                sp1 = 0.; sp2 = 0.;

                if(theInterface.IsPlastIncluded()){
                    s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                    s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                    if(wstep>0){
                        sp1 = elem.getINodeHier(1).getp(elem.getID())[wstep-1];
                        sp2 = elem.getINodeHier(2).getp(elem.getID())[wstep-1];
                    }

                    if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getpEFTable(elem.getID()) == dof)ds1=1.;
                    if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getpEFTable(elem.getID()) == dof)ds2=1.;
                }

                sign1=Math.signum(s1-sp1);
                sign2=Math.signum(s2-sp2);

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                sign1=1.; sign2=sign1;

                if(elem.getINodeHier(1).getID() == nodeID  && elem.getINodeHier(1).getpEFTable(elem.getID()) == dof)val*=(a2)*N1*sign1;
                if(elem.getINodeHier(2).getID() == nodeID  && elem.getINodeHier(2).getpEFTable(elem.getID()) == dof)val*=(a2)*N2*sign2;

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                dis+=val;
            }
        }
        return dis;
    }

    public double IntegrateIntGeneralisedSlipDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of spring energy with respect to slip variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;

        double kt = theInterface.getTangentSpring().getStiffness();
        double kt_0 = theInterface.getTangentSpring().getStiffness_0();
        double kh = theInterface.getTangentSpring().getIsotropicHardening();
        double ge = theInterface.getTangentSpring().get_GradientPlasticity();

        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        double z1 ,z2 ;
        if((elem.getINodeHier(1).getID() == nodeID && ((elem.getINodeHier(1).gets_posEFTable(elem.getID())==dof) || (elem.getINodeHier(1).gets_negEFTable(elem.getID())==dof)) )||
               (elem.getINodeHier(2).getID() == nodeID && ((elem.getINodeHier(2).gets_posEFTable(elem.getID())==dof) || (elem.getINodeHier(2).gets_negEFTable(elem.getID())==dof))) ){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }

            for(int i=0; i<numGP; i++){
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                z1=1.; z2=1.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsSlipIncluded())s1 = elem.getINodeHier(1).gets(elem.getID())[wstep];
                if(theInterface.IsSlipIncluded())s2 = elem.getINodeHier(2).gets(elem.getID())[wstep];

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];

                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                        //s1=-s1; s2=-s2;
                    }
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);

                if(elem.getINodeHier(1).getID() == nodeID && ((elem.getINodeHier(1).gets_posEFTable(elem.getID())==dof) || (elem.getINodeHier(1).gets_negEFTable(elem.getID())==dof)) ){
                    val*=kt*N1*(-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2)+(N1*s1+N2*s2))*(N1*z1+N2*z2)+
                            kt_0*N1*(-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2)+(N1*s1+N2*s2))+
                            kh*(N1*s1+N2*s2)*N1+ge*(dN1*s1+dN2*s2)*N1;
                }
                if(elem.getINodeHier(2).getID() == nodeID && ((elem.getINodeHier(2).gets_posEFTable(elem.getID())==dof) || (elem.getINodeHier(2).gets_negEFTable(elem.getID())==dof))){
                    val*=kt*N2*(-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2)+(N1*s1+N2*s2))*(N1*z1+N2*z2)+
                            kt_0*N2*(-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2)+(N1*s1+N2*s2))+
                            kh*(N1*s1+N2*s2)*N2+ge*(dN1*s1+dN2*s2)*N2;
                }


                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }

    public double IntegrateIntGeneralisedPlastDrivingForce(IELine elem, Interface theInterface, int wstep, int nodeID, int dof){
        // the partial derivative of spring energy with respect to plast variable
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double N1,N2;
        double dN1,dN2;

        double kn = theInterface.getNormalSpring().getStiffness();
        double kh = theInterface.getNormalSpring().getIsotropicHardening();
        double ge = theInterface.getNormalSpring().get_GradientPlasticity();

        double ua1=0.,ua2=0.;
        double ub1=0.,ub2=0.;
        double s1 ,s2 ;
        double z1 ,z2 ;
        if((elem.getINodeHier(1).getID() == nodeID && ((elem.getINodeHier(1).getp_posEFTable(elem.getID())==dof) || (elem.getINodeHier(1).getp_negEFTable(elem.getID())==dof)) )||
               (elem.getINodeHier(2).getID() == nodeID && ((elem.getINodeHier(2).getp_posEFTable(elem.getID())==dof) || (elem.getINodeHier(2).getp_negEFTable(elem.getID())==dof)) ) ){
            double[] xsis = new double[numGP];
            double[] weights = new double[numGP];
            switch(numGP){
                case 3:
                    xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                    weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                    break;
                case 4:
                    xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                    weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                    break;
                case 5:
                    xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                    weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                    break;
                default: System.exit(345); break;
            }


            for(int i=0; i<numGP; i++){
                ua1=0.;ua2=0.;
                ub1=0.;ub2=0.;
                s1=0.; s2=0.;
                z1=1.; z2=1.;
                val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                if(theInterface.IsPlastIncluded())s1 = elem.getINodeHier(1).getp(elem.getID())[wstep];
                if(theInterface.IsPlastIncluded())s2 = elem.getINodeHier(2).getp(elem.getID())[wstep];

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];

                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                        //s1=-s1; s2=-s2;
                    }
                }

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);


                if(elem.getINodeHier(1).getID() == nodeID && ((elem.getINodeHier(1).getp_posEFTable(elem.getID())==dof) || (elem.getINodeHier(1).getp_negEFTable(elem.getID())==dof)) ){
                    val*=kn*N1*(-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2)+(N1*s1+N2*s2))*(N1*z1+N2*z2)+
                            kh*(N1*s1+N2*s2)*N1*(N1*z1+N2*z2)+ge*(dN1*s1+dN2*s2)*N1*(N1*z1+N2*z2);
                }
                if(elem.getINodeHier(2).getID() == nodeID && ((elem.getINodeHier(2).getp_posEFTable(elem.getID())==dof) || (elem.getINodeHier(2).getp_negEFTable(elem.getID())==dof)) ){
                    val*=kn*N2*(-(N1*ua1+N2*ua2)+(N1*ub1+N2*ub2)+(N1*s1+N2*s2))*(N1*z1+N2*z2)+
                            kh*(N1*s1+N2*s2)*N2*(N1*z1+N2*z2)+ge*(dN1*s1+dN2*s2)*N2*(N1*z1+N2*z2);
                }

                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // END OF: methods for integrals on interfaces and interface elements
    /////////////////////////////////////////////////////////////////////////////////////////

    public void setGaussForEnergy(int n){
        switch(n){
            case 3: this.numGP=n; break;
            case 4: this.numGP=n; break;
            case 5: this.numGP=n; break;
            default: System.err.println(" Gauss Order for Integration of Energetics no 3 nor 4 neither 5, 5 selected by default"); this.numGP=5; break;
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // BEGIN OF: methods for integrals on interfaces and interface elements
    // 27 November 2011
    /////////////////////////////////////////////////////////////////////////////////////////
    
    public double IntegratePower(ELine elem, Domain theDomain, int step){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux=elem.getNodeHier(hierInelem_i).getv()[0][step][0];
                        uy=elem.getNodeHier(hierInelem_i).getv()[1][step][0];

                        px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
                        py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[1][step];

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }

    @Override
    double IntegrateM(Node colNode, Domain theDomain, Element elem_, Node elemNode) {
        double M=0.;
        ELine elem=(ELine) elem_;
        int n=theDomain.theFundSol.get_p_DOFs();

        int hierInelem=elem.getHierOfNode(elemNode.getID());
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];

        double xsi;

        for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){

                double val=1.0;
                CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
                weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

                xsi=CoordsOnElem[0];

                val*=elem.getShapeFunction(hierInelem, xsi);
                val*=elem.getJacobian(xsi);
                val*=weight[0];

                M+=val;
        }

        return M;
    }
    
    // **************************************************************************
    // NEW METHODS On 27/12/2011 Panagiotopoulos Christos
    // **************************************************************************
    public double IntegrateWork(ELine elem, Domain theDomain, int step, int state){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][state];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][state];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][state];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][state];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateWorkX(ELine elem, Domain theDomain, int step, int state){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];
                        
                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][state];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][state];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][state];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][state];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }
                        
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateWorkY(ELine elem, Domain theDomain, int step, int state){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];
                        
                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][state];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][state];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][state];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][state];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }
                        

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateWork(ELine elem, Domain theDomain, int DispStep, int TracStep, int DispState, int TracState){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        if(this.DomainAuxiliaryField){
                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][DispStep][DispState];
                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][DispStep][DispState];

                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[0][TracStep];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[1][TracStep];
                        }else{
                            ux=elem.getNodeHier(hierInelem_i).getu()[0][DispStep][DispState];
                            uy=elem.getNodeHier(hierInelem_i).getu()[1][DispStep][DispState];

                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[0][TracStep];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),TracState)[1][TracStep];
                        }

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateGeneralisedForce(ELine elem, Domain theDomain, int step, int nodeID, int dof, int state){
        double Length=0.;
        double ux,uy,px,py;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux=0.;
                        uy=0.;

                        if(theNode_i.getID() == nodeID){
                            switch(dof){
                                case 1: ux=1.; break;
                                case 2: uy=1.; break;
                            }
                        }

                        if(this.DomainAuxiliaryField){
                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }else{
                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
                        }
                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateNewCrack(IELine elem, int step){
        double Length=0.;
        double z,zp;
        if(step>0){
            for(Iterator<InterfaceNode>it=elem.getNodes().values().iterator(); it.hasNext();){
                InterfaceNode theNode = it.next();
                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;
                int hierInelem=elem.getHierOfINode(theNode.getID());
                z=theNode.getz(elem.getID())[step];
                zp=theNode.getz(elem.getID())[step-1];

                for(int i=0; i<SpaceQuadIntegrator.numofGauss; i++){

                        double val=1.0;
                        CoordsOnElem[0]=SpaceQuadIntegrator.theGaussData.getGaussCoordinate(i);
                        weight[0]=SpaceQuadIntegrator.theGaussData.getGaussWeight(i);

                        xsi=CoordsOnElem[0];


                        val*=elem.getShapeFunction(hierInelem, xsi)*(zp-z);
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=val;
                }
            }
        }
        
        return Length;
    }
    
    public double IntegrateDissipatedFrictionIEnergy(IELine elem, Interface theInterface, int wstep){
        // DRAFT
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double coulomb=theInterface.getEnergyDissipator().getCoulombCoef();

        double s_pos1 = 0.,s_pos2 = 0.;
        double s_neg1 = 0.,s_neg2 = 0.;
        double tn1 = 0., tn2 = 0.;
        double N1,N2;


        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }


        for(int i=0; i<numGP; i++){
            val=1.0;
            s_pos1 = 0. ; s_pos2 = 0.;
            s_neg1 = 0. ; s_neg2 = 0.;
            tn1=0.; tn2=0.;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            if(theInterface.IsFrictionIncluded()){
//                s_pos1 = elem.getINodeHier(1).getut()[wstep];
//                s_pos2 = elem.getINodeHier(2).getut()[wstep];
//                
//                s_neg1 = elem.getINodeHier(1).getut()[wstep-1];
//                s_neg2 = elem.getINodeHier(2).getut()[wstep-1];
            }
            
            if(wstep>0){
                tn1 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(1).getID(), wstep-1);
                tn2 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(2).getID(), wstep-1);
//                if(theInterface.getTangentSpring()==null){
//                    tn1 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(1).getID(), wstep-1);
//                    tn2 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(2).getID(), wstep-1);
//                }else{
//                    if(theInterface.IsDamageIncluded()){
//                        if(Math.abs(elem.getINodeHier(1).getz(elem.getID())[wstep])<=1.e-3)tn1 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(1).getID(), wstep-1);
//                        if(Math.abs(elem.getINodeHier(2).getz(elem.getID())[wstep])<=1.e-3)tn2 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(2).getID(), wstep-1);
//                    }
//                }
            }
            if(tn1>0)tn1=0; if(tn2>0)tn2=0.;
            
            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);

            val*=coulomb*(N1*tn1+N2*tn2);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }
    
    public double IntegrateDissipatedFrictionIForce(IELine elem, Interface theInterface, int wstep, int NodeID){
        // DRAFT
        double dis=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;

        double coulomb=theInterface.getEnergyDissipator().getCoulombCoef();

        double s_pos1 = 0.,s_pos2 = 0.;
        double s_neg1 = 0.,s_neg2 = 0.;
        double tn1 = 0., tn2 = 0.;
        double N1,N2;


        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }


        for(int i=0; i<numGP; i++){
            val=1.0;
            s_pos1 = 0. ; s_pos2 = 0.;
            s_neg1 = 0. ; s_neg2 = 0.;
            tn1=0.; tn2=0.;
            CoordsOnElem[0]=xsis[i];
            weight[0]=weights[i];
            xsi=CoordsOnElem[0];
            
            if(theInterface.IsFrictionIncluded()){
//                s_pos1 = elem.getINodeHier(1).getut_pos()[wstep];
//                s_pos2 = elem.getINodeHier(2).getut_pos()[wstep];
//                
//                s_neg1 = elem.getINodeHier(1).getut_neg()[wstep];
//                s_neg2 = elem.getINodeHier(2).getut_neg()[wstep];
            }
            
            if(wstep>0){
                tn1 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(1).getID(), wstep-1);
                tn2 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(2).getID(), wstep-1);
//                if(theInterface.getTangentSpring()==null){
//                    tn1 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(1).getID(), wstep-1);
//                    tn2 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(2).getID(), wstep-1);
//                }else{
//                    if(theInterface.IsDamageIncluded()){
//                        if(Math.abs(elem.getINodeHier(1).getz(elem.getID())[wstep])<=1.e-3)tn1 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(1).getID(), wstep-1);
//                        if(Math.abs(elem.getINodeHier(2).getz(elem.getID())[wstep])<=1.e-3)tn2 = theInterface.getNormalBEMTraction(elem.getID(), elem.getINodeHier(2).getID(), wstep-1);
//                    }
//                }
            }
            if(tn1>0)tn1=0; if(tn2>0)tn2=0.;
            
            N1=elem.getShapeFunction(1, xsi);
            N2=elem.getShapeFunction(2, xsi);
            
            if(NodeID==elem.getINodeHier(1).getID()){
                val*=coulomb*(N1*tn1+N2*tn2);
            }else{
                val*=coulomb*(N1*tn1+N2*tn2);
            }
            
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            dis+=val;
        }

        return dis;
    }
    
//    public double IntegrateEnergy(ELine elem, Domain theDomain, int step, int state, double nu, double tau){
//        double Length=0.;
//        double ux,uy,px,py;
//        Node theNode_i;
//        Node theNode_j;
//        int pstep;
//
//        double[] xsis = new double[numGP];
//        double[] weights = new double[numGP];
//        switch(numGP){
//            case 3:
//                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
//                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
//                break;
//            case 4:
//                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
//                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
//                break;
//            case 5:
//                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
//                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
//                break;
//           default: System.exit(345); break;
//        }
//
//        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
//            theNode_i= it.next();
//            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
//                theNode_j= jt.next();
//
//                double[] CoordsOnElem=new double[1];
//                double[] weight = new double[1];
//                double xsi;
//
//                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());
//
//                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());
//
//                for(int i=0; i<numGP; i++){
//
//                        double val=1.0;
//                        CoordsOnElem[0]=xsis[i];
//                        weight[0]=weights[i];
//
//                        xsi=CoordsOnElem[0];
//                        pstep=step-1;if(step==0)pstep=0;
//                        if(this.DomainAuxiliaryField){
//                            ux=elem.getNodeHier(hierInelem_i).getu_aux()[0][step][state];
//                            uy=elem.getNodeHier(hierInelem_i).getu_aux()[1][step][state];
//                            
//                            ux=(ux-elem.getNodeHier(hierInelem_i).getu_aux()[0][pstep][state])/tau;
//                            uy=(uy-elem.getNodeHier(hierInelem_i).getu_aux()[1][pstep][state])/tau;
//
//                            px=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
//                            py=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
//                        }else{
//                            ux=elem.getNodeHier(hierInelem_i).getu()[0][step][state];
//                            uy=elem.getNodeHier(hierInelem_i).getu()[1][step][state];
//                            
//                            ux=(ux-elem.getNodeHier(hierInelem_i).getu()[0][pstep][state])/tau;
//                            uy=(uy-elem.getNodeHier(hierInelem_i).getu()[1][pstep][state])/tau;
//
//                            px=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
//                            py=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[1][step];
//                        }
//                        
//                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*(ux*px+uy*py);
//                        val*=elem.getJacobian(xsi);
//                        val*=weight[0];
//
//                        Length+=val;
//                }
//            }
//        }
//        //System.err.println("2* Strain Energy Computed = "+Length);
//        return Length;
//    }

    public double IntegrateTIPower(ELine elem, Domain theDomain, int wstep){
        double Length=0.;
        double ux1,uy1,px1,py1;
        double ux2,uy2,px2,py2;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux1=elem.getNodeHier(hierInelem_i).getu()[0][wstep-1][3];
                        uy1=elem.getNodeHier(hierInelem_i).getu()[1][wstep-1][3];

                        px1=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[0][wstep-1];
                        py1=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[1][wstep-1];
                        
                        ux2=elem.getNodeHier(hierInelem_i).getu()[0][wstep][3];
                        uy2=elem.getNodeHier(hierInelem_i).getu()[1][wstep][3];

                        px2=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[0][wstep];
                        py2=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[1][wstep];

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*( (px2*ux2+px1*ux2-px2*ux1-px1*ux1)+(py2*uy2+py1*uy2-py2*uy1-py1*uy1));
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=0.5*val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    public double IntegrateTIPowerDifference(ELine elem, Domain theDomain, int wstep){
        // int_{t} (p-paux)*v dh
        double Length=0.;
        double ux1,uy1,px1,py1;
        double ux2,uy2,px2,py2;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
           default: System.exit(345); break;
        }

        for(Iterator<Node>it=elem.getNodes().values().iterator(); it.hasNext();){
            theNode_i= it.next();
            for(Iterator<Node>jt=elem.getNodes().values().iterator(); jt.hasNext();){
                theNode_j= jt.next();

                double[] CoordsOnElem=new double[1];
                double[] weight = new double[1];
                double xsi;

                int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

                int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

                for(int i=0; i<numGP; i++){

                        double val=1.0;
                        CoordsOnElem[0]=xsis[i];
                        weight[0]=weights[i];

                        xsi=CoordsOnElem[0];

                        ux1=elem.getNodeHier(hierInelem_i).getu()[0][wstep-1][3];
                        uy1=elem.getNodeHier(hierInelem_i).getu()[1][wstep-1][3];
                        
                        px1=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[0][wstep-1];
                        py1=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[1][wstep-1];
                        

                        px1-=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[0][wstep-1];
                        py1-=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[1][wstep-1];
                        
                        ux2=elem.getNodeHier(hierInelem_i).getu()[0][wstep][3];
                        uy2=elem.getNodeHier(hierInelem_i).getu()[1][wstep][3];

                        px2=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[0][wstep];
                        py2=elem.getNodeHier(hierInelem_j).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[1][wstep];
                        
                        px2-=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[0][wstep];
                        py2-=elem.getNodeHier(hierInelem_j).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),3)[1][wstep];

                        val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi)*( (px2*ux2+px1*ux2-px2*ux1-px1*ux1)+(py2*uy2+py1*uy2-py2*uy1-py1*uy1));
                        val*=elem.getJacobian(xsi);
                        val*=weight[0];

                        Length+=0.5*val;
                }
            }
        }
        //System.err.println("2* Strain Energy Computed = "+Length);
        return Length;
    }
    
    // for energy balance of viscoelasticity
    public double IntegrateTangentialStoredIEnergyAUX(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kt = 0.0;
        double k0 = 0.0; double r=1.;
        double a0 = 0.0;
        ViscousMaterial vMaterial_main;
        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsTangentialIncluded() && theInterface.getTangentSpring()!=null){
            kt = theInterface.getTangentSpring().getStiffness();
            k0 = theInterface.getTangentSpring().get_GradientZ();
            if(theInterface.getEnergyDissipator()!=null)a0 = theInterface.getEnergyDissipator().getAlpha0();
            if(Math.abs(k0)>1.e-10)r = theInterface.getTangentSpring().get_r();


            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double va1,va2;
            double vb1 = 0,vb2 = 0;
            double N1,N2;
            double dN1,dN2;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getut()[wstep];
                ua2 = elem.getINodeHier(2).getut()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getut()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ua2 = elem.getINodeHier(2).getut()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getut()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getut()[wstep];
                        ub2 = elem.getINodeHier(2).getut()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getut()[wstep];
                    }
                }
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getut_(vMaterial_main,tau)[wstep];
                va2 = elem.getINodeHier(2).getut_(vMaterial_main,tau)[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}

                val*=0.5*kt*((N1*va1+N2*va2)-(N1*vb1+N2*vb2))*((N1*va1+N2*va2)-(N1*vb1+N2*vb2))*(N1*z1+N2*z2)
                        +k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r
                        -a0*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }
    
    public double IntegrateNormalStoredIEnergyAUX(IELine elem, Interface theInterface, int wstep, double tau){
        // IT IS DEVELOPED FOR 2 NODE LINE INTERFACE ELEMENT
        double stored=0.;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;
        double val;
        double kn = 0.0;
        double k0 = 0.0; double r=1.;
        double a0 = 0.0;
        ViscousMaterial vMaterial_main;
        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        if(theInterface.IsNormalIncluded() && theInterface.getNormalSpring()!=null){
            kn = theInterface.getNormalSpring().getStiffness();
            k0 = theInterface.getNormalSpring().get_GradientZ();
            if(theInterface.getEnergyDissipator()!=null)a0 = theInterface.getEnergyDissipator().getAlpha0();
            if(Math.abs(k0)>1.e-10)r = theInterface.getNormalSpring().get_r();


            double z1 = 1.,z2 = 1.;
            double ua1,ua2;
            double ub1,ub2;
            double va1,va2;
            double vb1 = 0,vb2 = 0;
            double N1,N2;
            double dN1,dN2;

            for(int i=0; i<numGP; i++){
                val=1.0;
                ub1=0.;ub2=0.;
                z1 = 1.;z2 = 1.;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];
                xsi=CoordsOnElem[0];

                ua1 = elem.getINodeHier(1).getun()[wstep];
                ua2 = elem.getINodeHier(2).getun()[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    if(elem.getINodeHier(1).isMain()||elem.getINodeHier(2).isMain()){
                        ua1 = elem.getINodeHier(1).getun()[wstep];
                        ub1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ua2 = elem.getINodeHier(2).getun()[wstep];
                        ub2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }else{
                        ub1 = elem.getINodeHier(1).getun()[wstep];
                        ua1 = elem.getINodeHier(1).getTwin().getun()[wstep];
                        ub2 = elem.getINodeHier(2).getun()[wstep];
                        ua2 = elem.getINodeHier(2).getTwin().getun()[wstep];
                    }
                }
                vMaterial_main=(ViscousMaterial) theInterface.getDomain1().getMaterial();
                va1 = elem.getINodeHier(1).getun_(vMaterial_main,tau)[wstep];
                va2 = elem.getINodeHier(2).getun_(vMaterial_main,tau)[wstep];
                if(theInterface.getNumOfConnectedDomains()==2){
                    throw new UnsupportedOperationException("Not supported yet.");
                }

                if(theInterface.IsDamageIncluded())z1 = elem.getINodeHier(1).getz(elem.getID())[wstep];
                if(theInterface.IsDamageIncluded())z2 = elem.getINodeHier(2).getz(elem.getID())[wstep];

                N1=elem.getShapeFunction(1, xsi);
                N2=elem.getShapeFunction(2, xsi);

                dN1=elem.getShapeFunction_xsi(1, xsi);
                dN2=elem.getShapeFunction_xsi(2, xsi);
                //if(theInterface.zconstant){dN1=0; dN2=0;}

                val*=0.5*kn*((N1*va1+N2*va2)-(N1*vb1+N2*vb2))*((N1*va1+N2*va2)-(N1*vb1+N2*vb2))*(N1*z1+N2*z2)
                        +k0*Math.pow(Math.abs(dN1*z1+dN2*z2), r)/r
                        -a0*(N1*z1+N2*z2);
                val*=elem.getJacobian(xsi);
                val*=weight[0];
                stored+=val;
            }
        }
        return stored;
    }

    public double IntegrateInternalPointDisp(ELine elem, Node elemNode, Domain theDomain, ResultPoint theInternalPoint, int wdisp, int step, int wstate) {
        double disp=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){
            
            double val=1.0;
            CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
            weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

            xsi=CoordsOnElem[0];

            R=elem.getDistance(theInternalPoint, xsi);
            Normal=elem.getNormal(xsi);
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Normal);

            val*=elem.getShapeFunction(hierInelem, xsi);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            if(theDomain.getFundamentalSolution().ndofs==1){
            val*=theDomain.getFundamentalSolution().get_u_fund().get(wdisp, 0)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[0][step])
                    -theDomain.getFundamentalSolution().get_p_fund().get(wdisp, 0)*elem.getNodeHier(hierInelem).getu()[0][step][wstate]
                    ;    
            }else{
            val*=theDomain.getFundamentalSolution().get_u_fund().get(wdisp, 0)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[0][step]+((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[0])
                    +theDomain.getFundamentalSolution().get_u_fund().get(wdisp, 1)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[1][step]+((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[1])
                    -theDomain.getFundamentalSolution().get_p_fund().get(wdisp, 0)*elem.getNodeHier(hierInelem).getu()[0][step][wstate]
                    -theDomain.getFundamentalSolution().get_p_fund().get(wdisp, 1)*elem.getNodeHier(hierInelem).getu()[1][step][wstate]
                    ;
            }
            disp+=val;
        }
        return disp;
    }
    
    public double IntegrateInternalPointStress(ELine elem, Node elemNode, Domain theDomain, ResultPoint theInternalPoint, int wstress, int step, int wstate) {
        double disp=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){
            
            double val=1.0;
            CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
            weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

            xsi=CoordsOnElem[0];

            R=elem.getDistance(theInternalPoint, xsi);
            Normal=elem.getNormal(xsi);
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Normal);

            val*=elem.getShapeFunction(hierInelem, xsi);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            val*=    theDomain.getFundamentalSolution().get_s_fund().get(0, wstress)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[0][step]+
                                                                                        ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[0])
                    +theDomain.getFundamentalSolution().get_s_fund().get(1, wstress)*(elem.getNodeHier(hierInelem).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[1][step]+
                                                                                        ((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange()*Normal[1])
                    -theDomain.getFundamentalSolution().get_r_fund().get(0, wstress)*elem.getNodeHier(hierInelem).getu()[0][step][wstate]
                    -theDomain.getFundamentalSolution().get_r_fund().get(1, wstress)*elem.getNodeHier(hierInelem).getu()[1][step][wstate]
                    ;
            disp+=val;
        }
        return disp;
    }
    
    public double IntegrateInternalPointStressAux(ELine elem, Node elemNode, Domain theDomain, ResultPoint theInternalPoint, int wstress, int step, int wstate) {
        double disp=0;
        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double[] R;
        double[] Normal;
        double xsi,eta;
        int hierInelem=elem.getHierOfNode(elemNode.getID());
        for(int i=0; i<SpaceLineIntegrator.numofGauss; i++){
            
            double val=1.0;
            CoordsOnElem[0]=SpaceLineIntegrator.theGaussData.getGaussCoordinate(i);
            weight[0]=SpaceLineIntegrator.theGaussData.getGaussWeight(i);

            xsi=CoordsOnElem[0];

            R=elem.getDistance(theInternalPoint, xsi);
            Normal=elem.getNormal(xsi);
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Normal);

            val*=elem.getShapeFunction(hierInelem, xsi);
            val*=elem.getJacobian(xsi);
            val*=weight[0];
            val*=    theDomain.getFundamentalSolution().get_s_fund().get(0, wstress)*elem.getNodeHier(hierInelem).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[0][step]
                    +theDomain.getFundamentalSolution().get_s_fund().get(1, wstress)*elem.getNodeHier(hierInelem).getp_aux(elem,theDomain.getFundamentalSolution().get_p_DOFs(),wstate)[1][step]
                    -theDomain.getFundamentalSolution().get_r_fund().get(0, wstress)*elem.getNodeHier(hierInelem).getu_aux()[0][step][wstate]
                    -theDomain.getFundamentalSolution().get_r_fund().get(1, wstress)*elem.getNodeHier(hierInelem).getu_aux()[1][step][wstate]
                    ;
            disp+=val;
        }
        return disp;
    }
    
    public double IntegrateShapeFunctionProduct(ELine elem, int dof, int node_i, int node_j){
        double Length=0.;
        Node theNode_i;
        Node theNode_j;

        double[] xsis = new double[numGP];
        double[] weights = new double[numGP];
        switch(numGP){
            case 3:
                xsis[0]=0.774596669; xsis[1]=0.0;  xsis[2]=-xsis[0];
                weights[0]=0.555555555; weights[1]=0.888888888;  weights[2]=weights[0];
                break;
            case 4:
                xsis[0]=0.861136311; xsis[1]=0.339981043;  xsis[2]=-xsis[1];  xsis[3]=-xsis[0];
                weights[0]=0.347854845; weights[1]=0.652145154; weights[2]=weights[1];  weights[3]=weights[0];
                break;
            case 5:
                xsis[0]=0.906179845; xsis[1]=0.538469310;  xsis[2]=0.00; xsis[3]=-xsis[1];  xsis[4]=-xsis[0];
                weights[0]=0.236926885; weights[1]=0.478628670;  weights[2]=0.568888888; weights[3]=weights[1];  weights[4]=weights[0];
                break;
            default: System.exit(345); break;
        }

        
        theNode_i= elem.getNode(node_i);
        theNode_j= elem.getNode(node_j);

        double[] CoordsOnElem=new double[1];
        double[] weight = new double[1];
        double xsi;

        int hierInelem_i=elem.getHierOfNode(theNode_i.getID());

        int hierInelem_j=elem.getHierOfNode(theNode_j.getID());

        for(int i=0; i<numGP; i++){

                double val=1.0;
                CoordsOnElem[0]=xsis[i];
                weight[0]=weights[i];

                xsi=CoordsOnElem[0];

                val*=elem.getShapeFunctionProduct(hierInelem_i, hierInelem_j, xsi);
                val*=elem.getJacobian(xsi);
                val*=weight[0];

                Length+=val;
        }
        
        return Length;
    }
    
}


//IntegrateNormalStoredIEnergyLinear
//IntegrateNormalStoredIEnergyQuad
//IntegrateNormalStoredIEnergy
//
//IntegrateTangentialStoredIEnergy
//IntegrateTangentialStoredIEnergyQuad
//IntegrateTangentialStoredIEnergyLinear
//
//IntegrateIntGeneralisedNormalForce
//IntegrateIntGeneralisedNormalForceQuad
//IntegrateIntGeneralisedNormalForceLinear
//
//IntegrateIntGeneralisedTangentForce
//IntegrateIntGeneralisedTangentForceQuad
//IntegrateIntGeneralisedTangentForceLinear
//
//IntegrateIntGeneralisedNormalDrivingForce
//IntegrateIntGeneralisedTangentDrivingForce
//
//IntegrateTangentialStoredIEnergyAUX
//IntegrateNormalStoredIEnergyAUX