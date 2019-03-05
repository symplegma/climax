/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
/**
 *
 * @author pchr
 */
public class SpaceNodeIntegrator extends SpaceIntegrator{
    
    // constructor
    public SpaceNodeIntegrator(){
    }

    @Override
    AbstractMatrix IntegrateDu(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        AbstractMatrix H;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        double[] Out=new double[1];
        if(colNode.getID()!=elemNode.getID()){
            double[] R=new double[1];
            Domain1D aDomain = (Domain1D) theDomain;
            //R[0]=aDomain.getL();
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                R[0]=aDomain.getL();
                Out[0]=1.;
            }else{
                R[0]=-aDomain.getL();
                Out[0]=-1.;
            }
            FundamentalSolution.theFSdata.setR(R);
            H=theDomain.getFundamentalSolution().get_u_fund().times(Out[0]);
            //H=theDomain.getFundamentalSolution().get_u_fund().times(Out[0]);
        }else{
            double[] R=new double[1];
            Domain1D aDomain = (Domain1D) theDomain;
            //R[0]=aDomain.getL();
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                Out[0]=1.;
            }else{
                Out[0]=-1.;
            }
            R[0]=0.;
            FundamentalSolution.theFSdata.setR(R);
            //H=theDomain.getFundamentalSolution().get_u_fund();
            H=theDomain.getFundamentalSolution().get_u_fund().times(Out[0]);
        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDp(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        AbstractMatrix H;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        double[] R=new double[1];
        double[] Out=new double[1];
        if(colNode.getID()!=elemNode.getID()){
            Domain1D aDomain = (Domain1D) theDomain;
            //R[0]=aDomain.getL();
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                R[0]=aDomain.getL();
                Out[0]=1.;
            }else{
                R[0]=-aDomain.getL();
                Out[0]=-1.;
            }
            //Out[0]=1.;
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Out);
            H=theDomain.getFundamentalSolution().get_p_fund().times(Out[0]);
            //H=theDomain.getFundamentalSolution().get_p_fund().times(Out[0]);
        }else{
            Domain1D aDomain = (Domain1D) theDomain;
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                R[0]=0.;
                Out[0]=1.;
            }else{
                R[0]=0.;
                Out[0]=-1.;
            }
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Out);
            
            int m=FundamentalSolution.theFSdata.getMstep();
            int nn=FundamentalSolution.theFSdata.getNstep();
            if((m-nn)>1){
                H=theDomain.getFundamentalSolution().get_p_fund().times(Out[0]);
            }
        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDv(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        AbstractMatrix H;
        int n=theDomain.theFundSol.get_v_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        TimeFundamentalSolution tFunf;
        tFunf=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
        double[] Out=new double[1];
        if(colNode.getID()!=elemNode.getID()){
            double[] R=new double[1];
            Domain1D aDomain = (Domain1D) theDomain;
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                R[0]=aDomain.getL();
                Out[0]=1.;
            }else{
                R[0]=-aDomain.getL();
                Out[0]=-1.;
            }
            FundamentalSolution.theFSdata.setR(R);
            H=tFunf.get_v_fund().times(Out[0]);
            //H=tFunf.get_v_fund().times(Out[0]);
        }else{
            double[] R=new double[1];
            Domain1D aDomain = (Domain1D) theDomain;
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                Out[0]=1.;
            }else{
                Out[0]=-1.;
            }
            R[0]=0.;
            FundamentalSolution.theFSdata.setR(R);
            H=tFunf.get_v_fund().times(Out[0]);
            //H=tFunf.get_v_fund().times(Out[0]);
        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDp_res(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        AbstractMatrix H;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        if(colNode.getID()!=elemNode.getID()){
            double[] Out=new double[1];
            double[] R=new double[1];
            Domain1D aDomain = (Domain1D) theDomain;
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                R[0]=aDomain.getL();
                Out[0]=-1.;
            }else{
                R[0]=-aDomain.getL();
                //Out[0]=-1.;
                Out[0]=1.;
            }
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Out);
            TimeFundamentalSolution tFund;
            tFund=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
            H=tFund.get_pres_fund().times(Out[0]);
        }else{
        }
        return H;
    }

    @Override
    AbstractMatrix IntegrateDp_dif(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        AbstractMatrix H;
        int n=theDomain.theFundSol.get_u_DOFs();
        H=new AbstractMatrix(n,n); H.init();
        double[] R=new double[1];
        double[] Out=new double[1];
        if(colNode.getID()!=elemNode.getID()){
            
        }else{
            Domain1D aDomain = (Domain1D) theDomain;
            if(aDomain.getNodeHier(colNode.getID()) == 1){
                R[0]=0.;
                Out[0]=1.;
            }else{
                R[0]=0.;
                Out[0]=-1.;
            }
            FundamentalSolution.theFSdata.setR(R);
            FundamentalSolution.theFSdata.setOutwardNormal(Out);
            
            int m=FundamentalSolution.theFSdata.getMstep();
            int nn=FundamentalSolution.theFSdata.getNstep();
            if((m-nn)>1){
                TimeFundamentalSolution tFunf;
                tFunf=(TimeFundamentalSolution) theDomain.getFundamentalSolution();
                H=tFunf.get_pdif_fund().times(Out[0]);
            }
        }
        return H;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Energy, Power, etc..
    // 27 November 2011
    /////////////////////////////////////////////////////////////////////////////////////////
    public double IntegrateWork(ENode elem, Domain theDomain, int step){
        return elem.getNodeHier(1).getu()[0][step][0]*elem.getNodeHier(1).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
    }
    
    public double IntegrateWork(ENode elem, Domain theDomain, int step, int state){
        return elem.getNodeHier(1).getu()[0][step][state]*elem.getNodeHier(1).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs(),state)[0][step];
    }
    
    public double IntegratePower(ENode elem, Domain theDomain, int step){
        return elem.getNodeHier(1).getv()[0][step][0]*elem.getNodeHier(1).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
    }
    
    public double IntegrateTIPower(ENode elem, Domain theDomain, int step, double dt){
        double energy=0.0;
        double p1,p2,v1,v2;
        if(step>0){
            for(int i=1; i<=step; i++){
                p1=elem.getNodeHier(1).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][i-1];
                p2=elem.getNodeHier(1).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][i];
                v1=elem.getNodeHier(1).getv()[0][i-1][0];
                v2=elem.getNodeHier(1).getv()[0][i][0];
                energy+=(p1*(2.*v1+v2)+p2*(v1+2.*v2))*dt/6.;
            }
        }else{
            p2=elem.getNodeHier(1).getp(elem,theDomain.getFundamentalSolution().get_p_DOFs())[0][step];
            v2=elem.getNodeHier(1).getv()[0][step][0];
            energy+=p2*v2;
        }
        
        
        return energy;
    }

    @Override
    double IntegrateM(Node colNode, Domain theDomain, Element elem, Node elemNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
