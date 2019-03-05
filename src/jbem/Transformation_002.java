/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class Transformation_002 extends TransformationPlane{
 /**
 *
 * transformation given for mesh reasons 
 * to increase the accuracy of intgrations over quad elements
 */   
    // constructor
    public Transformation_002(int[] theNodes){
        this.theNodes=theNodes;
    }

    @Override
    public double getJacobian(double xsi, double eta) {
        double jac=0.;
        jac=1./(this.theNodes[2]*this.theNodes[3]);
        return jac;
    }

    @Override
    protected double getShapeFunction(int wSF, double xsi, double eta) {
        double N;
        switch(wSF){
            case 1: N=0.25*(1-xsi)*(1-eta); break;
            case 2: N=0.25*(1+xsi)*(1-eta); break;
            case 3: N=0.25*(1+xsi)*(1+eta); break;
            case 4: N=0.25*(1-xsi)*(1+eta); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    protected double getShapeFunction_xsi(int wSF, double xsi, double eta) {
        double N;
        switch(wSF){
            case 1: N=-0.25*(1-eta); break;
            case 2: N= 0.25*(1-eta); break;
            case 3: N= 0.25*(1+eta); break;
            case 4: N=-0.25*(1+eta); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    protected double getShapeFunction_eta(int wSF, double xsi, double eta) {
        double N;
        switch(wSF){
            case 1: N=-0.25*(1-xsi); break;
            case 2: N=-0.25*(1+xsi); break;
            case 3: N= 0.25*(1+xsi); break;
            case 4: N= 0.25*(1-xsi); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    public double getXsi(double xsi_, double eta_) {
        double xsi=0.;
        double[] c1,c2,c3,c4;
        c1=new double[2];c2=new double[2];c3=new double[2];c4=new double[2];
        
        c1[0]=-1.+2./this.theNodes[2]*(this.theNodes[0]-1);
        c1[1]=-1.+2./this.theNodes[3]*(this.theNodes[1]-1);
        
        c2[0]=c1[0]+2./this.theNodes[2];
        c2[1]=c1[1];
        
        c3[0]=c2[0];
        c3[1]=c1[1]+2./this.theNodes[3];
        
        c4[0]=c1[0];
        c4[1]=c3[1];
        xsi=getShapeFunction(1, xsi_, eta_)*c1[0]+
                getShapeFunction(2, xsi_, eta_)*c2[0]+
                getShapeFunction(3, xsi_, eta_)*c3[0]+
                getShapeFunction(4, xsi_, eta_)*c4[0];
        return xsi;
    }

    @Override
    public double getEta(double xsi_, double eta_) {
        double eta=0.;
        double[] c1,c2,c3,c4;
        c1=new double[2];c2=new double[2];c3=new double[2];c4=new double[2];
        
        c1[0]=-1.+2./this.theNodes[2]*(this.theNodes[0]-1);
        c1[1]=-1.+2./this.theNodes[3]*(this.theNodes[1]-1);
        
        c2[0]=c1[0]+2./this.theNodes[2];
        c2[1]=c1[1];
        
        c3[0]=c2[0];
        c3[1]=c1[1]+2./this.theNodes[3];
        
        c4[0]=c1[0];
        c4[1]=c3[1];
        eta=getShapeFunction(1, xsi_, eta_)*c1[1]+
                getShapeFunction(2, xsi_, eta_)*c2[1]+
                getShapeFunction(3, xsi_, eta_)*c3[1]+
                getShapeFunction(4, xsi_, eta_)*c4[1];
        return eta;
    }

}
