/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class Transformation_001 extends TransformationPlane{
 /**
 *
 * transformation given at 
 * "Programming the Boundary Element Method"
 * An introduction for Engineers
 * Gernot Beer
 * p.145
 */   
    // constructor
    public Transformation_001(int[] theNodes){
        this.theNodes=theNodes;
    }

    @Override
    public double getJacobian(double xsi, double eta) {
        double jac=0.;
        double j11,j12,j21,j22; 
        j11=0.; j12=0.; j21=0.; j22=0.;
        int nd; double xsi_L=0.;  double eta_L=0.;
        for(int i=1; i<=3; i++){
            nd=theNodes[i-1];
            switch(nd){
                case 1: xsi_L=-1.   ;   eta_L=-1.   ;  break;
                case 2: xsi_L=1.    ;   eta_L=-1.   ;  break;
                case 3: xsi_L=1.    ;   eta_L=1.    ;  break;
                case 4: xsi_L=-1.   ;   eta_L=1.    ;  break;
                case 5: xsi_L=0.    ;   eta_L=-1.   ;  break;
                case 6: xsi_L=1.    ;   eta_L=0.    ;  break;
                case 7: xsi_L=0.    ;   eta_L=1.    ;  break;
                case 8: xsi_L=-1.   ;   eta_L=0.    ;  break;
                case 9: xsi_L=0.    ;   eta_L=0.    ;  break;
            }
            j11+=getShapeFunction_xsi(i, xsi, eta)*xsi_L;
            j12+=getShapeFunction_xsi(i, xsi, eta)*eta_L;
            j21+=getShapeFunction_eta(i, xsi, eta)*xsi_L;
            j22+=getShapeFunction_eta(i, xsi, eta)*eta_L;
        }
        jac=j11*j22-j12*j21;
        return jac;
    }

    @Override
    protected double getShapeFunction(int wSF, double xsi, double eta) {
        double N;
        switch(wSF){
            case 1: N=0.25*(1+xsi)*(1-eta); break;
            case 2: N=0.25*(1+xsi)*(1+eta); break;
            case 3: N=0.50*(1-xsi); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    protected double getShapeFunction_xsi(int wSF, double xsi, double eta) {
        double N;
        switch(wSF){
            case 1: N=0.25*(1.-eta); break;
            case 2: N=0.25*(1.+eta); break;
            case 3: N=0.50*(-1.); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    protected double getShapeFunction_eta(int wSF, double xsi, double eta) {
        double N;
        switch(wSF){
            case 1: N=0.25*(-1.-xsi); break;
            case 2: N=0.25*(1.+xsi); break;
            case 3: N=0.50*(0.); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    public double getXsi(double xsi_, double eta_) {
        double xsi=0.;
        double xsi_L=0.;
        int nd;
        for(int i=1; i<=3; i++){
            nd=theNodes[i-1];
            switch(nd){
                case 1: xsi_L=-1.; break;
                case 2: xsi_L=1.; break;
                case 3: xsi_L=1.; break;
                case 4: xsi_L=-1.; break;
                case 5: xsi_L=0.; break;
                case 6: xsi_L=1.; break;
                case 7: xsi_L=0.; break;
                case 8: xsi_L=-1.; break;
                case 9: xsi_L=0.; break;
            }
            xsi+=getShapeFunction(i, xsi_, eta_)*xsi_L;
        }
        return xsi;
    }

    @Override
    public double getEta(double xsi_, double eta_) {
        double eta=0.;
        double eta_L=0.;
        int nd;
        for(int i=1; i<=3; i++){
            nd=theNodes[i-1];
            switch(nd){
                case 1: eta_L=-1.; break;
                case 2: eta_L=-1.; break;
                case 3: eta_L=1.; break;
                case 4: eta_L=1.; break;
                case 5: eta_L=-1.; break;
                case 6: eta_L=0.; break;
                case 7: eta_L=1.; break;
                case 8: eta_L=0.; break;
                case 9: eta_L=0.; break;
            }
            eta+=getShapeFunction(i, xsi_, eta_)*eta_L;
        }
        return eta;
    }

}
