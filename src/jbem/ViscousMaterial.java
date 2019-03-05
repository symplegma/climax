/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

/**
 *
 * @author pchr
 */
public class ViscousMaterial  extends ElasticMat{
    double DispRelaxationTime_0=1.0;
    double DispRelaxationTime_1=0.0;
    double DispRelaxationTime_2=0.0;
    double StressRelaxationTime_0=1.0;
    double StressRelaxationTime_1=0.0;
    double StressRelaxationTime_2=0.0;
    
    public ViscousMaterial(int id){
        this.id=id;
    }
    
    public void setDispRelaxationTime_0(double val){DispRelaxationTime_0=val;}
    
    public void setDispRelaxationTime_1(double val){DispRelaxationTime_1=val;}
    
    public void setDispRelaxationTime_2(double val){DispRelaxationTime_2=val;}
    
    public void setStressRelaxationTime_0(double val){StressRelaxationTime_0=val;}
    
    public void setStressRelaxationTime_1(double val){StressRelaxationTime_1=val;}
    
    public void setStressRelaxationTime_2(double val){StressRelaxationTime_2=val;}
    
    public double getDispRelaxationTime_0(){return this.DispRelaxationTime_0;}
    
    public double getDispRelaxationTime_1(){return this.DispRelaxationTime_1;}
    
    public double getDispRelaxationTime_2(){return this.DispRelaxationTime_2;}
    
    public double getStressRelaxationTime_0(){return this.StressRelaxationTime_0;}
    
    public double getStressRelaxationTime_1(){return this.StressRelaxationTime_1;}
    
    public double getStressRelaxationTime_2(){return this.StressRelaxationTime_2;}
    
    public double[] TransformBoundaryDisplacement(double[] disp, double tau) {
        double[] auxDisp;
        double x0=DispRelaxationTime_0;
        double x1=DispRelaxationTime_1;
        double x2=DispRelaxationTime_2;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        auxDisp = new double[disp.length];
        auxDisp[0] = disp[0];
        for (int it = 1; it < disp.length; it++) {
            if(it>1){
                auxDisp[it] = ((x2+tau*x1+tau*tau*x0)/(y2+tau*y1+tau*tau*y0))*disp[it]-
                        ((2*x2+tau*x1)/(y2+tau*y1+tau*tau*y0))*disp[it-1]+
                        ((x2)/(y2+tau*y1+tau*tau*y0))*disp[it-2];
            }else if(it>0){
                auxDisp[it] = ((x2+tau*x1+tau*tau*x0)/(y2+tau*y1+tau*tau*y0))*disp[it]-
                        ((2*x2+tau*x1)/(y2+tau*y1+tau*tau*y0))*disp[it-1];
            }else{
                auxDisp[it] = ((x2+tau*x1+tau*tau)/(y2+tau*y1+tau*tau*y0))*disp[it];
            }
            
        }
        return auxDisp;
    }
    
    public double[] TransformBoundaryTraction(double[] trac, double tau) {
        double[] auxTrac;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        auxTrac = new double[trac.length];
        for (int it = 0; it < trac.length; it++) {
            
            if(it>1){
                auxTrac[it] = trac[it]-((2*y2+y1*tau)/(y2+y1*tau+y0*tau*tau))*trac[it-1]+
                        ((y2)/(y2+y1*tau+y0*tau*tau))*trac[it-2];
            }else if(it>0){
                auxTrac[it] = trac[it]-((2*y2+y1*tau)/(y2+y1*tau+y0*tau*tau))*trac[it-1];
            }else{
                auxTrac[it] = trac[it];
            }
        }
        return auxTrac;
    }
    
    /*
     * Assumed that v(k)=Phi_1(tau)*u(k)+Phi_2(tau)*u(k-1)+Phi_3(tau)*u(k-2)
     * this function returns Phi_1(tau)
     */
    public double getPhi_1(double tau){
        double x0=DispRelaxationTime_0;
        double x1=DispRelaxationTime_1;
        double x2=DispRelaxationTime_2;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        return (x2+tau*x1+tau*tau*x0)/(y2+tau*y1+tau*tau*y0);
    }
    
    /*
     * Assumed that v(k)=Phi_1(tau)*u(k)+Phi_2(tau)*u(k-1)+Phi_3(tau)*u(k-2)
     * this function returns Phi_2(tau)
     */
    public double getPhi_2(double tau){
        double x0=DispRelaxationTime_0;
        double x1=DispRelaxationTime_1;
        double x2=DispRelaxationTime_2;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        return -(2.*x2+tau*x1)/(y2+tau*y1+tau*tau*y0);
    }
    
    /*
     * Assumed that v(k)=Phi_1(tau)*u(k)+Phi_2(tau)*u(k-1)+Phi_3(tau)*u(k-2)
     * this function returns Phi_3(tau)
     */
    public double getPhi_3(double tau){
        double x0=DispRelaxationTime_0;
        double x1=DispRelaxationTime_1;
        double x2=DispRelaxationTime_2;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        return (x2)/(y2+tau*y1+tau*tau*y0);
    }
    
    /*
     * Being p(k) the 'physical' total tractions in time steo k the auxiliary traction
     * which actually computed by BEM system is t(k) and it is given as
     * t(k)=p(k)+Psi_2(tau)*p(k-1)+Psi_3(tau)*p(k-2)
     * this function returns Psi_2(tau)
     */
    public double getPsi_2(double tau){
        double x0=DispRelaxationTime_0;
        double x1=DispRelaxationTime_1;
        double x2=DispRelaxationTime_2;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        return (2.*y2+tau*y1)/(y2+tau*y1+tau*tau*y0);
    }
    
    /*
     * Being p(k) the 'physical' total tractions in time steo k the auxiliary traction
     * which actually computed by BEM system is t(k) and it is given as
     * t(k)=p(k)+Psi_2(tau)*p(k-1)+Psi_3(tau)*p(k-2)
     * this function returns Psi_3(tau)
     */
    public double getPsi_3(double tau){
        double x0=DispRelaxationTime_0;
        double x1=DispRelaxationTime_1;
        double x2=DispRelaxationTime_2;
        double y0=StressRelaxationTime_0;
        double y1=StressRelaxationTime_1;
        double y2=StressRelaxationTime_2;
        return -(y2)/(y2+tau*y1+tau*tau*y0);
    }
}
