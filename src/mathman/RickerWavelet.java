/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mathman;

/**
 *
 * @author pchr
 */
public class RickerWavelet implements DoubleFunction{
    double f0=1.0;
    double t0=0.0;
    
    public RickerWavelet(){}
    
    public RickerWavelet(double f0){
        this.f0=f0;
    }
    
    public RickerWavelet(double f0,double t0){
        this.f0=f0;
        this.f0=t0;
    }
    
    public void setFreq(double val){this.f0=val;}
    
    public void setPoint(double val){this.t0=val;}

    @Override
    public double run(double x) {
        double term=Math.PI*f0*(x-t0);
        return (1.0-2*term*term)*Math.exp(-term*term);
    }

    @Override
    public double run(double x, double y) {
       double term=Math.PI*f0*(x-t0)*(y-t0);
        return (1.0-2*term*term)*Math.exp(-term*term);
    }
    
}
