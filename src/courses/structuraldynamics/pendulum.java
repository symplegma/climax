/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

// http://www.math.tamu.edu/~mpilant/math308/index.html
package courses.structuraldynamics;

import climax.contraption;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.GeneralPath;
import java.util.function.Function;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

/**
 *
 * @author pchr
 */
public class pendulum implements FirstOrderDifferentialEquations, contraption{
    protected static int num_pends=0;
    private final double[] masses;
    private final double[] lengths;
    private double[] th0;
    private double[] dth0;
    private double[][] th;
    private double[][] dth;
    public double g=9.81;//gravity 
    public double x0=0.0,y0=0.0;
    public double totLength;
    public double dx,dy, elen;
    protected int id,N;
    public Color SelfColor=Color.RED;
    public Color MotionColor=Color.ORANGE;
    public boolean linear=false;
    
    public pendulum(double[] masses, double[] lengths){
        num_pends++;
        id=num_pends;
        int len=Math.min(masses.length, lengths.length);
        this.masses=new double[len];
        this.lengths=new double[len];
        this.th0=new double[len];
        this.dth0=new double[len];
        totLength=0.0;
        elen=0.0;
        for(int i=0;i<len;i++){
            this.masses[i]=masses[i];
            this.lengths[i]=lengths[i];
            totLength+=lengths[i];
            this.th0[i]=0.0;
            this.dth0[i]=0.0;
            if(elen<lengths[i])elen=lengths[i];
        }
        dx = elen/20.0;
        dy=dx;
    }
    
    public pendulum(double[] masses){
        num_pends++;
        id=num_pends;
        int len=masses.length;
        this.masses=new double[len];
        this.lengths=new double[len];
        this.th0=new double[len];
        this.dth0=new double[len];
        totLength=0.0;
        elen=0.0;
        for(int i=0;i<len;i++){
            this.masses[i]=masses[i];
            this.lengths[i]=1;
            totLength+=lengths[i];
            this.th0[i]=0.0;
            this.dth0[i]=0.0;
            if(elen<lengths[i])elen=lengths[i];
        }
        dx = elen/20.0;
        dy=dx;
    }
    
    public pendulum(double leng){
        num_pends++;
        id=num_pends;
        int len=1;
        this.masses=new double[len];
        this.lengths=new double[len];
        this.th0=new double[len];
        this.dth0=new double[len];
        totLength=0.0;
        elen=0.0;
        for(int i=0;i<len;i++){
            this.masses[i]=1.0;
            this.lengths[i]=leng;
            totLength+=lengths[i];
            this.th0[i]=0.0;
            this.dth0[i]=0.0;
            if(elen<lengths[i])elen=lengths[i];
        }
        dx = elen/20.0;
        dy=dx;
    }
    
    public void setInitCond(double[] th0){
        th=new double[th0.length][1];
        this.th0=new double[th0.length]; 
        this.dth0=new double[th0.length];
        for(int i=0;i<th0.length;i++){th0[i]=th0[i];th[i][0]=th0[i];dth0[i]=0.0;}
    }
    
    public void setInitCond(double[] th0,double[] dth0){
        th=new double[th0.length][1];
        this.th0=new double[th0.length]; for(int i=0;i<th0.length;i++){this.th0[i]=th0[i];th[i][0]=th0[i];}
        this.dth0=new double[dth0.length]; for(int i=0;i<dth0.length;i++)this.dth0[i]=dth0[i];
    }
    
    @Override
    public int getID() {
        return this.id;
    }
    
    public double totLength(){
        return this.totLength;
    }
    
    public void setGravity(double grav){this.g=grav;}
    
    // methods relevant to contraption and graphical representation
    @Override
    public void SelfPortrait(Graphics2D g2, double ymin, double ymax, int ye, int ys, double xmin, double xmax, int xe, int xs) {
        Function<Double, Integer> transX = x -> {
            return ((int) ((int) ((-xmin * xe + x * (xe - xs) + xmax * xs) / (xmax - xmin))));
        };
        
        Function<Double, Integer> transY = y -> {
            //int yd;
            return ((int) ((int) ((-ymin * ye + y * (ye - ys) + ymax * ys) / (ymax - ymin))));
            //return yd;
        };
        Color origColor=g2.getColor();
        g2.setColor(this.SelfColor);
        
        // top horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0), transX.apply(x0+dx/2.0),  transY.apply(y0));
        
        // vertical
        g2.drawLine(transX.apply(x0), transY.apply(y0), transX.apply(x0),  transY.apply( (y0-this.totLength) ));
        
        double yend=y0;
        for(int i=0;i<this.lengths.length;i++){
            yend-=lengths[i];
            // mass
            GeneralPath gp = new GeneralPath();
            gp.moveTo(transX.apply(x0-dx), transY.apply(yend+dy));
            gp.lineTo(transX.apply(x0+dx), transY.apply(yend+dy));
            gp.lineTo(transX.apply(x0+dx), transY.apply(yend-dy));
            gp.lineTo(transX.apply(x0-dx), transY.apply(yend-dy));
            gp.closePath();
            g2.fill(gp);
            g2.draw(gp);
        }
        g2.setColor(origColor);
    }

    @Override
    public double maximum_X_coordinate() {
        return x0+dx;
    }

    @Override
    public double minimum_X_coordinate() {
        return x0-dx;
    }

    @Override
    public double maximum_Y_coordinate() {
        return y0;
    }

    @Override
    public double minimum_Y_coordinate() {
        return y0-dy-this.totLength;
    }

    @Override
    public void Motion(Graphics2D g2, double ymin, double ymax, int ye, int ys, double xmin, double xmax, int xe, int xs, int step, double scale) {
        Function<Double, Integer> transX = x -> {
            return ((int) ((int) ((-xmin * xe + x * (xe - xs) + xmax * xs) / (xmax - xmin))));
        };
        
        Function<Double, Integer> transY = y -> {
            //int yd;
            return ((int) ((int) ((-ymin * ye + y * (ye - ys) + ymax * ys) / (ymax - ymin))));
            //return yd;
        };
        Color origColor=g2.getColor();
        g2.setColor(this.MotionColor);
        
        double[] ux= new double[this.masses.length];
        double[] uy= new double[this.masses.length];
        
        for(int i=0;i<ux.length;i++){
            ux[i]=0.0;
            uy[i]=0.0;
            for(int j=0;j<=i;j++){
                ux[i]+=this.lengths[j]*Math.sin(scale*th[j][step]);
                uy[i]+=-this.lengths[j]*Math.cos(scale*th[j][step]);
                
            }
        }
        
        // top horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0), transX.apply(x0+dx/2.0),  transY.apply(y0));
        
        // vertical
        //g2.drawLine(transX.apply(x0), transY.apply(y0), transX.apply(x0+ux[this.masses.length-1]),  transY.apply( (y0+uy[this.masses.length-1]) ));
        
        
        for(int i=0;i<this.lengths.length;i++){
            // mass
            GeneralPath gp = new GeneralPath();
            gp.moveTo(transX.apply(x0-dx+ux[i]), transY.apply(y0+dy+uy[i]));
            gp.lineTo(transX.apply(x0+dx+ux[i]), transY.apply(y0+dy+uy[i]));
            gp.lineTo(transX.apply(x0+dx+ux[i]), transY.apply(y0-dy+uy[i]));
            gp.lineTo(transX.apply(x0-dx+ux[i]), transY.apply(y0-dy+uy[i]));
            gp.closePath();
            g2.fill(gp);
            g2.draw(gp);
            
            // vertical
            if(i==0){
                g2.drawLine(transX.apply(x0), transY.apply(y0), transX.apply(x0+ux[i]),  transY.apply( y0+uy[i] ));
            }else{
                g2.drawLine(transX.apply(x0+ux[i-1]), transY.apply(y0+uy[i-1]), transX.apply(x0+ux[i]),  transY.apply( y0+uy[i] ));
            }
        }
        g2.setColor(origColor);
    }
    
    public void setOrigin(double x0, double y0){
        this.x0=x0; this.y0=y0;
    }
    
    public void solve(double tot, double dt){
        if(this.masses.length<=2){
            ClassicalRungeKuttaIntegrator CRK = new ClassicalRungeKuttaIntegrator(dt);
            ContinuousOutputModel tracker = new ContinuousOutputModel();
            CRK.addStepHandler(tracker);
            N =(int) Math.round(tot/dt)+1;
            double[] y = new double[2*masses.length];
            for(int i=0;i<masses.length;i++){
                y[i]=this.th0[i];
                y[i+masses.length]=this.dth0[i];
            }
            CRK.integrate(this, 0.0, y, tot, y);
            double h=dt;
            th=new double[masses.length][N];
            dth=new double[masses.length][N];
            for(int it=0;it<N;it++){
                tracker.setInterpolatedTime(it*dt);
                double[] res=tracker.getInterpolatedState();
                //System.out.println("Solution step = "+it);
                for(int i=0;i<masses.length;i++){
                    th[i][it]=res[i];
                    dth[i][it]=res[i+masses.length];
                }
            }
        }else{
            throw new UnsupportedOperationException("Not supported yet for greater than two-mass pendulum."); //To change body of generated methods, choose Tools | Templates.
        }
    }

    @Override
    public int getDimension() {
        return 2*this.masses.length;
    }

    @Override
    public void computeDerivatives(double d, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        int dim=masses.length;
        if(linear){
            switch(dim){
                case 1: 
                    yDot[0] = y[1];
                    yDot[1] = -g*y[0]/this.totLength;
                    break;
                case 2:
                    yDot[0] = y[2];
                    yDot[1] = y[3];
                    double m11=(masses[0]+masses[1])*lengths[0]*lengths[0];
                    double m12=masses[1]*lengths[0]*lengths[1];
                    double m21=m12;
                    double m22=masses[1]*lengths[1]*lengths[1];
                    double den=m11*m22-m12*m21;
                    double im11=m22/den;
                    double im12=-m12/den;
                    double im21=-m21/den;
                    double im22=m11/den;
                    double k11=(masses[0]+masses[1])*lengths[0]*g;
                    double k22=masses[1]*lengths[1]*g;
                    double imk11=im11*k11;
                    double imk12=im12*k22;
                    double imk21=im21*k11;
                    double imk22=im22*k22;
                    yDot[2] = -imk11*y[0]-imk12*y[1];
                    yDot[3] = -imk21*y[0]-imk22*y[1];
            }
        }else{
            switch(dim){
                case 1: 
                    yDot[0] = y[1];
                    yDot[1] = -g*Math.sin(y[0])/this.totLength;
                    break;
                case 2:
                    yDot[0] = y[2];
                    yDot[1] = y[3];
                    double mu=1.0+(masses[0]+masses[1]);
                    double l1=lengths[0];
                    double l2=lengths[1];
                    double sin1=Math.sin(y[0]);
                    double sin2=Math.sin(y[1]);
                    double cos1m2=Math.cos(y[0]-y[1]);
                    double sin1m2=Math.sin(y[0]-y[1]);
                    // tesis of advisor kapitaniak
                    //yDot[2] = ( g*(sin2*cos1m2-mu*sin1)-(l2*y[3]*y[3]+l1*y[2]*y[2]*cos1m2)*sin1m2 )/( l1*(mu-cos1m2*cos1m2) );
                    //yDot[3] = ( g*mu*(sin1*cos1m2-sin2)+(mu*l1*y[2]*y[2]+l2*y[3]*y[3]*cos1m2)*sin1m2 )/( l2*(mu-cos1m2*cos1m2) );
                    double sin1m2t2=Math.sin(y[0]-2*y[1]);
                    double cos2t1m2t2=Math.cos(2*y[0]-2*y[1]);
                    double m1=masses[0];
                    double m2=masses[1];
                    double cos1=Math.cos(y[0]);
                    // myphysicslab.com
                    yDot[2] = ( -g*(2*m1+m2)*sin1-m2*g*sin1m2t2-2*sin1m2*m2*(y[3]*y[3]*l2+y[2]*y[2]*l1*cos1m2) )/( l1*(2*m1+m2-m2*cos2t1m2t2) );
                    yDot[3] = ( 2*sin1m2*(y[2]*y[2]*l1*(m1+m2)+g*(m1+m2)*cos1+y[3]*y[3]*l2*cos1m2) )/( l2*(2*m1+m2-m2*cos2t1m2t2) );
            }
        }
    }
    
    public double[] Theta(int dof){
        double[] val = new double[N];
        for(int i=0;i<N;i++)val[i]=th[dof-1][i];
        return val;
    }
    
    public double[] Theta(){
        return Theta(1);
    }
    
    public double[] DTheta(int dof){
        double[] val = new double[N];
        for(int i=0;i<N;i++)val[i]=dth[dof-1][i];
        return val;
    }
    
    public double[] DTheta(){
        return DTheta(1);
    }
    
    public double[] Xdef(int dof){
        double[] val = new double[N];
        for(int i=0;i<N;i++){
            val[i]=x0;
            for(int j=0;j<dof;j++){
                val[i]+=this.lengths[j]*Math.sin(th[j][i]);
            }
        }
        return val;
    }
    
    public double[] Ydef(int dof){
        double[] val = new double[N];
        for(int i=0;i<N;i++){
            val[i]=y0;
            for(int j=0;j<dof;j++){
                val[i]+=-this.lengths[j]*Math.cos(th[j][i]);
                val[i]+=this.lengths[j];
            }
        }
        return val;
    }
}
