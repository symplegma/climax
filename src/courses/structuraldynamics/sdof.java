/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package courses.structuraldynamics;


import climax.contraption;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.GeneralPath;
import java.util.function.Function;
import mathman.DoubleFunction;
import mathman.ZeroDoubleFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 *
 * @author pchr
 */
public class sdof  implements FirstOrderDifferentialEquations, contraption{
    protected static int num_sdof=0;
    double k,m,c;
    DoubleFunction df;
    double u0=0.0,v0=0.0;
    double[] u,v,a;
    double h,Atkm1=0.0,Btkm1=0.0;
    double xsi,omega,omegaD;
    double x0=0.0,y0=0.0, elen=10.0;
    double dx=2.0,dy=2.0;
    protected int id;
    protected Color SelfColor=Color.RED;
    protected Color MotionColor=Color.ORANGE;
    
    public sdof(double k, double m, double c){
        num_sdof++;
        id=num_sdof;
        this.k=k;
        this.m=m;
        this.c=c;
        this.omega=Math.sqrt(k/m);
        this.xsi=c/(2.0*m*omega);
        omegaD=omega*Math.sqrt(Math.abs(1-xsi*xsi));
        df=new ZeroDoubleFunction();
    }
    
    public sdof(double k, double m, double c, double u0, double v0){
        num_sdof++;
        id=num_sdof;
        this.k=k;
        this.m=m;
        this.c=c;
        this.u0=u0;
        this.v0=v0;
        this.omega=Math.sqrt(k/m);
        this.xsi=c/(2.0*m*omega);
        omegaD=omega*Math.sqrt(Math.abs(1-xsi*xsi));
        df=new ZeroDoubleFunction();
    }
    
    public sdof(double k, double m){
        num_sdof++;
        id=num_sdof;
        this.k=k;
        this.m=m;
        this.c=0.0;
        this.u0=0.0;
        this.v0=0.0;
        this.omega=Math.sqrt(k/m);
        this.xsi=c/(2.0*m*omega);
        omegaD=omega*Math.sqrt(Math.abs(1-xsi*xsi));
        df=new ZeroDoubleFunction();
    }
    
    public void setCritDampRatio(double xsi){
        this.xsi=xsi;
        c=2.0*m*omega*xsi;
        omegaD=omega*Math.sqrt(Math.abs(1-xsi*xsi));
    }
    
    public void setRHS(DoubleFunction df){this.df=df;}
    
    public void setInitCond(double u0){
        u=new double[1]; u[0]=u0;
        this.u0=u0;
        this.v0=0.0;
    }
    
    public void setInitCond(double u0, double v0){
        u=new double[1]; u[0]=u0;
        this.u0=u0;
        this.v0=v0;
    }

    @Override
    public int getDimension() {
        return 2; //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        yDot[0] = y[1];
        yDot[1] = -y[1]*c/m-y[0]*k/m+df.run(t)/m;
    }
    
    public void solve(){
        double T=2.0*Math.PI/omega;
        solve(T, T/100.0);
    }
    
    public void solve(int N){
        double T=2.0*Math.PI/omega;
        solve(N*T, T/100.0);
    }
    
    public void solve(int N, double dt){
        double T=2.0*Math.PI/omega;
        solve(N*T, dt);
    }
    
    public void solve(int N, int M){
        double T=2.0*Math.PI/omega;
        solve(N*T, T/M);
    }
    
    public void solve(double tot, double dt){
        ClassicalRungeKuttaIntegrator CRK = new ClassicalRungeKuttaIntegrator(dt);
        ContinuousOutputModel tracker = new ContinuousOutputModel();
        CRK.addStepHandler(tracker);
        int N =(int) Math.round(tot/dt)+1;
        double[] y = new double[2];
        y[0]=u0; y[1]=v0;
        CRK.integrate(this, 0.0, y, tot, y);
        this.h=dt;
        u=new double[N];
        v=new double[N];
        a=new double[N];
        for(int it=0;it<N;it++){
            tracker.setInterpolatedTime(it*dt);
            double[] res=tracker.getInterpolatedState();            
            u[it]=res[0];
            v[it]=res[1];
            res=tracker.getInterpolatedDerivatives();
            a[it]=res[1];
        }
    }
    
    public double[] Disp(){
        return this.u;
    }
    
    public double[] Velc(){
        return this.v;
    }
    
    public double[] Accl(){
        return this.a;
    }
    
    public double dt(){return this.h;}
    
    public double getCritDampRatio(){
        return xsi;
    }
    
    public double getNaturalFrequency(){return omega;}
    
    public double getNaturalPeriod(){return 2.0*Math.PI/omega;}
    
    public double maxDisp(){
        double val = Double.NEGATIVE_INFINITY;
        for(int i=0;i<this.u.length;i++){
            if(u[i]>val)val=u[i];
        }
        return val;
    }
    
    public double minDisp(){
        double val = Double.POSITIVE_INFINITY;
        for(int i=0;i<this.u.length;i++){
            if(u[i]<val)val=u[i];
        }
        return val;
    }
    
    public double maxVelc(){
        double val = Double.NEGATIVE_INFINITY;
        for(int i=0;i<this.v.length;i++){
            if(v[i]>val)val=v[i];
        }
        return val;
    }
    
    public double minVelc(){
        double val = Double.POSITIVE_INFINITY;
        for(int i=0;i<this.v.length;i++){
            if(v[i]<val)val=v[i];
        }
        return val;
    }
    
    public void duhamel(){
        double T=2.0*Math.PI/omega;
        duhamel(T, T/100.0);
    }
    
    public void duhamel(int N){
        double T=2.0*Math.PI/omega;
        duhamel(N*T, T/100.0);
    }
    
    public void duhamel(int N, double dt){
        double T=2.0*Math.PI/omega;
        duhamel(N*T, dt);
    }
    
    public void duhamel(int N, int M){
        double T=2.0*Math.PI/omega;
        duhamel(N*T, T/M);
    }
    
    public void duhamel(double tot, double dt){
        int N =(int) Math.round(tot/dt)+1;
        Atkm1=0.0; Btkm1=0.0;
        this.h=dt;
        u=new double[N]; u[0]=u0;
        v=new double[N]; v[0]=v0;
        double val=0; if(df!=null){val=df.run(0.0);}
        a=new double[N]; a[0]=(val-c*v0-k*u0)/m;
        for(int i=1;i<N;i++){
            val=0; if(df!=null){val=df.run(h*i);}
            double val0,dval0;
            if(xsi<1){
                double Ai=A(i);
                double Bi=B(i);
                val0=((v0+u0*xsi*omega)*Math.sin(omegaD*i*h)/omegaD+u0*Math.cos(omegaD*i*h))*Math.exp(-xsi*omega*i*h);
                dval0=-xsi*omega*((v0+u0*xsi*omega)*Math.sin(omegaD*i*h)/omegaD+u0*Math.cos(omegaD*i*h))*Math.exp(-xsi*omega*i*h)
                        +((v0+u0*xsi*omega)*Math.cos(omegaD*i*h)-u0*omegaD*Math.sin(omegaD*i*h))*Math.exp(-xsi*omega*i*h);
                
                
                u[i]=Math.exp(-xsi*omega*h*i)*(Ai*Math.sin(omegaD*h*i)-Bi*Math.cos(omegaD*h*i))/(m*omegaD);
                v[i]=-xsi*omega*u[i]+Math.exp(-xsi*omega*h*i)*(Ai*Math.cos(omegaD*h*i)+Bi*Math.sin(omegaD*h*i))/m+dval0;
                u[i]+=val0;
                a[i]=-2*xsi*omega*v[i]-omega*omega*u[i]+val/m;
            }
            if(xsi>=1){
                if(omegaD<=Double.MIN_VALUE){
                    double Aci=Ac(i);
                    double Bci=Bc(i);
                    val0=(u0*(1.0+omega*i*h)+v0*i*h)*Math.exp(-omega*i*h);
                    dval0=-omega*(u0*(1.0+omega*i*h)+v0*i*h)*Math.exp(-omega*i*h)+
                            (u0*omega+v0)*Math.exp(-omega*i*h);
                    
                    u[i]=Math.exp(-omega*i*h)*(Aci*i*h-Bci)/m;
                    v[i]=-omega*u[i]+Math.exp(-omega*i*h)*Aci/m+dval0;
                    u[i]+=val0;
                    a[i]=-2*xsi*omega*v[i]-omega*omega*u[i]+val/m;
                }else{
                    double Ahi=Ah(i);
                    double Bhi=Bh(i);
                    val0=((v0+u0*xsi*omega)*Math.sinh(omegaD*i*h)/omegaD+u0*Math.cosh(omegaD*i*h))*Math.exp(-xsi*omega*i*h);
                    dval0=-xsi*omega*((v0+u0*xsi*omega)*Math.sinh(omegaD*i*h)/omegaD+u0*Math.cosh(omegaD*i*h))*Math.exp(-xsi*omega*i*h)
                    +((v0+u0*xsi*omega)*Math.cosh(omegaD*i*h)+u0*omegaD*Math.sinh(omegaD*i*h))*Math.exp(-xsi*omega*i*h);
                    
                    u[i]=Math.exp(-xsi*omega*h*i)*(Ahi*Math.sinh(omegaD*h*i)-Bhi*Math.cosh(omegaD*h*i))/(m*omegaD);
                    v[i]=-xsi*omega*u[i]+Math.exp(-xsi*omega*h*i)*(Ahi*Math.cosh(omegaD*h*i)+Bhi*Math.sinh(omegaD*h*i))/m+dval0;
                    u[i]+=val0;
                    a[i]=-2*xsi*omega*v[i]-omega*omega*u[i]+val/m;
                }
            }
        }
    }
    
    private double A(int k){
        double p2km2=p(2*(k-1)*h/2.0);
        double p2km1=p((2*k-1)*h/2.0);
        double p2k=p((2*k)*h/2.0);
        Atkm1+=h*(p2km2+4.0*p2km1+p2k)/6.0;
        return Atkm1;
    }
    
    private double B(int k){
        double g2km2=g(2*(k-1)*h/2.0);
        double g2km1=g((2*k-1)*h/2.0);
        double g2k=g((2*k)*h/2.0);
        Btkm1+=h*(g2km2+4.0*g2km1+g2k)/6.0;
        return Btkm1;
    }
    
    private double p(double tk){
        double val=0; if(df!=null){val=df.run(tk);}
        return val*Math.exp(xsi*omega*tk)*Math.cos(omegaD*tk);
    }
    
    private double g(double tk){
        double val=0; if(df!=null){val=df.run(tk);}
        return val*Math.exp(xsi*omega*tk)*Math.sin(omegaD*tk);
    }
    
    private double Ac(int k){
        double p2km2=pc(2*(k-1)*h/2.0);
        double p2km1=pc((2*k-1)*h/2.0);
        double p2k=pc((2*k)*h/2.0);
        Atkm1+=h*(p2km2+4.0*p2km1+p2k)/6.0;
        return Atkm1;
    }
    
    private double Bc(int k){
        double g2km2=gc(2*(k-1)*h/2.0);
        double g2km1=gc((2*k-1)*h/2.0);
        double g2k=gc((2*k)*h/2.0);
        Btkm1+=h*(g2km2+4.0*g2km1+g2k)/6.0;
        return Btkm1;
    }
    
    private double pc(double tk){
        double val=0; if(df!=null){val=df.run(tk);}
        return val*Math.exp(omega*tk);
    }
    
    private double gc(double tk){
        double val=0; if(df!=null){val=df.run(tk);}
        return val*Math.exp(omega*tk)*tk;
    }
    
     private double Ah(int k){
        double p2km2=ph(2*(k-1)*h/2.0);
        double p2km1=ph((2*k-1)*h/2.0);
        double p2k=ph((2*k)*h/2.0);
        Atkm1+=h*(p2km2+4.0*p2km1+p2k)/6.0;
        return Atkm1;
    }
    
    private double Bh(int k){
        double g2km2=gh(2*(k-1)*h/2.0);
        double g2km1=gh((2*k-1)*h/2.0);
        double g2k=gh((2*k)*h/2.0);
        Btkm1+=h*(g2km2+4.0*g2km1+g2k)/6.0;
        return Btkm1;
    }
    
    private double ph(double tk){
        double val=0; if(df!=null){val=df.run(tk);}
        return val*Math.exp(xsi*omega*tk)*Math.cosh(omegaD*tk);
    }
    
    private double gh(double tk){
        double val=0; if(df!=null){val=df.run(tk);}
        return val*Math.exp(xsi*omega*tk)*Math.sinh(omegaD*tk);
    }
    
    public double getStiffness(){
        return this.k;
    }
    
    public double getMass(){
        return this.m;
    }
    
    public double getDamping(){
        return this.c;
    }
    
    // methods relevant to contraption and graphical representation
    public void setOrigin(double x0, double y0){
        this.x0=x0; this.y0=y0;
    }
    
    public void setLength(double l){this.elen=l;}
    
    public void setWidth(double w){this.dx=w/2.0;}
    
    public void setHeight(double h){this.dy=h/2.0;}

    public int getID() {
        return this.id;
    }

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
        
        // top vertical dashpot and spring
        g2.drawLine(transX.apply(x0-dx/4.0), transY.apply(y0), transX.apply(x0-dx/4.0),  transY.apply( (y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/4.0), transY.apply(y0), transX.apply(x0+dx/4.0),  transY.apply( (y0-elen)/3.0 ));
        
        // dashpot
        g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( (y0-elen)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( (y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( (y0-elen)/3.0 ), transX.apply(x0-dx/2.0+dx/8.0),  transY.apply( (y0-elen)/3.0+(y0-elen)/6.0 ));
        g2.drawLine(transX.apply(x0-dx/8.0), transY.apply( (y0-elen)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( (y0-elen)/3.0+(y0-elen)/6.0 ));
        
        g2.drawLine(transX.apply(x0-dx/4.0), transY.apply((y0-elen)/3.0+(y0-elen)/10.0 ), transX.apply(x0-dx/4.0),  transY.apply(2.0*(y0-elen)/3.0));
        g2.drawLine(transX.apply(x0-dx/4.0-dx/10.0), transY.apply((y0-elen)/3.0+(y0-elen)/10.0 ), transX.apply(x0-dx/4.0+dx/10.0),  transY.apply((y0-elen)/3.0+(y0-elen)/10.0 ));
        
        //spring
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply((y0-elen)/3.0),                       transX.apply(x0+dx/8.0),         transY.apply((y0-elen)/3.0));
        
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply((y0-elen)/3.0),                       transX.apply(x0+dx/8.0),         transY.apply( (y0-elen)/3.0 + 0.1*(y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( (y0-elen)/3.0 + 0.1*(y0-elen)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( (y0-elen)/3.0 + 0.2*(y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( (y0-elen)/3.0 + 0.2*(y0-elen)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( (y0-elen)/3.0 + 0.3*(y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( (y0-elen)/3.0 + 0.3*(y0-elen)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( (y0-elen)/3.0 + 0.4*(y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( (y0-elen)/3.0 + 0.4*(y0-elen)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( (y0-elen)/3.0 + 0.5*(y0-elen)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( (y0-elen)/3.0 + 0.5*(y0-elen)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( (y0-elen)/3.0 + 0.6*(y0-elen)/3.0 ));
        
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( (y0-elen)/3.0 + 0.6*(y0-elen)/3.0 ),                       transX.apply(x0+dx/8.0),         transY.apply( (y0-elen)/3.0 + 0.6*(y0-elen)/3.0 ));
        
        // bottom vertical dashpot and spring 
        g2.drawLine(transX.apply(x0-dx/4.0), transY.apply( (y0-elen)/3.0 + 0.6*(y0-elen)/3.0 ), transX.apply(x0-dx/4.0),  transY.apply( (y0-elen) ));
        g2.drawLine(transX.apply(x0+dx/4.0), transY.apply( (y0-elen)/3.0 + 0.6*(y0-elen)/3.0 ), transX.apply(x0+dx/4.0),  transY.apply( (y0-elen) ));
        
        // mass
        GeneralPath gp = new GeneralPath();
        gp.moveTo(transX.apply(x0-dx), transY.apply(y0-elen));
        gp.lineTo(transX.apply(x0+dx), transY.apply(y0-elen));
        gp.lineTo(transX.apply(x0+dx), transY.apply(y0-2*dy-elen));
        gp.lineTo(transX.apply(x0-dx), transY.apply(y0-2*dy-elen));
        gp.closePath();
        g2.fill(gp);
        g2.draw(gp);
        
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
        return y0-2*dy-elen;
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
        double disp=-this.u[step]*scale;
        // top horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0), transX.apply(x0+dx/2.0),  transY.apply(y0));
        
        // top vertical dashpot and spring
        g2.drawLine(transX.apply(x0-dx/4.0), transY.apply(y0), transX.apply(x0-dx/4.0),  transY.apply( (y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/4.0), transY.apply(y0), transX.apply(x0+dx/4.0),  transY.apply( (y0-elen+disp)/3.0 ));
        
        // dashpot
        g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( (y0-elen+disp)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( (y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( (y0-elen+disp)/3.0 ), transX.apply(x0-dx/2.0+dx/8.0),  transY.apply( (y0-elen+disp)/3.0+(y0-elen+disp)/6.0 ));
        g2.drawLine(transX.apply(x0-dx/8.0), transY.apply( (y0-elen+disp)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( (y0-elen+disp)/3.0+(y0-elen+disp)/6.0 ));
        
        g2.drawLine(transX.apply(x0-dx/4.0), transY.apply((y0-elen+disp)/3.0+(y0-elen+disp)/10.0 ), transX.apply(x0-dx/4.0),  transY.apply(2.0*(y0-elen+disp)/3.0));
        g2.drawLine(transX.apply(x0-dx/4.0-dx/10.0), transY.apply((y0-elen+disp)/3.0+(y0-elen+disp)/10.0 ), transX.apply(x0-dx/4.0+dx/10.0),  transY.apply((y0-elen+disp)/3.0+(y0-elen+disp)/10.0 ));
        
        //spring
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply((y0-elen+disp)/3.0),                       transX.apply(x0+dx/8.0),         transY.apply((y0-elen+disp)/3.0));
        
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply((y0-elen+disp)/3.0),                       transX.apply(x0+dx/8.0),         transY.apply( (y0-elen+disp)/3.0 + 0.1*(y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( (y0-elen+disp)/3.0 + 0.1*(y0-elen+disp)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( (y0-elen+disp)/3.0 + 0.2*(y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( (y0-elen+disp)/3.0 + 0.2*(y0-elen+disp)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( (y0-elen+disp)/3.0 + 0.3*(y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( (y0-elen+disp)/3.0 + 0.3*(y0-elen+disp)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( (y0-elen+disp)/3.0 + 0.4*(y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( (y0-elen+disp)/3.0 + 0.4*(y0-elen+disp)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( (y0-elen+disp)/3.0 + 0.5*(y0-elen+disp)/3.0 ));
        g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( (y0-elen+disp)/3.0 + 0.5*(y0-elen+disp)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( (y0-elen+disp)/3.0 + 0.6*(y0-elen+disp)/3.0 ));
        
        g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( (y0-elen+disp)/3.0 + 0.6*(y0-elen+disp)/3.0 ),                       transX.apply(x0+dx/8.0),         transY.apply( (y0-elen+disp)/3.0 + 0.6*(y0-elen+disp)/3.0 ));
        
        // bottom vertical dashpot and spring 
        g2.drawLine(transX.apply(x0-dx/4.0), transY.apply( (y0-elen+disp)/3.0 + 0.6*(y0-elen+disp)/3.0 ), transX.apply(x0-dx/4.0),  transY.apply( (y0-elen) +disp));
        g2.drawLine(transX.apply(x0+dx/4.0), transY.apply( (y0-elen+disp)/3.0 + 0.6*(y0-elen+disp)/3.0 ), transX.apply(x0+dx/4.0),  transY.apply( (y0-elen) +disp));
        
        GeneralPath gp = new GeneralPath();
        gp.moveTo(transX.apply(x0-dx), transY.apply(y0-elen+disp));
        gp.lineTo(transX.apply(x0+dx), transY.apply(y0-elen+disp));
        gp.lineTo(transX.apply(x0+dx), transY.apply(y0-2*dy-elen+disp));
        gp.lineTo(transX.apply(x0-dx), transY.apply(y0-2*dy-elen+disp));
        gp.closePath();
        g2.fill(gp);
        g2.draw(gp);
        
        g2.setColor(origColor);
    }
    
    Complex TransferFunction(double omext){
        Complex H = new Complex(k-m*omext*omext, c*omext);
        //Complex H = new Complex(m*omega*omega-m*omext*omext, 2.0*omega*xsi*omext*m);
        Complex UnitReal = new Complex(1.0, 0.0);
        return UnitReal.divide(H);
    }
}
