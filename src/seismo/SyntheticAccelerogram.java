/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seismo;

/**
 *
 * @author pchr
 */
public class SyntheticAccelerogram {
    private double So=150; // m2/s3
    private double om2=5.; // K-T (Kanai-Tajimi)
    private double xsi2=0.7;
    private double om1=0.8; // H-P (high pass)
    private double xsi1=1.0;
    private double om3=10.; // L-P (low pass) actually not used here
    private double xsi3=1.4;
    private int N=500;
    private double T=2.5;
    private double t1=0.;
    private double t2=2.5;
    private double Om_a=0.1;
    private double Om_b=100.;
    private boolean K_T=true;
    private boolean H_P=true;
    private boolean L_P=true;
    private double[] acc,vel,dis;
    
    // constructors
    public SyntheticAccelerogram(){}
    
    public SyntheticAccelerogram(double So, double om1, double xsi1, double om2, double xsi2){
        this.So=So;
        this.om1=om1;
        this.om2=om2;
        this.xsi1=xsi1;
        this.xsi2=xsi2;
    }
    
    public void setN(int N){this.N=N;}
    
    public void setT(double T){
        this.T=T;
        t1=0.;
        t2=T;
    }
    
    public void setSo(double So){this.So=So;}
    
    public void set_t1(double t1){this.t1=t1;}
    
    public void set_t2(double t2){this.t2=t2;}
    
    public void setOm_a(double Om_a){this.Om_a=Om_a;}
    
    public void setOm_b(double Om_b){this.Om_b=Om_b;}
    
    public void setK_Tfilter(boolean bool){this.K_T=bool;}
    
    public void setH_Pfilter(boolean bool){this.H_P=bool;}
    
    public void setL_Pfilter(boolean bool){this.L_P=bool;}
    
    private double PSD(double Om){
        double val=So;
        if(K_T)val*=((1+4*xsi2*(Om/om2)*(Om/om2))/((1-(Om/om2)*(Om/om2))*(1-(Om/om2)*(Om/om2)) + 4*xsi2*xsi2*(Om/om2)*(Om/om2)));
        if(H_P)val*=((Om/om1)*(Om/om1)*(Om/om1)*(Om/om1)/((1-(Om/om1)*(Om/om1))*(1-(Om/om1)*(Om/om1)) + 4*xsi1*xsi1*(Om/om1)*(Om/om1)));
        if(L_P)val*=(1./((1-(Om/om3)*(Om/om3))*(1-(Om/om3)*(Om/om3)) + 4*xsi3*xsi3*(Om/om3)*(Om/om3)));
        return val;
    }
    
    public void Calc(int steps){
        double dt=T/((double)steps);
        acc =new double[steps+1];
        vel =new double[steps+1];
        dis =new double[steps+1];
        for(int i=0;i<acc.length;i++){
            acc[i]=0.0;
            vel[i]=0.0;
            dis[i]=0.0;
        }
        for(int n=1;n<=N;n++){
            double dom=(Om_b-Om_a)/N;
            double Omn=(n-0.5)*dom;
            double An=Math.sqrt(2*PSD(Omn)*dom);
            java.util.Random r = new java.util.Random();
            double phin=r.nextGaussian()*2.0*Math.PI;
            for(int i=0;i<acc.length;i++){
                acc[i]+=An*Math.sin(Omn*i*dt+phin)*psi(i*dt);
            }
        }
        
        for(int i=1;i<acc.length;i++){
            vel[i]=vel[i-1]+(1-0.5)*dt*acc[i-1]+0.5*dt*acc[i];
            dis[i]=dis[i-1]+dt*vel[i-1]+(0.5-1./6.)*dt*dt*acc[i-1]+dt*dt*acc[i]/6.;
        }
    }
    
    public double[] getSynthAcc(){
        return acc;
    }
    
    public double[] getSynthDis(){
        return dis;
    }
    
    public double[] getSynthVel(){
        return vel;
    }
    
    private double psi(double t){
        double val=0.;
        double b=2./(t2-t1);
        if(t2<t1){
            System.out.println("t1>t2, changed to t1=0 and t2=T.");
            t1=0.;t2=T;
        }
        if(t<t1){
            val=(t/t1)*(t/t1);
        }else if(t>=t1 && t<=t2){
            val=1.;
        }else{
            val=Math.exp(-b*(t-t2));
        }
        return val;
    }
}
