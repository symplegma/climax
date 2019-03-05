/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

//import mathman.Complex;
import org.apache.commons.math3.complex.Complex;
/**
 *
 * @author pchr
 */
public class OrthotropicMaterial2D extends Material{
    private double ElasticModulus_1;
    private double ElasticModulus_2;
    private double ElasticModulus_3;
    private double PoissonRatio_12;
    private double PoissonRatio_21;
    private double PoissonRatio_13;
    private double PoissonRatio_31;
    private double PoissonRatio_23;
    private double PoissonRatio_32;
    private double ShearModulus_12;
    
    private boolean PlaneStress=false;
    
    // constructor
    public OrthotropicMaterial2D(){}
    
    public OrthotropicMaterial2D(int id, double E1, double E2, double E3, double G12, double v12, double v13, double v23){
        this.id=id;
        this.ElasticModulus_1=E1;
        this.ElasticModulus_2=E2;
        this.ElasticModulus_3=E3;
        this.ShearModulus_12=G12;
        
        
        this.PoissonRatio_12=v12;
        this.PoissonRatio_13=v13;
        this.PoissonRatio_23=v23;
        
        if(PlaneStress){
            this.PoissonRatio_21=this.PoissonRatio_12*E2/E1;
            this.PoissonRatio_31=this.PoissonRatio_13*E3/E1;
            this.PoissonRatio_32=this.PoissonRatio_23*E3/E2;
            if(Math.abs((-PoissonRatio_12/ElasticModulus_1)-(-PoissonRatio_21/ElasticModulus_2))>=jmat.MachinePrecision.getMachinePrecision()){
                System.out.println("Equation I.2.5 is not valid");
            }
        }else{
            this.PoissonRatio_21=this.PoissonRatio_12*E2/E1;
            this.PoissonRatio_31=this.PoissonRatio_13*E3/E1;
            this.PoissonRatio_32=this.PoissonRatio_23*E3/E2;
//            v=v/(1.+v);
//            this.PoissonRatio_21=(E1*E2-E2*E3*v13*v13-E2*E3*v12*v13*v23-E1*E3*v23*v23)/(E1*(E2*v12+E3*v13*v23));
//            this.PoissonRatio_31=(-E1*E1*E2+E1*E2*E2*v12*v12+2.*E1*E2*E3*v13*v13-E2*E3*E3*v13*v13*v13*v13+3.*E1*E2*E3*v12*v13*v23-E2*E3*E3*v12*v13*v13*v13*v23+E1*E1*E3*v23*v23)/(E1*(E2*v12*v13+E1*v23)*(E2*v12+E3*v13*v23));
//            this.PoissonRatio_32=(E1-E2*v12*v12-E3*v13*v13-E3*v12*v13*v23)/(E2*v12*v13+E1*v23);
            if(Math.abs((-(PoissonRatio_12+PoissonRatio_13*PoissonRatio_32)/ElasticModulus_1)-(-(PoissonRatio_21+PoissonRatio_23*PoissonRatio_31)/ElasticModulus_2))>=jmat.MachinePrecision.getMachinePrecision()){
                System.out.println("Equation I.2.6 is not valid");
            }
        }
    }
    
    public OrthotropicMaterial2D(int id, double E, double v){
        this.id=id;
        this.ElasticModulus_1=E;
        this.ElasticModulus_2=E;
        this.ElasticModulus_3=E;
        this.ShearModulus_12=E/(2.+2.*v);
        
        
        this.PoissonRatio_12=v;
        this.PoissonRatio_13=v;
        this.PoissonRatio_23=v;
        
        this.PoissonRatio_21=v;
        this.PoissonRatio_31=v;
        this.PoissonRatio_32=v;
    }
    
    public void setPlaneStress(boolean val){this.PlaneStress=val;}
    
    public boolean getPlaneStress(){return this.PlaneStress;}
    
    public double[] getBetas(){
        double[] betas = new double[4];
        if(this.PlaneStress){
            betas[0]=1./ElasticModulus_1; // b11
            betas[1]=1./ElasticModulus_2; // b22
            betas[2]=-PoissonRatio_12/ElasticModulus_1; // b12
            betas[3]=1./ShearModulus_12; //b66
        }else{
            betas[0]=(1.-PoissonRatio_13*PoissonRatio_31)/ElasticModulus_1; // b11
            betas[1]=(1.-PoissonRatio_23*PoissonRatio_32)/ElasticModulus_2; // b22
            betas[2]=-(PoissonRatio_12+PoissonRatio_13*PoissonRatio_32)/ElasticModulus_1; // b12
            betas[3]=1./ShearModulus_12; //b66
        }
        
        return betas;
    }
    
    
    public Complex getBetaMinus(){
        Complex beta;
        double gamma;
        gamma=getBetas()[3]/2.;
        double val=Math.sqrt(getBetas()[0]*getBetas()[1])-(getBetas()[2]+gamma);
        if(val>=0.){
        beta = new Complex(Math.sqrt(val),0.0);
        }else{
            beta = new Complex(0.0, Math.sqrt(-val));
        }
        return beta;
    }
    
    public Complex getBetaPlus(){
        Complex beta;
        double gamma=getBetas()[3]/2.;
        double val=Math.sqrt(getBetas()[0]*getBetas()[1])+(getBetas()[2]+gamma);
        if(val>=0.){
            beta = new Complex(Math.sqrt(val),0.0);
        }else{
            beta = new Complex(0.0, Math.sqrt(val));
        }
        return beta;
    }
    
    public Complex getRoot_1(){
        Complex beta;
        beta=(this.getBetaMinus().add(this.getBetaPlus().multiply(new Complex(0.,1.)))).divide(new Complex(Math.sqrt(2.*this.getBetas()[0]),0.));
        return beta;
    }
    
    public Complex getRoot_2(){
        Complex beta;
        beta=((this.getBetaMinus().multiply(new Complex(-1.0,0.))).add(this.getBetaPlus().multiply(new Complex(0.,1.)))).divide(new Complex(Math.sqrt(2.*this.getBetas()[0]),0.));
        return beta;
    }
    
    public enum OrthotropicMaterialCase {
        REAL,IMAGINARY,ZERO
    }
    
    public OrthotropicMaterialCase getMaterialCase(){
        OrthotropicMaterialCase theCase = null;
        double tol=jmat.MachinePrecision.getMachinePrecision();
        if(Math.signum(getRoot_1().getReal()*getRoot_2().getReal())==-1){
            if(Math.abs(this.getBetaMinus().getImaginary())>=tol){
                System.err.println("orthotropic material case REAL with imaginary part:");
                System.err.println(this.getBetaMinus().getImaginary());
            }
            theCase=OrthotropicMaterialCase.REAL;
        }else{
            if( Math.abs(getRoot_1().getImaginary()-getRoot_2().getImaginary())<=tol){
                theCase=OrthotropicMaterialCase.ZERO;
            }else{
                if(Math.abs(this.getRoot_1().getReal())>=tol){
                    System.err.println("orthotropic material case IMAGINARY with real part:");
                    System.err.println(this.getBetaMinus().getReal());
                }
                if(Math.abs(this.getRoot_1().getImaginary())<=0.0){
                    System.err.println("orthotropic material case IMAGINARY with non-positive imaginary part:");
                    System.err.println(this.getBetaMinus().getImaginary());
                }
                if(Math.abs(this.getRoot_2().getReal())>=tol){
                    System.err.println("orthotropic material case IMAGINARY with real part:");
                    System.err.println(this.getBetaMinus().getReal());
                }
                if(Math.abs(this.getRoot_2().getImaginary())<=0.0){
                    System.err.println("orthotropic material case IMAGINARY with non-positive imaginary part:");
                    System.err.println(this.getBetaMinus().getImaginary());
                }
                theCase=OrthotropicMaterialCase.IMAGINARY;
            }
        }
        return theCase;
    }
    
    public void print(){
        double[] vec=this.getBetas();
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println(this.getMaterialCase());
        System.out.println("Orthotropic material data");
        System.out.println("E1 = "+this.ElasticModulus_1);
        System.out.println("E2 = "+this.ElasticModulus_2);
        System.out.println("E3 = "+this.ElasticModulus_3);
        System.out.println("G12 = "+this.ShearModulus_12);
        System.out.println("v12 = "+this.PoissonRatio_12);
        System.out.println("v13 = "+this.PoissonRatio_13);
        System.out.println("v23 = "+this.PoissonRatio_23);
        System.out.println("v21 = "+this.PoissonRatio_21);
        System.out.println("v31 = "+this.PoissonRatio_31);
        System.out.println("v32 = "+this.PoissonRatio_32);
        System.out.println("-------------------------");
        System.out.println("β11 = "+vec[0]);
        System.out.println("β22 = "+vec[1]);
        System.out.println("β12 = "+vec[2]);
        System.out.println("β66 = "+vec[3]);
        System.out.println("-------------------------");
        System.out.println("β_+ = "+this.getBetaPlus().toString());
        if(Math.abs(this.getBetaPlus().getImaginary())>=jmat.MachinePrecision.getMachinePrecision())System.out.println("β_+ IS EXPECTED TO BE REAL POSITIVE");
        if(this.getBetaPlus().getReal()<0.0)System.out.println("β_+ IS EXPECTED TO BE REAL POSITIVE");
        System.out.println("β_- = "+this.getBetaMinus().toString());
        System.out.println("-------------------------");
        System.out.println("μ_1 = "+this.getRoot_1().toString());
        System.out.println("μ_2 = "+this.getRoot_2().toString());
        System.out.println("-------------------------");
        System.out.println("Angle="+this.angle1x+" rad");
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++");
    }
}
