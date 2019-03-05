/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class EnergyDissipator {
    protected int id;
    protected double alphaI;
    protected double alpha2;
    protected double alpha3;
    protected double alpha2_1=0.;
    protected double alpha0=0.;
    protected double epsilon=0.0;
    protected double CoulombCoef;
    protected double viscosityCoef;
    protected int model=0;
    // 0 : constant, 
    // 1 : engineering mixed mode (1.6.2 Chapter Tomas, Book Vlado), it refers to model of: 
    //     J.W. Hutchinson and Z. Suo. Mixed mode cracking in layered materials, volume 29 of Advances in Applied Mechanics. Academic Press: New York, 1992.
    //     implemented for Kruzik+Roubicek (?)
    // 2 : non-associative model based on plasticity
    // - Vodicka R, Mantic V (2011) An energetic approach to mixed mode delamination problem an SGBEM implementation. 
    // - Panagiotopoulos CG, Mantic V, Roubicek T (2012) BEM solution of delamination problems using an interface damage and plasticity model
    //   Computational Mechanics
    // 3 : explicit form of non-associative model based on plasticity
    // - see paper in ECF20, Trondheim, Norway, 2014

    protected double lamdba=1.;
    protected int  fracture_mode_mixity_angle=0;
    // fracture_mode_mixity_angle=0: ψ_G, c=1: ψ_U, c=2: ψ_S

    // constructor
    public EnergyDissipator(){}

    public EnergyDissipator(int id){
        this.id=id;
    }

    public EnergyDissipator(int id, double alphaI){
        this.id=id;
        this.alphaI=alphaI;
    }

    public EnergyDissipator(int id, double alphaI, double alpha2){
        this.id=id;
        this.alphaI=alphaI;
        this.alpha2=alpha2;
    }
    
    public EnergyDissipator(int id, double alphaI, double alpha2, double alpha3){
        this.id=id;
        this.alphaI=alphaI;
        this.alpha2=alpha2;
        this.alpha3=alpha3;
    }

    public void setAlphaI(double v){this.alphaI=v;}

    public void setAlpha2(double v){this.alpha2=v;}
    
    public void setAlpha2_1(double v){this.alpha2_1=v;}
    
    public void setAlpha3(double v){this.alpha3=v;}
    
    public void setAlpha0(double v){this.alpha0=v;}
    
    public void setEpsilon(double v){this.epsilon=v;}

    public double getAlphaI(){return this.alphaI;}

    public double getAlpha2(){return this.alpha2;}
    
    public double getAlpha2_1(){return this.alpha2_1;}
    
    public double getAlpha3(){return this.alpha3;}
    
    public double getAlpha0(){return this.alpha0;}
    
    public double getEpsilon(){return this.epsilon;}
    
    public void setCoulombCoef(double m){this.CoulombCoef=m;}
    
    public double getCoulombCoef(){return this.CoulombCoef;}
    
    public void setViscosityCoef(double m){this.viscosityCoef=m;}
    
    public double getViscosityCoef(){return this.viscosityCoef;}
    
    public void setModel(int mod){this.model=mod;}
    
    public int getModel(){return this.model;}
    
    public void setFractureAngle(int mod){this.fracture_mode_mixity_angle=mod;}
    
    public int getFractureAngle(){return this.fracture_mode_mixity_angle;}
    
    public void setLambda(double val){this.lamdba=val;}
    
    public double getLambda(){return this.lamdba;}
    
    public double f(double phi, double kn, double kt, double kh, double a1, double sy, double tanpsi){
        switch(this.fracture_mode_mixity_angle){
            case 1:
                return Math.tan(phi)*Math.sqrt(kn/kt)*(kt+kh)/kh-Math.sqrt(kn/(2.*a1))*sy/(kh*Math.cos(phi))-tanpsi;
            case 2:
                return (kt/kn)*(Math.tan(phi)*Math.sqrt(kn/kt)*(kt+kh)/kh-Math.sqrt(kn/(2.*a1))*sy/(kh*Math.cos(phi)))-tanpsi;
            default:
                return Math.sqrt(kt/kn)*(Math.tan(phi)*Math.sqrt(kn/kt)*(kt+kh)/kh-Math.sqrt(kn/(2.*a1))*sy/(kh*Math.cos(phi)))-tanpsi;
        }
    }
    
    public double getAlpha(Spring TangentSpring, Spring NormalSpring, int model, double psi){
        double val=0.;
        double sy = TangentSpring.get_YieldStress();
        double kt = TangentSpring.getStiffness();
        double kn = NormalSpring.getStiffness();
        double kh = TangentSpring.getIsotropicHardening();
        switch(model){
            case 2:
                throw new UnsupportedOperationException("Not supported yet.");
            case 3:
                if(psi<=Math.asin(sy/Math.sqrt(2.*kt*alphaI))){
                    val = alphaI;
                }else{
                    val=(2.*alphaI*(kt+kh)-sy*sy)*(1.+Math.tan(psi)*Math.tan(psi))/(2.*(kt+kh+kh*Math.tan(psi)*Math.tan(psi)));
                }
                break;
                
        }
        return val;
    }
}
