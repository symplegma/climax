/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class Spring {
    protected int id;
    protected double stiffness;
    protected double IsotropicHardening;
    protected double gradientZ;
    protected double r;
    protected double stiffness_0=0.0;
    protected double gradientPlasticity;
    protected double DamageViscosity=0.0;
    protected Double CriticalStress=null;
    protected Double YieldStress=null;

    // constructor
    public Spring(){}

    public Spring(int id, double stiffness){
        this.id=id;
        this.stiffness=stiffness;
        this.IsotropicHardening=0.;
    }

    // methods
    public int getID(){
        return this.id;
    }

    public double getStiffness(){
        return this.stiffness;
    }
    
    public double getStiffness_0(){
        return this.stiffness_0;
    }

    public void setIsotropicHardening(double v){
        this.IsotropicHardening=v;
    }
    
    public double getIsotropicHardening(){return this.IsotropicHardening;}
    
    public void setGradientZ(double k0){this.gradientZ=k0;}
    
    public void setStiffness_0(double k0){this.stiffness_0=k0;}
    
    public void setGradient_r(double r){this.r=r;}
    
    public void setGradientZ_r(double k0, double r){this.gradientZ=k0; this.r=r;}
    
    public double get_r(){return this.r;}
    
    public double get_GradientZ(){return this.gradientZ;}
    
    public void setGradientPlasticity(double e){this.gradientPlasticity=e;}
    
    public double get_GradientPlasticity(){return this.gradientPlasticity;}
    
    public void setDamageViscosity(double e){this.DamageViscosity=e;}
    
    public double get_DamageViscosity(){return this.DamageViscosity;}
    
    public void setCriticalStress(double e){this.CriticalStress=e;}
    
    public double get_CriticalStress(){return this.CriticalStress;}
    
    public void setYieldStress(double e){this.YieldStress=e;}
    
    public double get_YieldStress(){return this.YieldStress;}

}
