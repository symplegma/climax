/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class ElasticMat extends Material{
    private double ElasticModulus;
    private double PoissonRatio;
    private double Lamem;
    private double Lamel;
    private double density;
    private double Area; // used only in rod problems
    private double ThermalExpansionCoef; // also found as dilatation coef. (Brebia & Dominguez 1989)
    
    // constructor
    public ElasticMat(){}
    
    public ElasticMat(int id){
        this.id=id;
    }
    
    public ElasticMat(int id, double density){
        this.id=id;
        this.density=density;
    }
    
    // methods
    public void setE_v(double E, double v){
        this.ElasticModulus=E;
        this.PoissonRatio=v;
        Lamem=E/(2.+2.*v);
        Lamel=E*v/((1.+v)*(1-2.*v));
    }
    
    
    public void setLameM_L(double M, double L){
        this.Lamem=M;
        this.Lamel=L;
        this.ElasticModulus=M*(3.*L+2.*M)/(L+M);
        this.PoissonRatio=L/(2.*(L+M));
    }
    
    public void setLameM_v(double M, double v){
        this.Lamem=M;
        this.PoissonRatio=v;
        this.ElasticModulus=2.*M*(1+v);
        this.Lamel=ElasticModulus*v/((1.+v)*(1-2.*v));
    }
    
    public void setArea(double A){
        this.Area=A;
    }
    
    public void setDensity(double density){
        this.density=density;
    }
    
    public void setThermalExpansionCoef(double ThermalExpansionCoef){
        this.ThermalExpansionCoef=ThermalExpansionCoef;
    }
    
    public double getElasticModulus(){
        return this.ElasticModulus;
    }
    
    public double getPoissonRatio(){
        return this.PoissonRatio;
    }
    
    public double getLame_L(){
        return this.Lamel;
    }
    
    public double getLame_M(){
        return this.Lamem;
    }
    
    public double getDensity(){
        return this.density;
    }
    
    public double getArea(){
        return this.Area;
    }
    
//    public double getThermalExpansionCoef(){
//        return this.ThermalExpansionCoef;
//    }
    
    public double getExtendedThermalCoef(){
        // taken from Brebia & Dominguez second edition p.189
        // also found in Thermal Stress Analysis of Multi-layer Thin Films and Coatings 
        // by an Advanced Boundary Element Method, Xiaolin Chen and Yijun Liu,
        // CMES, vol.2, pp. 337-349, 2001
        return 2.*(1.+this.PoissonRatio)*this.Lamem*this.ThermalExpansionCoef/(1.-2.*this.PoissonRatio);
    }
}
