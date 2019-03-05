/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jfem;

/**
 *
 * @author pchr
 */
public class AcousticMaterial extends Material{
    protected double MatConstant=0.0;

    public AcousticMaterial(){
        type=4;
    }

    public AcousticMaterial(int id, double Elasticity){
        type=4;
        this.id =id;
        this.MatConstant=Elasticity;
        ++numberOfMaterials;
    }
    
    public double getMatConstant(){return this.MatConstant;}
}
