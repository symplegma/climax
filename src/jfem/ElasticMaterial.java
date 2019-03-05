/*******************************************************************************
* Climax.                                                                      *
* Copyright (C) 2009-2017 C.G. Panagiotopoulos [http://www.symplegma.org]      *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): C.G. Panagiotopoulos (pchr76@gmail.com)
// *****************************************************************************

package jfem;

/**
 *
 * @author pchr
 */
public class ElasticMaterial extends Material{
    protected double Elasticity;
    protected double Poisson;
    protected double density;
    protected double FracTough;
    protected double DamageStress=0.0;

    public ElasticMaterial(){}

    public ElasticMaterial(int id, double Elasticity){
        this.id =id;
        this.Elasticity=Elasticity;
        type=0;
        ++numberOfMaterials;
    }

    public ElasticMaterial(int id, double Elasticity, double Poisson){
        this.id =id;
        this.Elasticity=Elasticity;
        this.Poisson=Poisson;
        type=0;
        ++numberOfMaterials;
    }

    public ElasticMaterial(int id, double Elasticity, double Poisson, double density){
        this.id =id;
        this.Elasticity=Elasticity;
        this.Poisson=Poisson;
        this.density=density;
        type=0;
        ++numberOfMaterials;
    }

    public double getElasticity(){
        return this.Elasticity;
    }

    public double getG(){
        double G=this.Elasticity/(2.+2.*this.Poisson);
        return G;
    }

    public double getDensity(){
        return this.density;
    }

    public double getPoisson(){
        return this.Poisson;
    }

    public void setElasticity(double Elasticity){
        this.Elasticity=Elasticity;
    }

    public void setDensity(double density){
        this.density=density;
    }

    public void setPoisson(double Poisson){
        this.Poisson=Poisson;
    }

    public void setFracTough(double FracTough){
        this.FracTough=FracTough;
    }

    public double getFracTough(){
        return this.FracTough;
    }
    
    public double getLameLambda(){
        return Elasticity*Poisson/((1.+Poisson)*(1-2.*Poisson));
    }
    
    public void setDamageStress(double DS){this.DamageStress=DS;}

    public double getDamageStress(){
        return this.DamageStress;
    }
}

