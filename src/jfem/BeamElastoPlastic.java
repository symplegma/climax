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
public class BeamElastoPlastic extends ElastoPlastic{
    private double  coef=0.;
    
    // constructor
    public BeamElastoPlastic(int id){
        this.id =id;
        this.type=1; // means that material is UniaxialElastoPlastic
        YieldStress = new double[2];
        ++numberOfMaterials;
    }
    
    public BeamElastoPlastic(int id, double Elasticity){
        this.id =id;
        this.type=1; // means that material is UniaxialElastoPlastic
        this.Elasticity=Elasticity;
        YieldStress = new double[2];
        ++numberOfMaterials;
    }
    
    public BeamElastoPlastic(int id, double Elasticity, double Poisson){
        this.id =id;
        this.type=1; // means that material is UniaxialElastoPlastic
        this.Elasticity=Elasticity;
        this.Poisson=Poisson;
        YieldStress = new double[2];
        ++numberOfMaterials;
    }
    
    public BeamElastoPlastic(int id, double Elasticity, double Poisson, double density){
        this.id =id;
        this.type=1; // means that material is UniaxialElastoPlastic
        this.Elasticity=Elasticity;
        this.Poisson=Poisson;
        this.density=density;
        YieldStress = new double[2];
        ++numberOfMaterials;
    }
    
    // methods
    @Override
    public void setKisotropic(double Kisotropic){
        this.Kisotropic=Kisotropic;
    }
    
    @Override
    public void setKkinematic(double Kkinematic){
        this.Kkinematic=Kkinematic;
    }
    
    @Override
    public void setYieldStress(double[] YieldStress){
        this.YieldStress[0]=YieldStress[0];
        this.YieldStress[1]=YieldStress[1];
    }
    
    public void setC(double c){
        this.coef=c;
    }
    
    public double[] getYieldStress(){
        return this.YieldStress;
    }
    
    public double getC(){
        return this.coef;
    }

}
