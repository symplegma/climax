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
public class MaterialElasticPlasticPoint extends MaterialPoint{
    boolean isPlastified=false;
    private double[] Stress;
    private double[]  Eplastic;
    private double[]  Etotal;
    private double  Xsi;
    private double  YieldFunction;
    private double[]  aPar;
    private double[]  qPar;

    private double[]  EplasticTrial;
    private double[]  aParTrial;
    private double[]  qParTrial;
    private double[]  EplasticIncr;

    public MaterialElasticPlasticPoint(){
        this.id =++numOfMaterialPoints;
    }

    @Override
    void initMatPoint(int params) {
        Stress = new double[params];
        Eplastic = new double[params];
        Etotal = new double[params];
        aPar = new double[params];
        qPar = new double[params];

        EplasticTrial = new double[params];
        aParTrial = new double[params];
        qParTrial = new double[params];
        
        EplasticIncr = new double[params];
    }

    public boolean AskIfIsPlastified(){
        return this.isPlastified;
    }

    public double[] giveStrainTakeStress(double[] Etotal){
        for(int i=0;i<Stress.length;i++){
            if(this.aMaterial.getType()==0){
                Stress[i]=( (ElasticMaterial) aMaterial).getElasticity()*Etotal[i];
            }else{
                Stress[i]=( (ElastoPlastic) aMaterial).getElasticity()*(Etotal[i]-Eplastic[i]);
                this.Xsi=this.Stress[i]-qPar[i];
                YieldFunction=Math.abs(this.Xsi)
                                  - (((ElastoPlastic) aMaterial).getYieldStress()[0]
                                    +((ElastoPlastic) aMaterial).getKisotropic()*this.aPar[i]);
                //this.isPlastified=false;
                if(this.YieldFunction>0.){ //then return mapping from Simo & Hughes pp.45 Box 1.5
                    this.isPlastified=true;
                    double Dg=(this.YieldFunction)/( ((ElastoPlastic) aMaterial).getElasticity()*z
                                                    +((ElastoPlastic) aMaterial).getKisotropic()
                                                    +((ElastoPlastic) aMaterial).getKkinematic());
                    this.Stress[i]-=Dg*((ElastoPlastic) aMaterial).getElasticity()*z*Math.signum(this.Xsi);
                    this.EplasticTrial[i]=this.Eplastic[i]+Dg*Math.signum(this.Xsi);
                    this.qParTrial[i]=this.qPar[i]+Dg*((ElastoPlastic) aMaterial).getKkinematic()*Math.signum(this.Xsi);
                    this.aParTrial[i]=this.aPar[i]+Dg;
                }

            }
        }
        return this.Stress;
    }
    
    public double[] giveStraingetPlasticStrain(double[] Etotal){
        double[] val = null;
        val = new double[Etotal.length];
        for(int i=0;i<Etotal.length;i++){
            if(this.aMaterial.getType()==0){
                val[i]=0.0;
            }else{
                Stress[i]=( (ElastoPlastic) aMaterial).getElasticity()*(Etotal[i]-Eplastic[i]);
                this.Xsi=this.Stress[i]-qPar[i];
                YieldFunction=Math.abs(this.Xsi)
                                  - (((ElastoPlastic) aMaterial).getYieldStress()[0]
                                    +((ElastoPlastic) aMaterial).getKisotropic()*this.aPar[i]);
                //this.isPlastified=false;
                if(this.YieldFunction>0.){ //then return mapping from Simo & Hughes pp.45 Box 1.5
                    double Dg=(this.YieldFunction)/( ((ElastoPlastic) aMaterial).getElasticity()*z
                                                    +((ElastoPlastic) aMaterial).getKisotropic()
                                                    +((ElastoPlastic) aMaterial).getKkinematic());
                    this.Stress[i]-=Dg*((ElastoPlastic) aMaterial).getElasticity()*z*Math.signum(this.Xsi);
                    val[i]=this.Eplastic[i]+Dg*Math.signum(this.Xsi);
                }
            }
        }
        return val;
    }
    
    public double[] getPlasticStrain(){return this.Eplastic;}
    
    public double[] getPlasticStrainInc(){return this.EplasticIncr;}

    @Override
    void clear() {
        for(int i=0; i<Eplastic.length;i++){
            Eplastic[i]=0.;
            Etotal[i]=0.;
            Stress[i]=0.;
            aPar[i]=0.;
            qPar[i]=0.;
            EplasticTrial[i]=0.;
            aParTrial[i]=0.;
            qParTrial[i]=0.;
        }

        this.Xsi=0.;
        this.YieldFunction=0.;
        this.isPlastified=false;
    }

    @Override
    void commit(){
        for(int i=0; i<Eplastic.length;i++){
            EplasticIncr[i]=EplasticTrial[i]-Eplastic[i];
            Eplastic[i]=EplasticTrial[i];
            aPar[i]=aParTrial[i];
            qPar[i]=qParTrial[i];
        }
    }
}
