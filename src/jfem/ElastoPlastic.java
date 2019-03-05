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
abstract public class ElastoPlastic extends ElasticMaterial{
    protected double  Kisotropic;
    protected double  Kkinematic;
    protected double[]  YieldStress;

    public ElastoPlastic(){}

    // methods that subclasses are going to produce
    public void setKkinematic(double Kkinematic){
    }

    public void setKisotropic(double Kisotropic){
    }

    //public void setYieldStress(double YieldStress){
    //}

    public void setYieldStress(double[] YieldStress){
    }

    public double getKisotropic(){
        return 0.;
    }

    public double getKkinematic(){
        return 0.;
    }

    abstract public double[] getYieldStress();
}
