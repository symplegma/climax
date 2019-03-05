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
abstract public class MaterialPoint {
    protected int id;
    protected double[] coords;
    protected double weight;
    protected Material aMaterial;
    static int numOfMaterialPoints=0;
    protected double z=1.;

    // constructor
    public MaterialPoint(){
    }

    public void setCoords(double[] d, double weight){
        this.coords=new double[d.length];
        System.arraycopy(d, 0, coords, 0, d.length);
        this.weight=weight;
    }

    abstract void initMatPoint(int params);

    public void setMaterial(Material aMaterial){
        this.aMaterial=aMaterial;
    }

    abstract void clear();

    abstract void commit();

    public int getID(){
        return this.id;
    }

    public double getGaussCoord(int wcoord){
        return this.coords[wcoord];
    }

    public double getGaussWeight(){
        return this.weight;
    }
    
    public void setDamageVar(double z){this.z=z;}

}
