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
abstract public class BCImposer {
    protected int TypeOfImposer;
    protected double coef_k;
    protected double coef_m;
    
    // constructor
    public BCImposer(){
        this.coef_k=1.0;
        this.coef_m=1.0;
    }
    
    // methods
    public int getType(){
        return TypeOfImposer;
    }
    
    public void setCoefficient_k(double val){
        this.coef_k=val;
    }
    
    public double getCoefficient_k(){
        return this.coef_k;
    }
    
    public void setCoefficient_m(double val){
        this.coef_m=val;
    }
    
    public double getCoefficient_m(){
        return this.coef_m;
    }
    
    abstract void imposeLeft(Analysis theAnalysis);
    abstract void imposeF_ext(Analysis theAnalysis, int time);
    abstract void imposeF_int(Analysis theAnalysis);
}
