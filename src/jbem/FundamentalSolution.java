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

package jbem;

import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
abstract public class FundamentalSolution {
    protected int ndofs;
    public static FundamentalSolutionData theFSdata = new FundamentalSolutionData();
    protected boolean isLogFunction=false;
    protected boolean isComplexFunction=false;    
    protected static AbstractMatrix Rot;
    protected int SpaceDimension;
    // constructor 
    public FundamentalSolution(){}
    
    public int getSpaceDimension(){return this.SpaceDimension;}
    
    // abstract methods
    abstract public int get_u_DOFs();
    abstract public int get_v_DOFs();
    abstract public int get_p_DOFs();
    
    abstract public AbstractMatrix get_u_fund();
    abstract public AbstractMatrix get_p_fund();
    abstract public AbstractMatrix get_s_fund();
    abstract public AbstractMatrix get_r_fund();
    
    // methods
    public boolean isLogFunction(){
        return this.isLogFunction;
    }
    
    public boolean isComplexFunction(){
        return this.isComplexFunction;
    }
    

}
