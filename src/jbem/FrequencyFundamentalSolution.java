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
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;

/**
 *
 * @author pchr
 */
abstract public class FrequencyFundamentalSolution extends FundamentalSolution{
    protected FundamentalSolution RespectiveStatic;
    // constructor
    public FrequencyFundamentalSolution(){}
    
    public AbstractMatrix get_pres_fund() {
        return this.RespectiveStatic.get_p_fund();
    }

    @Override
    public AbstractMatrix get_u_fund() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix get_p_fund() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix get_s_fund() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix get_r_fund() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    abstract public Array2DRowFieldMatrix<Complex> get_u_fundC();
    abstract public Array2DRowFieldMatrix<Complex> get_p_fundC();
    abstract public Array2DRowFieldMatrix<Complex> get_s_fundC();
    abstract public Array2DRowFieldMatrix<Complex> get_r_fundC();
}
