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
import java.util.Iterator;
import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
abstract public class PenaltyImposer extends BCImposer {
    
    // constructor 
    public PenaltyImposer(){
    }
    
    protected void imposeK(Analysis theAnalysis){
        for(Iterator<ConstraintElement> it=theAnalysis.getDomain().getConstraintElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getK();
            int[] aInt ;
            aInt=elem.getFtable();
            theAnalysis.getSOE().addToK(mat, aInt, coef_k);
        }
    }
    
    protected void imposeM(Analysis theAnalysis){
        for(Iterator<ConstraintElement> it=theAnalysis.getDomain().getConstraintElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM();
            int[] aInt ;
            aInt=elem.getFtable();
            theAnalysis.getSOE().addToM(mat, aInt, coef_m);
        }
    }

    void imposeF_ext(Analysis theAnalysis, int time){
        Domain aDomain=theAnalysis.getDomain();
        SOE theSOE=theAnalysis.getSOE();
        for(Iterator<ConstraintElement> it=aDomain.getConstraintElements().values().iterator(); it.hasNext();){
            ConstraintElement aConstraint = it.next();
            AbstractMatrix mat ;
            mat=aConstraint.getFloading_stif(time);
            int[] aInt ;
            aInt=aConstraint.getFtable();
            theSOE.addToFext(mat, aInt, 1.0);
        }
    }
    
    void imposeF_int(Analysis theAnalysis){
        Domain aDomain=theAnalysis.getDomain();
        SOE theSOE=theAnalysis.getSOE();
        for(Iterator<ConstraintElement> it=aDomain.getConstraintElements().values().iterator(); it.hasNext();){
            ConstraintElement aConstraint = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=aConstraint.getFtable();
            mat=aConstraint.getF();
            theSOE.addToFint(mat, aInt, 1.0);
        }
    }
    
    @Override
    abstract void imposeLeft(Analysis theAnalysis) ;
}

