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
import jmat.AbstractMatrix;
/**
 *
 * @author pchr
 */
public class Linear extends Algorithm{
    private double tolerance;
    
    // constructor
    public Linear(Analysis theAnalysis, double tolerance){
        this.theAnalysis=theAnalysis;
        this.tolerance=tolerance;
        this.theType=1;
    }
    
    public Linear(double tolerance){
        this.tolerance=tolerance;
        this.theType=1;
    }
    
    public Linear(){
        this.tolerance=jmat.MachinePrecision.getMachinePrecision();
        this.theType=1;
    }
    
    //methods
    public int solve(){
        int log=0;

        //this.theAnalysis.theSOE.getResidual();
        int order=theAnalysis.theDomain.getNdofs();
        double[] Res = new double[order];
        AbstractMatrix r= new AbstractMatrix(order,1);

        theAnalysis.run();
        theAnalysis.update();
        theAnalysis.formRight();
        Res= theAnalysis.theSOE.getResidual();
        for(int j=0; j<Res.length; j++){
            r.addVal(j, 0, Res[j]);
        }
        if(r.norm1()/theAnalysis.theSOE.getNormFext()>tolerance){log=1;}

        theAnalysis.commit();
        return log;
    }
    
    public void setAnalysis(Analysis theAnalysis){
        this.theAnalysis=theAnalysis;
    }

}
