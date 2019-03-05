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
public class InitialNewtonRaphson extends Algorithm{
    private double tolerance;
    private int maxSteps;
    
    // constructor
    public InitialNewtonRaphson(Analysis theAnalysis, double tolerance,int maxSteps){
        this.theAnalysis=theAnalysis;
        this.tolerance=tolerance;
        this.maxSteps=maxSteps;
        this.theType=4;
    }
    
    public InitialNewtonRaphson(double tolerance,int maxSteps){
        this.tolerance=tolerance;
        this.maxSteps=maxSteps;
        this.theType=4;
    }
    
    //methods
    public int solve(){
        int log=0;
        
        //this.theAnalysis.theSOE.getResidual();
        int order=theAnalysis.theDomain.getNdofs();
        double[] Res = new double[order];
        Res= theAnalysis.theSOE.getResidual();
        AbstractMatrix r= new AbstractMatrix(order,1);
        for(int j=0; j<Res.length; j++){
            r.addVal(j, 0, Res[j]);
        }
        
        iterations=0;
        while( (r.norm1()/theAnalysis.theSOE.getNormFext()>tolerance) && (iterations<maxSteps)  ){
            iterations++;
            theAnalysis.run();
            theAnalysis.update();
            theAnalysis.formRight();
            Res= theAnalysis.theSOE.getResidual();
            for(int j=0; j<Res.length; j++){
                r.set(j, 0, Res[j]);
            }
            if(iterations>=maxSteps&&r.norm1()/theAnalysis.theSOE.getNormFext()>tolerance){log=1;}
        }

        theAnalysis.commit();
        return log;
    }
    
        public void setAnalysis(Analysis theAnalysis){
        this.theAnalysis=theAnalysis;
    }
}
