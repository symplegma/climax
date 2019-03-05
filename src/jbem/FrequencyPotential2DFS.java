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

import static mathman.BesselK.BesselKs;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;

/**
 *
 * @author pchr
 */
public class FrequencyPotential2DFS extends FrequencyFundamentalSolution{
    
    // constructor
    public FrequencyPotential2DFS(){
        this.isComplexFunction=true;
        this.isLogFunction=true;
        this.ndofs=1;
        RespectiveStatic = new SteadyStatePotential2DFS();
        SpaceDimension=2;
    }

    @Override
    public int get_u_DOFs() {
        return this.ndofs;
    }

    @Override
    public int get_v_DOFs() {
        return 0;
    }

    @Override
    public int get_p_DOFs() {
        return this.ndofs;
    }

    @Override
    public Array2DRowFieldMatrix<Complex> get_u_fundC() {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential2DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential2DFS.theFSdata.getR();
        double r=SteadyStatePotential2DFS.theFSdata.getAbsR();
        double omg=SteadyStatePotential2DFS.theFSdata.getCurrentFrequency();
        double k = theMaterial.getConductivity();
        
        double vlc=theMaterial.getPropagationVelc();
        Complex[][] dmat = new Complex[ndofs][ndofs];
        
        Complex cdat= new Complex(0.0,omg*r/vlc);
        cdat=BesselKs(cdat)[0];cdat=cdat.divide(2.*Math.PI*k);
        dmat[0][0]=new Complex(cdat.getReal(),cdat.getImaginary());
        Array2DRowFieldMatrix<Complex>  F = new Array2DRowFieldMatrix<Complex>(dmat);
        return F;
    }

    @Override
    public Array2DRowFieldMatrix<Complex> get_p_fundC() {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential2DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential2DFS.theFSdata.getR();
        double r=SteadyStatePotential2DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential2DFS.theFSdata.getDR();
        double[] n=SteadyStatePotential2DFS.theFSdata.getOutwardNormal();
        double omg=SteadyStatePotential2DFS.theFSdata.getCurrentFrequency();
        double k = theMaterial.getConductivity();
        
        double vlc=theMaterial.getPropagationVelc();
        Complex[][] dmat = new Complex[ndofs][ndofs];
        
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        Complex cdat= new Complex(0.0,omg*r/vlc);
        cdat=BesselKs(cdat)[1];cdat=cdat.multiply(new Complex(0.0,-omg*costh/(2.*Math.PI*vlc)));
        dmat[0][0]=new Complex(cdat.getReal(),cdat.getImaginary());
        Array2DRowFieldMatrix<Complex>  F = new Array2DRowFieldMatrix<Complex>(dmat);
        
        //F.set(0, 0, -costh/(2.*Math.PI*r));
        return F;
    }

    @Override
    public Array2DRowFieldMatrix<Complex> get_s_fundC() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Array2DRowFieldMatrix<Complex> get_r_fundC() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
