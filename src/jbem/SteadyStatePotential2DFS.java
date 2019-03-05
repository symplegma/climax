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
public class SteadyStatePotential2DFS extends FundamentalSolution implements LogFunction{
    
    // constructor
    public SteadyStatePotential2DFS(){
        this.ndofs=1;
        this.isLogFunction=true;
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
    public AbstractMatrix get_u_fund() {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential2DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential2DFS.theFSdata.getR();
        double r=SteadyStatePotential2DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential2DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        double k = theMaterial.getConductivity();
        
        F.set(0, 0, Math.log(1./r)/(2.*Math.PI*k));
        
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential2DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential2DFS.theFSdata.getR();
        double r=SteadyStatePotential2DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential2DFS.theFSdata.getDR();
        double[] n=SteadyStatePotential2DFS.theFSdata.getOutwardNormal();
        
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        
        F.set(0, 0, -costh/(2.*Math.PI*r));
        
        return F;
    }

    public AbstractMatrix getuLogPart() {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential2DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential2DFS.theFSdata.getR();
        double r=SteadyStatePotential2DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential2DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        double k = theMaterial.getConductivity();
        
        F.set(0, 0, 1./(2.*Math.PI*k));
        
        return F;
    }

    public AbstractMatrix getuNonLogPart(double val) {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential2DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential2DFS.theFSdata.getR();
        double r=SteadyStatePotential2DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential2DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        double k = theMaterial.getConductivity();
        
        F.set(0, 0, val/(2.*Math.PI*k));
        
        return F;
    }

    public AbstractMatrix getvLogPart() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public AbstractMatrix getvNonLogPart() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_s_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_r_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
