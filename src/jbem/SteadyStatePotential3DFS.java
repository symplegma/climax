/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class SteadyStatePotential3DFS extends FundamentalSolution{
    
    // constructor
    public SteadyStatePotential3DFS(){
        this.ndofs=1;
        SpaceDimension=3;
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
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential3DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential3DFS.theFSdata.getR();
        double r=SteadyStatePotential3DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential3DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        double k = theMaterial.getConductivity();
        
        F.set(0, 0, 1./(4.*Math.PI*k*r));
        
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        PotentialMat theMaterial = (PotentialMat) SteadyStatePotential3DFS.theFSdata.getMaterial();
        double[] R=SteadyStatePotential3DFS.theFSdata.getR();
        double r=SteadyStatePotential3DFS.theFSdata.getAbsR();
        
        double[] dr=SteadyStatePotential3DFS.theFSdata.getDR();
        double[] n=SteadyStatePotential3DFS.theFSdata.getOutwardNormal();
        
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        
        F.set(0, 0, -costh/(4.*Math.PI*r*r));
        
        return F;
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
