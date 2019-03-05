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
public class StaticRodFS extends FundamentalSolution{
    
    // constructor
    public StaticRodFS(){
        this.ndofs=1;
        SpaceDimension=1;
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
    public int get_p_DOFs(){ 
        return this.ndofs;
    }

    @Override
    public AbstractMatrix get_u_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticRodFS.theFSdata.getMaterial();
        double[] R=StaticRodFS.theFSdata.getR();
        double r=0.0;
        for (int i=0; i<R.length; i++){
            r+=R[i]*R[i];
        }
        r=Math.sqrt(r);
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double E=anElasticMat.getElasticModulus();
        double A=anElasticMat.getArea();
        
        F.init();
        
        F.set(0, 0, r/(2.*A*E));
        
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticRodFS.theFSdata.getMaterial();
        double[] R=StaticRodFS.theFSdata.getR();
        double[] n=StaticElasticity3DFS.theFSdata.getOutwardNormal();
        double r=0.0;
        for (int i=0; i<R.length; i++){
            r+=R[i]*R[i];
        }
        r=Math.sqrt(r);
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        
        F.init();
        
        //F.set(0, 0, r*n[0]/(2.*R[0]));
        if(R[0]!=0.){
            F.set(0, 0, r/(2.*R[0]));
        }else{
            F.set(0, 0, r*n[0]/(2.*R[0]));
        }
        
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
