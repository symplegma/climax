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
public class StaticElasticity3DFS extends FundamentalSolution{
    
    // constructor
    public StaticElasticity3DFS(){
        this.ndofs=3;
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
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity3DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity3DFS.theFSdata.getR();
        double r=StaticElasticity3DFS.theFSdata.getAbsR();
        
        double[] dr=StaticElasticity3DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double c=1./( 16.0*Math.PI*anElasticMat.getLame_M()*(1.-anElasticMat.getPoissonRatio()) );
        double c1=3.-4.*anElasticMat.getPoissonRatio();
        
        F.init();
        
        F.set(0, 0, c/r*(c1+dr[0]*dr[0]));
        F.set(0, 1, c/r*(dr[0]*dr[1]));
        F.set(0, 2, c/r*(dr[0]*dr[2]));
        
        F.set(1, 0, c/r*(dr[0]*dr[1]));
        F.set(1, 1, c/r*(c1+dr[1]*dr[1]));
        F.set(1, 2, c/r*(dr[1]*dr[2]));
        
        F.set(2, 0, c/r*(dr[0]*dr[2]));
        F.set(2, 1, c/r*(dr[1]*dr[2]));
        F.set(2, 2, c/r*(c1+dr[2]*dr[2]));
        
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity3DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity3DFS.theFSdata.getR();
        double r=StaticElasticity3DFS.theFSdata.getAbsR();
        
        double[] dr=StaticElasticity3DFS.theFSdata.getDR();
        double[] n=StaticElasticity3DFS.theFSdata.getOutwardNormal();
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double c2=1./( 8.0*Math.PI*(1.-anElasticMat.getPoissonRatio()) );
        double c3=1.-2.*anElasticMat.getPoissonRatio();
        //double G=anElasticMat.getLame_M();
        
        F.init();
        /*
        F.set(0, 0, -c2/(r*r)*(c3+3.0*dr[0]*dr[0])*costh);
        F.set(0, 1, -c2/(r*r)*(3.0*dr[0]*dr[1]*costh-c3*(n[1]*dr[0]-n[0]*dr[1])));
        F.set(0, 2, -c2/(r*r)*(3.0*dr[0]*dr[2]*costh-c3*(n[2]*dr[0]-n[0]*dr[2])));
        
        F.set(1, 0, -c2/(r*r)*(3.0*dr[0]*dr[1]*costh-c3*(n[0]*dr[1]-n[1]*dr[0])));
        F.set(1, 1, -c2/(r*r)*(c3+3.0*dr[1]*dr[1])*costh);
        F.set(1, 2, -c2/(r*r)*(3.0*dr[1]*dr[2]*costh-c3*(n[2]*dr[1]-n[1]*dr[2])));
        
        F.set(2, 0, -c2/(r*r)*(3.0*dr[0]*dr[2]*costh-c3*(n[0]*dr[2]-n[2]*dr[0])));
        F.set(2, 1, -c2/(r*r)*(3.0*dr[1]*dr[2]*costh-c3*(n[1]*dr[2]-n[2]*dr[1])));
        F.set(2, 2, -c2/(r*r)*(c3+3.0*dr[2]*dr[2])*costh);
        */
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                F.set(i, j, 
                        -c2/(r*r)*(c3*(n[i]*dr[j]-n[j]*dr[i])+(3.*dr[i]*dr[j]+c3*delta(i,j))*costh)
                        );
            }
        }
        
        return F;
    }
    
    private double delta(int i, int j){
        double d=0.;
        if(i==j){
            d=1.;
        }
        return d;
    }

    @Override
    public AbstractMatrix get_s_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity3DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity3DFS.theFSdata.getR();
        double r=StaticElasticity3DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        double G=anElasticMat.getLame_M();
        double d=2.;
        double C2=1./(8.*Math.PI*(1.-v));
        double C3=1.-2.*v;
        double C5=G/(4.*Math.PI*(1.-v));
        double C6=5.;
        double C7=1.-4.*v;
        
        double[] dr=StaticElasticity3DFS.theFSdata.getDR();
        double[] n=StaticElasticity3DFS.theFSdata.getOutwardNormal();
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,2*this.ndofs);
        
        F.init();
        int i,k,j;
        int x,y,z; x=0; y=1; z=2;
        
        i=x;j=x;k=x;
        F.set(0, 0, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxxx
        i=x;j=x;k=y;
        F.set(1, 0, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxxy
        i=x;j=x;k=z;
        F.set(2, 0, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxxz
        
        i=y;j=y;k=x;
        F.set(0, 1, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Syyx
        i=y;j=y;k=y;
        F.set(1, 1, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Syyy
        i=y;j=y;k=z;
        F.set(2, 1, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Syyz
        
        i=z;j=z;k=x;
        F.set(0, 2, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Szzx
        i=z;j=z;k=y;
        F.set(1, 2, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Szzy
        i=z;j=z;k=z;
        F.set(2, 2, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Szzz
        
        i=x;j=y;k=x;
        F.set(0, 3, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxyx
        i=x;j=y;k=y;
        F.set(1, 3, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxyy
        i=x;j=y;k=z;
        F.set(2, 3, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxyz
        
        i=y;j=z;k=x;
        F.set(0, 4, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Syzx
        i=y;j=z;k=y;
        F.set(1, 4, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Syzy
        i=y;j=z;k=z;
        F.set(2, 4, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Syzz
        
        i=x;j=z;k=x;
        F.set(0, 5, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxzx
        i=x;j=z;k=y;
        F.set(1, 5, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxzy
        i=x;j=z;k=z;
        F.set(2, 5, (C2/(r*r))*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+3.*dr[i]*dr[j]*dr[k])); // Sxzz
        
        return F;
    }

    @Override
    public AbstractMatrix get_r_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity3DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity3DFS.theFSdata.getR();
        double r=StaticElasticity3DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        double G=anElasticMat.getLame_M();
        double d=2.;
        double C2=1./(8.*Math.PI*(1.-v));
        double C3=1.-2.*v;
        double C5=G/(4.*Math.PI*(1.-v));
        double C6=5.;
        double C7=1.-4.*v;
        
        double[] dr=StaticElasticity2DFS.theFSdata.getDR();
        double[] n=StaticElasticity2DFS.theFSdata.getOutwardNormal();
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,2*this.ndofs);
        
        F.init();
        
        int i,k,j;
        int x,y,z; x=0; y=1; z=2;
        
        i=x;j=x;k=x;
        F.set(0, 0, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxxx
        i=x;j=x;k=y;
        F.set(1, 0, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxxy
        
        i=x;j=x;k=z;
        F.set(2, 0, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxxz
        //-----------------------------------------------------------------------------------------------
        i=y;j=y;k=x;
        F.set(0, 1, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryyx
        i=y;j=y;k=y;
        F.set(1, 1, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryyy
        
        i=y;j=y;k=z;
        F.set(2, 1, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryyz
        //-----------------------------------------------------------------------------------------------
        i=z;j=z;k=x;
        F.set(0, 2, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rzzx
        
        i=z;j=z;k=y;
        F.set(1, 2, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rzzy
        
        i=z;j=z;k=z;
        F.set(2, 2, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rzzz
        //-----------------------------------------------------------------------------------------------
        i=x;j=y;k=x;
        F.set(0, 3, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxyx
        i=x;j=y;k=y;
        F.set(1, 3, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxyy
        
        i=x;j=y;k=z;
        F.set(2, 3, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxyz
        //-----------------------------------------------------------------------------------------------
        i=y;j=z;k=x;
        F.set(0, 4, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryzx
        i=y;j=z;k=y;
        F.set(1, 4, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryzy
        
        i=y;j=z;k=z;
        F.set(2, 4, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryzz
        //-----------------------------------------------------------------------------------------------
        i=x;j=z;k=x;
        F.set(0, 5, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxzx
        i=x;j=z;k=y;
        F.set(1, 5, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxzy
        
        i=x;j=z;k=z;
        F.set(2, 5, C5/(r*r*r)*(3.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +3.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(3.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxzz
        //-----------------------------------------------------------------------------------------------
        return F;
    }

}
