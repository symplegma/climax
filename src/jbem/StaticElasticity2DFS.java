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
public class StaticElasticity2DFS extends FundamentalSolution implements LogFunction{
    private boolean PlaneStress=false;
    
    // constructor
    public StaticElasticity2DFS(){
        this.ndofs=2;
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
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity2DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity2DFS.theFSdata.getR();
        double r=StaticElasticity2DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        if(this.PlaneStress)v=v/(1.+v);
        double angle = FundamentalSolution.theFSdata.getMaterial().getAngle();
        double cos=Math.cos(angle); double sin=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, cos); Rot.set(0, 1, -sin);
        Rot.set(1, 0, sin); Rot.set(1, 1, cos);
        
        double[] dr=StaticElasticity2DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double c=1./( 8.0*Math.PI*anElasticMat.getLame_M()*(1.-v) );
        double c1=3.-4.*v;
        
        F.init();
        
        F.set(0, 0, c*(c1*Math.log(1./r)+dr[0]*dr[0]));
        F.set(0, 1, c*dr[0]*dr[1]);
        
        F.set(1, 0, c*dr[0]*dr[1]);
        F.set(1, 1, c*(c1*Math.log(1./r)+dr[1]*dr[1]));
        
//        F=Rot.times(F.times(Rot.transpose()));
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity2DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity2DFS.theFSdata.getR();
        double r=StaticElasticity2DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        if(this.PlaneStress)v=v/(1.+v);
        
        double[] dr=StaticElasticity2DFS.theFSdata.getDR();
        double[] n=StaticElasticity2DFS.theFSdata.getOutwardNormal();
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double angle = FundamentalSolution.theFSdata.getMaterial().getAngle();
        double cos=Math.cos(angle); double sin=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, cos); Rot.set(0, 1, -sin);
        Rot.set(1, 0, sin); Rot.set(1, 1, cos);
        
        double c2=1./( 4.0*Math.PI*(1.-v) );
        double c3=1.-2.*v;
        
        F.init();
        
        F.set(0, 0, -c2/r*costh*(c3+2.0*dr[0]*dr[0]));
        F.set(0, 1, -c2/r*(2.0*dr[0]*dr[1]*costh-c3*(n[1]*dr[0]-n[0]*dr[1])));
        
        F.set(1, 0,  -c2/r*(2.0*dr[0]*dr[1]*costh-c3*(n[0]*dr[1]-n[1]*dr[0])));
        F.set(1, 1, -c2/r*costh*(c3+2.0*dr[1]*dr[1]));
//        F=Rot.times(F.times(Rot.transpose()));
        return F;
    }
    
    @Override
    public AbstractMatrix get_s_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity2DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity2DFS.theFSdata.getR();
        double r=StaticElasticity2DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        double G=anElasticMat.getLame_M();
        double d=1.;
        double C2=1./(4.*Math.PI*(1.-v));
        double C3=1.-2.*v;
        double C5=G/(2.*Math.PI*(1.-v));
        double C6=4.;
        double C7=1.-4.*v;
        if(this.PlaneStress){
            d=1.;
            C2=(1.+v)/(4.*Math.PI);
            C3=(1.-v)/(1.+v);
            C5=(1.+v)*G/(2.*Math.PI);
            C6=4.;
            C7=(1.-3*v)/(1.+v);
        }
        
        double[] dr=StaticElasticity3DFS.theFSdata.getDR();
        double[] n=StaticElasticity3DFS.theFSdata.getOutwardNormal();
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs+1);
        
        F.init();
        int i,k,j;
        int x,y; x=0; y=1;
        
        i=x;j=x;k=x;
        F.set(0, 0, C2/r*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+2.*dr[i]*dr[j]*dr[k])); // Sxxx
        i=x;j=x;k=y;
        F.set(1, 0, C2/r*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+2.*dr[i]*dr[j]*dr[k])); // Sxxy
        
        i=y;j=y;k=x;
        F.set(0, 1, C2/r*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+2.*dr[i]*dr[j]*dr[k])); // Syyx
        i=y;j=y;k=y;
        F.set(1, 1, C2/r*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+2.*dr[i]*dr[j]*dr[k])); // Syyy
        
        i=x;j=y;k=x;
        F.set(0, 2, C2/r*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+2.*dr[i]*dr[j]*dr[k])); // Sxyx
        i=x;j=y;k=y;
        F.set(1, 2, C2/r*(C3*(delta(k,i)*dr[j]+delta(k,j)*dr[i]-delta(i,j)*dr[k])+2.*dr[i]*dr[j]*dr[k])); // Sxyy
        
        return F;
    }
    
    @Override
    public AbstractMatrix get_r_fund() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity2DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity2DFS.theFSdata.getR();
        double r=StaticElasticity2DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        double G=anElasticMat.getLame_M();
        double d=1.;
        double C2=1./(4.*Math.PI*(1.-v));
        double C3=1.-2.*v;
        double C5=G/(2.*Math.PI*(1.-v));
        double C6=4.;
        double C7=1.-4.*v;
        if(this.PlaneStress){
            d=1.;
            C2=(1.+v)/(4.*Math.PI);
            C3=(1.-v)/(1.+v);
            C5=(1.+v)*G/(2.*Math.PI);
            C6=4.;
            C7=(1.-3*v)/(1.+v);
        }
        
        double[] dr=StaticElasticity2DFS.theFSdata.getDR();
        double[] n=StaticElasticity2DFS.theFSdata.getOutwardNormal();
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs+1);
        
        F.init();
        
        int i,k,j;
        int x,y; x=0; y=1;
        
        i=x;j=x;k=x;
        F.set(0, 0, C5/(r*r)*(2.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +2.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(2.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxxx
        i=x;j=x;k=y;
        F.set(1, 0, C5/(r*r)*(2.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +2.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(2.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxxy
        
        i=y;j=y;k=x;
        F.set(0, 1, C5/(r*r)*(2.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +2.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(2.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryyx
        i=y;j=y;k=y;
        F.set(1, 1, C5/(r*r)*(2.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +2.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(2.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Ryyy
        
        i=x;j=y;k=x;
        F.set(0, 2, C5/(r*r)*(2.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +2.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(2.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxyx
        i=x;j=y;k=y;
        F.set(1, 2, C5/(r*r)*(2.*costh*(C3*delta(i,j)*dr[k]+v*(delta(i,k)*dr[j]+delta(j,k)*dr[i])-C6*dr[i]*dr[j]*dr[k])
                             +2.*v*(n[i]*dr[j]*dr[k]+n[j]*dr[i]*dr[k])
                             +C3*(2.*n[k]*dr[i]*dr[j]+n[j]*delta(i,k)+n[i]*delta(j,k))
                             -C7*n[k]*delta(i,j))); // Rxyy
        
        return F;
    }
    
    private double delta(int i, int j){
        double val=0.; if(i==j)val=1.;
        return val;
    }

    public AbstractMatrix getuLogPart() {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity2DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity2DFS.theFSdata.getR();
        double r=StaticElasticity2DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        if(this.PlaneStress)v=v/(1.+v);
        
        double[] dr=StaticElasticity2DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double c=1./( 8.0*Math.PI*anElasticMat.getLame_M()*(1.-v) );
        double c1=3.-4.*v;
        
        double angle = FundamentalSolution.theFSdata.getMaterial().getAngle();
        double cos=Math.cos(angle); double sin=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, cos); Rot.set(0, 1, -sin);
        Rot.set(1, 0, sin); Rot.set(1, 1, cos);
        
        F.init();
        
        F.set(0, 0, c*(c1));
        F.set(0, 1, 0.);
        
        F.set(1, 0, 0.);
        F.set(1, 1, c*(c1));
//        F=Rot.times(F.times(Rot.transpose()));
        return F;
    }

    public AbstractMatrix getuNonLogPart(double val) {
        ElasticMat anElasticMat=(ElasticMat) StaticElasticity2DFS.theFSdata.getMaterial();
        double[] R=StaticElasticity2DFS.theFSdata.getR();
        double r=StaticElasticity2DFS.theFSdata.getAbsR();
        double v=anElasticMat.getPoissonRatio();
        if(this.PlaneStress)v=v/(1.+v);
        
        double[] dr=StaticElasticity2DFS.theFSdata.getDR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double c=1./( 8.0*Math.PI*anElasticMat.getLame_M()*(1.-v) );
        double c1=3.-4.*v;
        
        double angle = FundamentalSolution.theFSdata.getMaterial().getAngle();
        double cos=Math.cos(angle); double sin=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, cos); Rot.set(0, 1, -sin);
        Rot.set(1, 0, sin); Rot.set(1, 1, cos);
        
        F.init();
        
        F.set(0, 0, c*c1*val+c*(dr[0]*dr[0]));
        F.set(0, 1, c*dr[0]*dr[1]);
        
        F.set(1, 0, c*dr[0]*dr[1]);
        F.set(1, 1, c*c1*val+c*(dr[1]*dr[1]));
        
//        F=Rot.times(F.times(Rot.transpose()));
        
        return F;
    }

    public AbstractMatrix getvLogPart() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public AbstractMatrix getvNonLogPart() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public void setPlaneStress(boolean val){this.PlaneStress=val;}
    
    public boolean getPlaneStress(){return this.PlaneStress;}

}
