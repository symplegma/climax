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
public class TransientRodFS extends TimeFundamentalSolution{
    
    // constructor
    public TransientRodFS(){
        this.ndofs=1;
        RespectiveStatic = new StaticRodFS();
        SpaceDimension=1;
    }

    @Override
    public int get_u_DOFs() {
        return this.ndofs;
    }

    @Override
    public int get_v_DOFs() {
        return this.ndofs;
    }

    @Override
    public int get_p_DOFs() {
        return this.ndofs;
    }

    @Override
    public AbstractMatrix get_u_fund() {
        ElasticMat anElasticMat=(ElasticMat) FundamentalSolution.theFSdata.getMaterial();
        double[] R=TransientRodFS.theFSdata.getR();
        double r=0.0;
        for (int i=0; i<R.length; i++){
            r+=R[i]*R[i];
        }
        r=Math.sqrt(r);
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double E=anElasticMat.getElasticModulus();
        double p=anElasticMat.getDensity();
        double vel=Math.sqrt(E/p);
        double A=anElasticMat.getArea();
        
        int m=TransientRodFS.theFSdata.getMstep();
        int n=TransientRodFS.theFSdata.getNstep();
        //if(m>n){System.out.println("m>n  "+"m= "+m+", n= "+n);}else if(n>m){System.out.println("n>m  "+"m= "+m+", n= "+n);}else{System.out.println("n=m  "+"m= "+m+", n= "+n);}
        double dt=TransientRodFS.theFSdata.getTimeStep();
        
        double m1=(m-n)*dt-r/vel+dt;
        double m2=(m-n)*dt-dt-r/vel;
        
        double k1=0.; double k2=0.;
        
        if( (n-1)*dt>m*dt-r/vel ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/vel ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/vel;
        }
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getpTimeDist();
        
        F.init();
        
        if(dist==1){
            switch(pos){
                case 0: F.set(0, 0, (n*(k1-k2)-(k1*k1-k2*k2)/(2.*dt))*vel/(2.*E*A) ); break;
                case 1: F.set(0, 0, (n*(k2-k1)+k1-k2+(k1*k1-k2*k2)/(2.*dt))*vel/(2.*E*A) ); break;
            }
        }else{
            F.set(0, 0, -(-k1+k2)*vel/(4.*E*A) );
        }
        k1=1.;
        if(r!=0.)k1=-1.;
        return F.times(k1*E*A);
    }

    @Override
    public AbstractMatrix get_p_fund() {
        ElasticMat anElasticMat=(ElasticMat) FundamentalSolution.theFSdata.getMaterial();
        double[] R=TransientRodFS.theFSdata.getR();
        double[] nv=TransientRodFS.theFSdata.getOutwardNormal();
        double r=0.0;
        for (int i=0; i<R.length; i++){
            r+=R[i]*R[i];
        }
        r=Math.sqrt(r);
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double E=anElasticMat.getElasticModulus();
        double p=anElasticMat.getDensity();
        double vel=Math.sqrt(E/p);
        double A=anElasticMat.getArea();
        
        int m=TransientRodFS.theFSdata.getMstep();
        int n=TransientRodFS.theFSdata.getNstep();
        double dt=TransientRodFS.theFSdata.getTimeStep();
        
        double m1=m*dt-r/vel;
        double k1=0.0;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getuTimeDist();
        
        F.init();
        
        if(dist==1){
            switch(pos){
                case 0: 
                    if(R[0]!=0.){F.set(0, 0, r/(R[0])*(n*dt-m1)/(dt*2.*E*A));}
                    else{F.set(0, 0, -nv[0]*(n*dt-m1)/(dt*2.*E*A));}
                break;
                case 1: 
                    if(R[0]!=0.){F.set(0, 0, r/(R[0])*(m1-(n-1)*dt)/(dt*2.*E*A));}
                    else{F.set(0, 0, -nv[0]*(m1-(n-1)*dt)/(dt*2.*E*A));}
                break;
            }
        }else{
            if(R[0]!=0.){F.set(0, 0,  r/(R[0])/(4.*E*A));}
            else{F.set(0, 0, -nv[0]/(4.*E*A));}
        }
        
        return F.times(-k1*E*A);
    }

    @Override
    public AbstractMatrix get_v_fund() {
        ElasticMat anElasticMat=(ElasticMat) FundamentalSolution.theFSdata.getMaterial();
        double[] R=TransientRodFS.theFSdata.getR();
        double[] nv=TransientRodFS.theFSdata.getOutwardNormal();
        double r=0.0;
        for (int i=0; i<R.length; i++){
            r+=R[i]*R[i];
        }
        r=Math.sqrt(r);
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double E=anElasticMat.getElasticModulus();
        double p=anElasticMat.getDensity();
        double vel=Math.sqrt(E/p);
        double A=anElasticMat.getArea();
        
        int m=TransientRodFS.theFSdata.getMstep();
        int n=TransientRodFS.theFSdata.getNstep();
        double dt=TransientRodFS.theFSdata.getTimeStep();
        
        double m1=m*dt-r/vel;
        double k1=0.0;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getpTimeDist();
        
        F.init();
        
        if(dist==1){
            switch(pos){
                case 0: 
                    F.set(0, 0, -vel*(n*dt-m1)/(dt*2.*E*A));
                break;
                case 1: 
                    F.set(0, 0, -vel*(m1-(n-1)*dt)/(dt*2.*E*A));
                break;
            }
        }else{
            F.set(0, 0, -vel/(4.*E*A));
        }
        if(r!=0.)k1*=-1.;
        return F.times(k1*E*A);
    }

    @Override
    public AbstractMatrix get_pdif_fund() {
        ElasticMat anElasticMat=(ElasticMat) FundamentalSolution.theFSdata.getMaterial();
        double[] R=FundamentalSolution.theFSdata.getR();
        double[] nv=StaticElasticity3DFS.theFSdata.getOutwardNormal();
        double r=0.0;
        for (int i=0; i<R.length; i++){
            r+=R[i]*R[i];
        }
        r=Math.sqrt(r);
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        
        double E=anElasticMat.getElasticModulus();
        double p=anElasticMat.getDensity();
        double vel=Math.sqrt(E/p);
        double A=anElasticMat.getArea();
        
        int m=TransientRodFS.theFSdata.getMstep();
        int n=TransientRodFS.theFSdata.getNstep();
        double dt=TransientRodFS.theFSdata.getTimeStep();
        
        double m1=m*dt-r/vel;
        double k1=0.0;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getuTimeDist();
        
        F.init();
        
        if(dist==1){
            switch(pos){
                case 0: 
                    F.set(0, 0, r/(R[0])*(n*dt-m1)/(dt*2.*E*A));
                break;
                case 1: 
                    F.set(0, 0, r/(R[0])*(m1-(n-1)*dt)/(dt*2.*E*A));
                break;
            }
        }else{
            F.set(0, 0,  r/(R[0])/(4.*E*A));
        }
        F=F.times(-k1);
        
        F=F.minus(this.RespectiveStatic.get_p_fund());

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
