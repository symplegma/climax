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
public class TransientElasticity3DFS extends TimeFundamentalSolution{
    private static AbstractMatrix MyMat;
    private double aparam;
    
    // constructor
    public TransientElasticity3DFS(){
        this.ndofs=3;
        RespectiveStatic = new StaticElasticity3DFS();
        MyMat = new AbstractMatrix(3,3);
        aparam=0.6;
        SpaceDimension=3;
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
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        int m=TransientElasticity3DFS.theFSdata.getMstep();
        int n=TransientElasticity3DFS.theFSdata.getNstep();
        //if(m>n){System.out.println("m>n  "+"m= "+m+", n= "+n);}else if(n>m){System.out.println("n>m  "+"m= "+m+", n= "+n);}else{System.out.println("n=m  "+"m= "+m+", n= "+n);}
        double dt=TransientElasticity3DFS.theFSdata.getTimeStep();
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getpTimeDist();
        
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        
        double tcoef=0.;
        double m1,m2;
        double k1,k2=0.;
        
        // aij contribution
        //m1=(m-n)*dt-r/c1;
        //k1=1.0;
        //if(m1<0. || m1>dt){k1=0.0;}
        m1=m*dt-r/c1;
        k1=0.0;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=(n*dt-m1)/(dt); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(m1-(n-1)*dt)/(dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }else{
            //tcoef=1./2.;
            switch(pos){
                case 0: 
                    tcoef=aparam; break;
                case 1:
                    tcoef=1.-aparam; break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Ap().times(tcoef*k1));
        
        // bij contribution
        //m1=(m-n)*dt-r/c2;
        //k1=1.0;
        //if(m1<0. || m1>dt){k1=0.0;}
        m1=m*dt-r/c2;
        k1=0.0;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=(n*dt-m1)/(dt); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(m1-(n-1)*dt)/(dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(2); break;
            }
        }else{
            //tcoef=1./2.;
            switch(pos){
                case 0: 
                    tcoef=aparam; break;
                case 1:
                    tcoef=1.-aparam; break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Bp().times(tcoef*k1));
        
        // cij contribution
        // with c1
        //EDO EIMAI ......... 29/6/2009 
        /*m1=(m-n)*dt-r/c1;
        m2=m1-dt;
        if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c1 ){
            k1=0; k2=0;
            //System.out.println("found in part 1"+" with m= "+m+", n= "+n+", r= "+r+", c1= "+c1+", dt= "+dt);
        }else if( n*dt<m*dt-r/c1 ){
            k1=(n-1)*dt; k2=n*dt;
            //System.out.println("found in part 2"+" with m= "+m+", n= "+n+", r= "+r+", c1= "+c1+", dt= "+dt);
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c1;
            //System.out.println("found in part 3"+" with m= "+m+", n= "+n+", r= "+r+", c1= "+c1+", dt= "+dt);
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=-(k1-k2)/(6.*dt*r*r)*( 2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*(m-n)-3.*dt*(k1+k2)*(1+m-n) ); 
                    tcoef=-(k1-k2)*(2.*k1*k1+2.*k2*k2+6.*dt*dt*m*n-3.*dt*k2*(m+n)+k1*(2.*k2-3.*dt*(m+n)))/(6.*dt*r*r);
                    break;
                case 1:
                    //tcoef=1./(6.*dt*r*r)*( 2.*k1*k1*k1+k2*k2*(-2.*k2+3.*dt*(m-n))+3.*dt*k1*k1*(n-m) ); break;
                    tcoef=(k1-k2)*(2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*m*(n-1)-3.*dt*(k1+k2)*(m+n-1))/(6.*dt*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(3); break;
            }
        }else{
            //tcoef=(k1-k2)*(k1+k2-2.*dt*m)/(4.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Cp().times(tcoef));
        // with c2
        /*m1=(m-n)*dt-r/c2;
        m2=m1-dt;
        if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c2 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c2 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c2;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=-(k1-k2)/(6.*dt*r*r)*( 2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*(m-n)-3.*dt*(k1+k2)*(1+m-n) );
                    tcoef=-(k1-k2)*(2.*k1*k1+2.*k2*k2+6.*dt*dt*m*n-3.*dt*k2*(m+n)+k1*(2.*k2-3.*dt*(m+n)))/(6.*dt*r*r);
                    break;
                case 1:
                    //tcoef=1./(6.*dt*r*r)*( 2.*k1*k1*k1+k2*k2*(-2.*k2+3.*dt*(m-n))+3.*dt*k1*k1*(n-m) ); break;
                    tcoef=(k1-k2)*(2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*m*(n-1)-3.*dt*(k1+k2)*(m+n-1))/(6.*dt*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(4); break;
            }
        }else{
            //tcoef=(k1-k2)*(k1+k2-2.*dt*m)/(4.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        tcoef=-tcoef;
        F=F.plus(this.Cp().times(tcoef));
        
        
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        int m=TransientElasticity3DFS.theFSdata.getMstep();
        int n=TransientElasticity3DFS.theFSdata.getNstep();
        double dt=TransientElasticity3DFS.theFSdata.getTimeStep();
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getuTimeDist();
        //int dist=TransientRodFS.theFSdata.getvTimeDist();
        
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        
        double tcoef=0.;
        double m1,m2;
        double k1,k2=0.;
        
        // dij contribution
        /*m1=(m-n)*dt-r/c1;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c1;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=(n*dt-m1)/dt; break;
                case 1:
                    tcoef=(m1-(n-1)*dt)/(dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(5); break;
            }
        }else{
            //tcoef=1./2.;
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam; break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Dp().times(tcoef*k1));
        // eij contribution
        /*m1=(m-n)*dt-r/c2;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c2;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=(n*dt-m1)/dt; break;
                case 1:
                    tcoef=(m1-(n-1)*dt)/(dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(6); break;
            }
        }else{
            //tcoef=1./2.;
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam; break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Ep().times(tcoef*k1));
        // fij contribution
        // with c1
        m1=(m-n)*dt-r/c1;
        m2=m1-dt;
        /*if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c1 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c1 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c1;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=-(k1-k2)/(6.*dt*r*r)*( 2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*(m-n)-3.*dt*(k1+k2)*(1+m-n) ); 
                    tcoef=-(k1-k2)*(2.*k1*k1+2.*k2*k2+6.*dt*dt*m*n-3.*dt*k2*(m+n)+k1*(2.*k2-3.*dt*(m+n)))/(6.*dt*r*r);
                    break;
                case 1:
                    //tcoef=1./(6.*dt*r*r)*( 2.*k1*k1*k1+k2*k2*(-2.*k2+3.*dt*(m-n))+3.*dt*k1*k1*(n-m) ); break;
                    tcoef=(k1-k2)*(2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*m*(n-1)-3.*dt*(k1+k2)*(m+n-1))/(6.*dt*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(7); break;
            }
        }else{
            //tcoef=(k1-k2)*(k1+k2-2.*dt*m)/(4.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Fp().times(tcoef));
        // with c2
        m1=(m-n)*dt-r/c2;
        m2=m1-dt;
        /*if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c2 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c2 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c2;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=-(k1-k2)/(6.*dt*r*r)*( 2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*(m-n)-3.*dt*(k1+k2)*(1+m-n) ); 
                    tcoef=-(k1-k2)*(2.*k1*k1+2.*k2*k2+6.*dt*dt*m*n-3.*dt*k2*(m+n)+k1*(2.*k2-3.*dt*(m+n)))/(6.*dt*r*r);
                    break;
                case 1:
                    //tcoef=1./(6.*dt*r*r)*( 2.*k1*k1*k1+k2*k2*(-2.*k2+3.*dt*(m-n))+3.*dt*k1*k1*(n-m) ); break;
                    tcoef=(k1-k2)*(2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*m*(n-1)-3.*dt*(k1+k2)*(m+n-1))/(6.*dt*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(8); break;
            }
        }else{
            //tcoef=(k1-k2)*(k1+k2-2.*dt*m)/(4.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        tcoef=-tcoef;
        F=F.plus(this.Fp().times(tcoef));
        // gij contribution
        /*m1=(m-n)*dt-r/c1;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c1;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=-1./dt; break;
                case 1:
                    tcoef=1./dt; break;
                default: System.out.println("error in get_u_fund"); System.exit(9); break;
            }
        }else{
            tcoef=0.;
        }
        F=F.plus(this.Gp().times(tcoef*k1));
        
        // hij contribution
        /*m1=(m-n)*dt-r/c2;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c2;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=-1./dt; break;
                case 1:
                    tcoef=1./dt; break;
                default: System.out.println("error in get_u_fund"); System.exit(10); break;
            }
        }else{
            tcoef=0.;
        }
        F=F.plus(this.Hp().times(tcoef*k1));
        
        
        return F;
    }
    
    public AbstractMatrix get_v_fund() {
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        int m=TransientElasticity3DFS.theFSdata.getMstep();
        int n=TransientElasticity3DFS.theFSdata.getNstep();
        double dt=TransientElasticity3DFS.theFSdata.getTimeStep();
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getpTimeDist();
        
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        
        double tcoef=0.;
        double m1,m2;
        double k1,k2=0.;
        
        // aij contribution
        /*m1=(m-n)*dt-r/c1;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c1;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=-1./dt; break;
                case 1:
                    tcoef=1./dt; break;
                default: System.out.println("error in get_u_fund"); System.exit(11); break;
            }
        }else{
            tcoef=0.;
        }
        F=F.plus(this.Ap().times(tcoef*k1));
        
        // bij contribution
        /*m1=(m-n)*dt-r/c2;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c2;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=-1./dt; break;
                case 1:
                    tcoef=1./dt; break;
                default: System.out.println("error in get_u_fund"); System.exit(12); break;
            }
        }else{
            tcoef=0.;
        }
        F=F.plus(this.Bp().times(tcoef*k1));
        
        // cij contribution
        // dirac part, c1
        /*m1=(m-n)*dt-r/c1;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c1;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/(dt*r*c1); break;
                    tcoef=(n*dt-m1)/dt*(m*dt-m1)/(r*r); break;
                case 1:
                    //tcoef=m1/(dt*r*c1); break;
                    tcoef=(m1-(n-1)*dt)/dt*(m*dt-m1)/(r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(13); break;
            }
        }else{
            //tcoef=1./(2.*r*c1);
            //tcoef=1./2.*(m*dt-m1)/(r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(m*dt-m1)/(r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(m*dt-m1)/(r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Cp().times(tcoef*k1));
        // dirac part, c2
        /*m1=(m-n)*dt-r/c2;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c2;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/(dt*r*c1); break;
                    tcoef=(n*dt-m1)/dt*(m*dt-m1)/(r*r); break;
                case 1:
                    //tcoef=m1/(dt*r*c1); break;
                    tcoef=(m1-(n-1)*dt)/dt*(m*dt-m1)/(r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(14); break;
            }
        }else{
            //tcoef=1./(2.*r*c1);
            //tcoef=1./2.*(m*dt-m1)/(r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(m*dt-m1)/(r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(m*dt-m1)/(r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        tcoef=-tcoef;
        F=F.plus(this.Cp().times(tcoef*k1));
        // Heaviside part, c1
        m1=(m-n)*dt-r/c1;
        m2=m1-dt;
        /*if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c1 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c1 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c1;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=(k1-k2)*(k1+k2-2.*dt)/(r*r*2.*dt);
                    tcoef=(k1-k2)*(k1+k2-2.*dt*n)/(r*r*2.*dt);
                    break;
                case 1:
                    //tcoef=(k2*k2-k1*k1)/(r*r*2.*dt); break;
                    tcoef=-(k1-k2)*(k1+k2-2.*dt*(n-1))/(r*r*2.*dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(15); break;
            }
        }else{
            //tcoef=(k2-k1)/(2.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k2-k1)/(r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k2-k1)/(r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Cp().times(tcoef));
        // Heaviside part, c2
        /*m1=(m-n)*dt-r/c2;
        m2=m1-dt;
        if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c2 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c2 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c2;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=(k1-k2)*(k1+k2-2.*dt)/(r*r*2.*dt);
                    tcoef=(k1-k2)*(k1+k2-2.*dt*n)/(r*r*2.*dt);
                    break;
                case 1:
                    //tcoef=(k2*k2-k1*k1)/(r*r*2.*dt); break;
                    tcoef=-(k1-k2)*(k1+k2-2.*dt*(n-1))/(r*r*2.*dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(16); break;
            }
        }else{
            //tcoef=(k2-k1)/(2.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k2-k1)/(r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k2-k1)/(r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        tcoef=-tcoef;
        F=F.plus(this.Cp().times(tcoef));
        
        return F;
    }
    
    @Override
    public AbstractMatrix get_pdif_fund() {
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        F.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        int m=TransientElasticity3DFS.theFSdata.getMstep();
        int n=TransientElasticity3DFS.theFSdata.getNstep();
        double dt=TransientElasticity3DFS.theFSdata.getTimeStep();
        
        int pos =TransientRodFS.theFSdata.getTimePos();
        int dist=TransientRodFS.theFSdata.getuTimeDist();
        //int dist=TransientRodFS.theFSdata.getvTimeDist();
        
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double r=TransientElasticity3DFS.theFSdata.getAbsR();

        double tcoef=0.;
        double m1,m2;
        double k1,k2=0.;
        
        // dij contribution
        /*m1=(m-n)*dt-r/c1;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c1;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=(n*dt-m1)/dt; break;
                case 1:
                    tcoef=(m1-(n-1)*dt)/(dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(17); break;
            }
        }else{
            //tcoef=1./2.;
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam; break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Dp().times(tcoef*k1));
        // eij contribution
        /*m1=(m-n)*dt-r/c2;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c2;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=(n*dt-m1)/dt; break;
                case 1:
                    tcoef=(m1-(n-1)*dt)/(dt); break;
                default: System.out.println("error in get_u_fund"); System.exit(18); break;
            }
        }else{
            //tcoef=1./2.;
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam; break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Ep().times(tcoef*k1));
        // fij contribution
        // with c1
        m1=(m-n)*dt-r/c1;
        m2=m1-dt;
        /*if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c1 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c1 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c1;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=-(k1-k2)/(6.*dt*r*r)*( 2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*(m-n)-3.*dt*(k1+k2)*(1+m-n) ); 
                    tcoef=-(k1-k2)*(2.*k1*k1+2.*k2*k2+6.*dt*dt*m*n-3.*dt*k2*(m+n)+k1*(2.*k2-3.*dt*(m+n)))/(6.*dt*r*r);
                    break;
                case 1:
                    //tcoef=1./(6.*dt*r*r)*( 2.*k1*k1*k1+k2*k2*(-2.*k2+3.*dt*(m-n))+3.*dt*k1*k1*(n-m) ); break;
                    tcoef=(k1-k2)*(2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*m*(n-1)-3.*dt*(k1+k2)*(m+n-1))/(6.*dt*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(19); break;
            }
        }else{
            //tcoef=(k1-k2)*(k1+k2-2.*dt*m)/(4.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        F=F.plus(this.Fp().times(tcoef));
        // with c2
        m1=(m-n)*dt-r/c2;
        m2=m1-dt;
        /*if( (m1>=0) && (m2>=0) ){
            k1=0.; k2=dt;
        }else if( (m1>=0) && (m2<0) ){
            k1=0.; k2=m1;
        }else if( (m1<0) && (m2>=0) ){
            k1=m2; k2=dt;
        }else if( (m1<0) && (m2<0) ){
            k1=0.; k2=0.;
        }*/
        k1=0.; k2=0.;
        if( (n-1)*dt>m*dt-r/c2 ){
            k1=0; k2=0;
        }else if( n*dt<m*dt-r/c2 ){
            k1=(n-1)*dt; k2=n*dt;
        }else {
            k1=(n-1)*dt; k2=m*dt-r/c2;
        }
        
        if(dist==1){
            switch(pos){
                case 0: 
                    //tcoef=-(k1-k2)/(6.*dt*r*r)*( 2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*(m-n)-3.*dt*(k1+k2)*(1+m-n) ); 
                    tcoef=-(k1-k2)*(2.*k1*k1+2.*k2*k2+6.*dt*dt*m*n-3.*dt*k2*(m+n)+k1*(2.*k2-3.*dt*(m+n)))/(6.*dt*r*r);
                    break;
                case 1:
                    //tcoef=1./(6.*dt*r*r)*( 2.*k1*k1*k1+k2*k2*(-2.*k2+3.*dt*(m-n))+3.*dt*k1*k1*(n-m) ); break;
                    tcoef=(k1-k2)*(2.*(k1*k1+k1*k2+k2*k2)+6.*dt*dt*m*(n-1)-3.*dt*(k1+k2)*(m+n-1))/(6.*dt*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(20); break;
            }
        }else{
            //tcoef=(k1-k2)*(k1+k2-2.*dt*m)/(4.*r*r);
            switch(pos){
                case 0: 
                    //tcoef=(dt-m1)/dt; break;
                    tcoef=aparam*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                case 1:
                    //tcoef=m1/dt; break;
                    tcoef=(1.-aparam)*(k1-k2)*(k1+k2-2.*dt*m)/(2.*r*r); break;
                default: System.out.println("error in get_u_fund"); System.exit(1); break;
            }
        }
        tcoef=-tcoef;
        F=F.plus(this.Fp().times(tcoef));
        // gij contribution
        /*m1=(m-n)*dt-r/c1;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c1;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=-1./dt; break;
                case 1:
                    tcoef=1./dt; break;
                default: System.out.println("error in get_u_fund"); System.exit(21); break;
            }
        }else{
            tcoef=0.;
        }
        F=F.plus(this.Gp().times(tcoef*k1));
        
        // hij contribution
        /*m1=(m-n)*dt-r/c2;
        k1=1.0;
        if(m1<0. || m1>dt){k1=0.0;}*/
        m1=m*dt-r/c2;
        k1=0.; k2=0.;
        if((m1>=(n-1)*dt)&&(m1<=n*dt)){k1=1.0;}
        
        if(dist==1){
            switch(pos){
                case 0: 
                    tcoef=-1./dt; break;
                case 1:
                    tcoef=1./dt; break;
                default: System.out.println("error in get_u_fund"); System.exit(22); break;
            }
        }else{
            tcoef=0.;
        }
        F=F.plus(this.Hp().times(tcoef*k1));
        
        F=F.minus(this.RespectiveStatic.get_p_fund());
        return F;
    }
    
    private AbstractMatrix Ap(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double p=anElasticMat.getDensity();
        
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 1./(4.*pi*p)*1./(c1*c1)*dr[i]*dr[j]/r);
            }
        }
        return MyMat;
    }
    
    private AbstractMatrix Bp(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double p=anElasticMat.getDensity();
        
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(i, j, 1./(4.*pi*p)*1./(c2*c2)*(delta(i,j)/r-dr[i]*dr[j]/r));
            }
        }
        return MyMat;
    }
    
    private AbstractMatrix Cp(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double p=anElasticMat.getDensity();
        
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 1./(4.*pi*p)*(3.*dr[i]*dr[j]/r-delta(i,j)/r));
            }
        }
        return MyMat;
    }
    
    private AbstractMatrix Dp(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double[] n=TransientElasticity3DFS.theFSdata.getOutwardNormal();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        
        //Manolis
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 
                        1./(4*pi)*(2.*c2*c2/(c1*c1)-1.)*dr[j]*n[i]/(r*r)-
                        1./(4*pi)*c2*c2/(c1*c1)*(12.*dr[i]*dr[j]*costh/(r*r)-2.*(dr[j]*n[i]+dr[i]*n[j]+delta(i,j)*costh)/(r*r))
                        );
            }
        }
        
        /* // dominguez
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(i, j, 
                        1./(4.*pi)*2./(r*r)*(c2/c1)*(c2/c1)*(costh*delta(i,j)+dr[i]*n[j])-
                        1./(4.*pi)*12./(r*r)*(c2/c1)*(c2/c1)*dr[i]*dr[j]*costh+
                        1./(4.*pi)*2./(r*r)*(c2/c1)*(c2/c1)*dr[i]*n[j]
                        );
            }
        }
        */
        return MyMat;
    }
    
    private AbstractMatrix Ep(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double[] n=TransientElasticity3DFS.theFSdata.getOutwardNormal();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        
        //Manolis
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 
                        1./(4*pi)*12.*dr[i]*dr[j]*costh/(r*r)-
                        1./(4*pi)*2.*dr[j]*n[i]/(r*r)-
                        1./(4*pi)*3.*(dr[i]*n[j]+delta(i,j)*costh)/(r*r)
                        );
            }
        }
        /* // dominguez
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(i, j, 
                        -1./(4.*pi)*2./(r*r)*(costh*delta(i,j)+dr[j]*n[i])-
                        1./(4.*pi)*1./(r*r)*(costh*delta(i,j)+dr[j]*n[i])+
                        1./(4.*pi)*12./(r*r)*dr[i]*dr[j]*costh-
                        1./(4.*pi)*2./(r*r)*dr[i]*n[j]
                        );
            }
        }
        */
        return MyMat;
    }
    
    private AbstractMatrix Fp(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double[] n=TransientElasticity3DFS.theFSdata.getOutwardNormal();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        //manolis
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 
                        -1./(4*pi)*6.*c2*c2*5.*dr[i]*dr[j]*costh/(r*r)+
                        1./(4*pi)*6.*c2*c2*(dr[j]*n[i]+dr[i]*n[j]+delta(i,j)*costh)/(r*r)
                        );
            }
        }
        
        /* // dominguez
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(i, j, 
                        1./(4.*pi)*6.*c2*c2/(r*r)*(costh*delta(i,j)+dr[j]*n[i])-
                        1./(4.*pi)*30.*c2*c2/(r*r)*dr[i]*dr[j]*costh+
                        1./(4.*pi)*6.*c2*c2/(r*r)*dr[i]*n[j]
                        );
            }
        }
        */
        return MyMat;
    }
    
    private AbstractMatrix Gp(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double[] n=TransientElasticity3DFS.theFSdata.getOutwardNormal();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        //manolis
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 
                        -1./(4*pi)*2.*c2*c2/(c1*c1*c1)*(dr[i]*dr[j]*costh)/r+
                        1./(4*pi)*1./c1*(2.*c2*c2/(c1*c1)-1.)*dr[j]*n[i]/r
                        );
            }
        }
        
        /* // dominguez
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(i, j, 
                        -1./(4.*pi)*2./(r*c2)*c2*c2*c2/(c1*c1*c1)*dr[i]*dr[j]*costh-
                        1./(4.*pi)*1./(r)*(1.-2.*c2*c2/(c1*c1))*1./c1*dr[i]*n[j]
                        );
            }
        }
        */
        return MyMat;
    }
    
    private AbstractMatrix Hp(){
        MyMat.init();
        ElasticMat anElasticMat=(ElasticMat) TransientElasticity3DFS.theFSdata.getMaterial();
        double[] R=TransientElasticity3DFS.theFSdata.getR();
        double r=TransientElasticity3DFS.theFSdata.getAbsR();
        double[] dr=TransientElasticity3DFS.theFSdata.getDR();
        double[] n=TransientElasticity3DFS.theFSdata.getOutwardNormal();
        double c1=(anElasticMat.getLame_L()+2.*anElasticMat.getLame_M())/(anElasticMat.getDensity());
        c1=Math.sqrt(c1);
        double c2=anElasticMat.getLame_M()/anElasticMat.getDensity();
        c2=Math.sqrt(c2);
        double pi=Math.PI;
        double costh=0.0;
        for(int i=0; i<dr.length; i++){
            costh+=dr[i]*n[i];
        }
        //manolis
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(j, i, 
                        1./(4*pi)*2./c2*dr[i]*dr[j]*costh/r-
                        1./(4*pi)*1./c2*(dr[i]*n[j]+delta(i,j)*costh)/r
                        );
            }
        }
        
        /* // dominguez
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                MyMat.set(i, j, 
                        1./(4.*pi)*1./(r*c2)*(costh*delta(i,j)+dr[j]*n[i])-
                        1./(4.*pi)*2./(r*c2)*c2*c2*c2/(c1*c1*c1)*dr[i]*dr[j]*costh
                        );
            }
        }
        */
        return MyMat;
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
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_r_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
