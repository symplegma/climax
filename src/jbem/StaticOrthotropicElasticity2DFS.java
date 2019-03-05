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
public class StaticOrthotropicElasticity2DFS extends FundamentalSolution implements LogFunction{
    private boolean PlaneStress=false;
    
    // constructor
    public StaticOrthotropicElasticity2DFS(){
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
        OrthotropicMaterial2D anElasticMat=(OrthotropicMaterial2D) StaticElasticity3DFS.theFSdata.getMaterial();
        anElasticMat.setPlaneStress(this.PlaneStress);
        double[] betas = anElasticMat.getBetas();
        double[] R=FundamentalSolution.theFSdata.getR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        double angle = ((OrthotropicMaterial2D) FundamentalSolution.theFSdata.getMaterial()).getAngle();
        double c=Math.cos(angle); double s=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, c); Rot.set(0, 1, -s);
        Rot.set(1, 0, s); Rot.set(1, 1, c);
        double y11,y12;
        double y21,y22;
        
        double m = betas[0]*Math.pow(R[0], 4)+2.*(betas[2]+betas[3]/2.)*Math.pow(R[0], 2)*Math.pow(R[1], 2)+
                betas[1]*Math.pow(R[1], 4);
        double m1 = (betas[2]+betas[3]/2.)*Math.pow(R[0], 2)+betas[1]*Math.pow(R[1], 2);
        double m2 = betas[0]*Math.pow(R[0], 2)+(betas[2]+betas[3]/2.)*Math.pow(R[1], 2);
        double l = (Math.sqrt(2.)*R[0]*R[1])/(Math.sqrt(betas[0])*Math.pow(R[0], 2)+Math.sqrt(betas[1])*Math.pow(R[1], 2));
        double l1 = (anElasticMat.getBetaPlus().getReal()*Math.pow(R[0], 2))/(m1);
        double l2 = (anElasticMat.getBetaPlus().getReal()*Math.pow(R[1], 2))/(m2);
        
        F.init();
        
        if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.REAL){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = (Math.pow(anElasticMat.getBetaPlus().getReal(), 2)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+Math.pow(betas[3]/2., 2))/(8.*Math.PI*anElasticMat.getBetaPlus().getReal()*anElasticMat.getBetaMinus().getReal());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)-(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (4.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getReal()*Math.sqrt(betas[1]));
            double C22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)-(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (4.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getReal()*Math.sqrt(betas[0]));
            
            y11=-K11*Math.log(m/betas[1])-C11*Math.atan(anElasticMat.getBetaMinus().getReal()*l1);
            y12=K12*Math.log((1.+anElasticMat.getBetaMinus().getReal()*l)/(1.-anElasticMat.getBetaMinus().getReal()*l));
            y21=K12*Math.log((1.+anElasticMat.getBetaMinus().getReal()*l)/(1.-anElasticMat.getBetaMinus().getReal()*l));
            y22=-K22*Math.log(m/betas[0])-C22*Math.atan(anElasticMat.getBetaMinus().getReal()*l2);

//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.IMAGINARY){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = -(-Math.pow(anElasticMat.getBetaPlus().getReal(), 2)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+Math.pow(betas[3]/2., 2))/(4.*Math.PI*anElasticMat.getBetaPlus().getReal()*anElasticMat.getBetaMinus().getImaginary());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = -(-(Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)-(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getImaginary()*Math.sqrt(betas[1]));
            double C22 = -(-(Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)-(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getImaginary()*Math.sqrt(betas[0]));
            
            y11=-K11*Math.log(m/betas[1])-C11*Math.log((1.-anElasticMat.getBetaMinus().getImaginary()*l1)/(1.+anElasticMat.getBetaMinus().getImaginary()*l1));
            y12=K12*Math.atan(-anElasticMat.getBetaMinus().getImaginary()*l);
            y21=K12*Math.atan(-anElasticMat.getBetaMinus().getImaginary()*l);
            y22=-K22*Math.log(m/betas[0])-C22*Math.log((1.-anElasticMat.getBetaMinus().getImaginary()*l2)/(1.+anElasticMat.getBetaMinus().getImaginary()*l2));

//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
            
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.ZERO){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = (Math.pow(betas[3]/2., 2))/(4.*Math.PI*anElasticMat.getBetaPlus().getReal());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = (Math.pow(betas[3]/2., 2))/
                    (4.*Math.sqrt(2.)*Math.PI*Math.sqrt(betas[1]));
            double C22 = (Math.pow(betas[3]/2., 2))/
                    (4.*Math.sqrt(2.)*Math.PI*Math.sqrt(betas[0]));
            
            y11=-K11*Math.log(m/betas[1])+C11*l1;
            y12=K12*l;
            y21=K12*l;
            y22=-K22*Math.log(m/betas[0])+C22*l2;

//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        }else{
            System.err.println("Problem in get_u_fund of class "+this.getClass());
            System.exit(ndofs);
        }
        return F;
    }

    @Override
    public AbstractMatrix get_p_fund() {
        OrthotropicMaterial2D anElasticMat=(OrthotropicMaterial2D) StaticElasticity3DFS.theFSdata.getMaterial();
        anElasticMat.setPlaneStress(this.PlaneStress);
        double[] betas = anElasticMat.getBetas();
        double[] R=FundamentalSolution.theFSdata.getR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        double[] n=StaticOrthotropicElasticity2DFS.theFSdata.getOutwardNormal();
        double angle = ((OrthotropicMaterial2D) FundamentalSolution.theFSdata.getMaterial()).getAngle();
        double c=Math.cos(angle); double s=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, c); Rot.set(0, 1, -s);
        Rot.set(1, 0, s); Rot.set(1, 1, c);
        double y11,y12;
        double y21,y22;
        double A1 = 0,A2 = 0,B1 = 0,B2 = 0,C1 = 0,C2 = 0,D1 = 0,D2 = 0;
        double m = betas[0]*Math.pow(R[0], 4)+2.*(betas[2]+betas[3]/2.)*Math.pow(R[0], 2)*Math.pow(R[1], 2)+
                betas[1]*Math.pow(R[1], 4);
        double bplus = anElasticMat.getBetaPlus().getReal();
        double bminus2 = 0.;
        if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.REAL){
            bminus2 = anElasticMat.getBetaMinus().getReal()*anElasticMat.getBetaMinus().getReal();
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.IMAGINARY){
            bminus2 = -anElasticMat.getBetaMinus().getImaginary()*anElasticMat.getBetaMinus().getImaginary();
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.ZERO){
            bminus2 = 0.;
        }else{
            System.err.println("Problem in get_u_fund of class "+this.getClass());
            System.exit(ndofs);
        }
        A1=Math.sqrt(betas[1])*(bplus*bplus-betas[3]/2.)/(2.*Math.sqrt(2.)*Math.PI*bplus);
        A2=Math.sqrt(betas[0])*(bplus*bplus-betas[3]/2.)/(2.*Math.sqrt(2.)*Math.PI*bplus);
        B1=Math.sqrt(betas[0])*(bplus*bplus+betas[3]/2.)/(2.*Math.sqrt(2.)*Math.PI*bplus);
        B2=Math.sqrt(betas[1])*(bplus*bplus+betas[3]/2.)/(2.*Math.sqrt(2.)*Math.PI*bplus);
        C1=((bplus*bplus+betas[3]/2.)*Math.sqrt(betas[0]*betas[1])-2.*(bminus2+betas[3]/2.)*bplus*bplus)/(2.*Math.sqrt(2.)*Math.PI*bplus*Math.sqrt(betas[0]));
        C2=((bplus*bplus-betas[3]/2.)*betas[0])/(2.*Math.sqrt(2.)*Math.PI*bplus*Math.sqrt(betas[1]));
        D1=((bplus*bplus-betas[3]/2.)*betas[1])/(2.*Math.sqrt(2.)*Math.PI*bplus*Math.sqrt(betas[0]));
        D2=((bplus*bplus+betas[3]/2.)*Math.sqrt(betas[0]*betas[1])-2.*(bminus2+betas[3]/2.)*bplus*bplus)/(2.*Math.sqrt(2.)*Math.PI*bplus*Math.sqrt(betas[1]));
        
        F.init();
        
        y11=-(A1*Math.pow(R[1], 2)+B1*Math.pow(R[0], 2))*(R[0]*n[0]+R[1]*n[1])/m;
        y12=(C1*Math.pow(R[0], 2)+D1*Math.pow(R[1], 2))*R[1]*n[0]/m-
                 (A2*Math.pow(R[0], 2)+B2*Math.pow(R[1], 2))*R[0]*n[1]/m;
        y21=-(A1*Math.pow(R[1], 2)+B1*Math.pow(R[0], 2))*R[1]*n[0]/m+
                 (C2*Math.pow(R[0], 2)+D2*Math.pow(R[1], 2))*R[0]*n[1]/m;
        y22=-(A2*Math.pow(R[0], 2)+B2*Math.pow(R[1], 2))*(R[0]*n[0]+R[1]*n[1])/m;

//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        
        return F;
    }

    public AbstractMatrix getuLogPart() {
        OrthotropicMaterial2D anElasticMat=(OrthotropicMaterial2D) StaticElasticity3DFS.theFSdata.getMaterial();
        anElasticMat.setPlaneStress(this.PlaneStress);
        double[] betas = anElasticMat.getBetas();
        double[] R=FundamentalSolution.theFSdata.getR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        double angle = ((OrthotropicMaterial2D) FundamentalSolution.theFSdata.getMaterial()).getAngle();
        double c=Math.cos(angle); double s=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, c); Rot.set(0, 1, -s);
        Rot.set(1, 0, s); Rot.set(1, 1, c);
        double y11,y12;
        double y21,y22;
        
        double m = betas[0]*Math.pow(R[0], 4)+2.*(betas[2]+betas[3]/2.)*Math.pow(R[0], 2)*Math.pow(R[1], 2)+
                betas[1]*Math.pow(R[1], 4);
        double m1 = (betas[2]+betas[3]/2.)*Math.pow(R[0], 2)+betas[1]*Math.pow(R[1], 2);
        double m2 = betas[0]*Math.pow(R[0], 2)+(betas[2]+betas[3]/2.)*Math.pow(R[1], 2);
        double l = (Math.sqrt(2.)*R[0]*R[1])/(Math.sqrt(betas[0])*Math.pow(R[0], 2)+Math.sqrt(betas[1])*Math.pow(R[1], 2));
        double l1 = (anElasticMat.getBetaPlus().getReal()*Math.pow(R[0], 2))/(m1);
        double l2 = (anElasticMat.getBetaPlus().getReal()*Math.pow(R[1], 2))/(m2);

        
        F.init();
        
        if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.REAL){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = (Math.pow(anElasticMat.getBetaPlus().getReal(), 2)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+Math.pow(betas[3]/2., 2))/(8.*Math.PI*anElasticMat.getBetaPlus().getReal()*anElasticMat.getBetaMinus().getReal());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)-(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (4.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getReal()*Math.sqrt(betas[1]));
            double C22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)-(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (4.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getReal()*Math.sqrt(betas[0]));
            
            y11=K11;
            y12=0.0;
            y21=0.0;
            y22=K22;
            
//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.IMAGINARY){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = -(-Math.pow(anElasticMat.getBetaPlus().getReal(), 2)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+Math.pow(betas[3]/2., 2))/(4.*Math.PI*anElasticMat.getBetaPlus().getReal()*anElasticMat.getBetaMinus().getImaginary());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = -(-(Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)-(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getImaginary()*Math.sqrt(betas[1]));
            double C22 = -(-(Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)-(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getImaginary()*Math.sqrt(betas[0]));
            
            y11=K11;
            y12=0.0;
            y21=0.0;
            y22=K22;
            
//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
            
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.ZERO){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = (Math.pow(betas[3]/2., 2))/(4.*Math.PI*anElasticMat.getBetaPlus().getReal());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = (Math.pow(betas[3]/2., 2))/
                    (4.*Math.sqrt(2.)*Math.PI*Math.sqrt(betas[1]));
            double C22 = (Math.pow(betas[3]/2., 2))/
                    (4.*Math.sqrt(2.)*Math.PI*Math.sqrt(betas[0]));
            
            y11=K11;
            y12=0.0;
            y21=0.0;
            y22=K22;
            
//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        }else{
            System.err.println("Problem in get_u_fund of class "+this.getClass());
            System.exit(ndofs);
        }
        
        return F;
    }

    public AbstractMatrix getuNonLogPart(double val) {
        OrthotropicMaterial2D anElasticMat=(OrthotropicMaterial2D) FundamentalSolution.theFSdata.getMaterial();
        anElasticMat.setPlaneStress(this.PlaneStress);
        double[] betas = anElasticMat.getBetas();
        double[] R=FundamentalSolution.theFSdata.getR();
        AbstractMatrix F=new AbstractMatrix(this.ndofs,this.ndofs);
        double r=FundamentalSolution.theFSdata.getAbsR();
        double angle = ((OrthotropicMaterial2D) FundamentalSolution.theFSdata.getMaterial()).getAngle();
        double c=Math.cos(angle); double s=Math.sin(angle);
        Rot = new AbstractMatrix(2,2);
        Rot.set(0, 0, c); Rot.set(0, 1, -s);
        Rot.set(1, 0, s); Rot.set(1, 1, c);
        double y11,y12;
        double y21,y22;
        
        double m = betas[0]*Math.pow(R[0], 4)+2.*(betas[2]+betas[3]/2.)*Math.pow(R[0], 2)*Math.pow(R[1], 2)+
                betas[1]*Math.pow(R[1], 4);
        double m1 = (betas[2]+betas[3]/2.)*Math.pow(R[0], 2)+betas[1]*Math.pow(R[1], 2);
        double m2 = betas[0]*Math.pow(R[0], 2)+(betas[2]+betas[3]/2.)*Math.pow(R[1], 2);
        double l = (Math.sqrt(2.)*R[0]*R[1])/(Math.sqrt(betas[0])*Math.pow(R[0], 2)+Math.sqrt(betas[1])*Math.pow(R[1], 2));
        double l1 = (anElasticMat.getBetaPlus().getReal()*Math.pow(R[0], 2))/(m1);
        double l2 = (anElasticMat.getBetaPlus().getReal()*Math.pow(R[1], 2))/(m2);
        
        F.init();
        
        if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.REAL){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = (Math.pow(anElasticMat.getBetaPlus().getReal(), 2)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+Math.pow(betas[3]/2., 2))/(8.*Math.PI*anElasticMat.getBetaPlus().getReal()*anElasticMat.getBetaMinus().getReal());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)-(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (4.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getReal()*Math.sqrt(betas[1]));
            double C22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getReal(), 2)-(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (4.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getReal()*Math.sqrt(betas[0]));
            
            y11=K11*val-K11*Math.log(m/(betas[1]*r))-C11*Math.atan(anElasticMat.getBetaMinus().getReal()*l1);
            y12=K12*Math.log((1.+anElasticMat.getBetaMinus().getReal()*l)/(1.-anElasticMat.getBetaMinus().getReal()*l));
            y21=K12*Math.log((1.+anElasticMat.getBetaMinus().getReal()*l)/(1.-anElasticMat.getBetaMinus().getReal()*l));
            y22=K22*val-K22*Math.log(m/(betas[0]*r))-C22*Math.atan(anElasticMat.getBetaMinus().getReal()*l2);
            
                        
//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.IMAGINARY){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = -(-Math.pow(anElasticMat.getBetaPlus().getReal(), 2)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+Math.pow(betas[3]/2., 2))/(4.*Math.PI*anElasticMat.getBetaPlus().getReal()*anElasticMat.getBetaMinus().getImaginary());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = -(-(Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)-(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getImaginary()*Math.sqrt(betas[1]));
            double C22 = -(-(Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)-(-Math.pow(anElasticMat.getBetaMinus().getImaginary(), 2)+betas[3]/2.)*betas[3]/2.)/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaMinus().getImaginary()*Math.sqrt(betas[0]));
            
            y11=K11*val-K11*Math.log(m/(betas[1]*r))-C11*Math.log((1.-anElasticMat.getBetaMinus().getImaginary()*l1)/(1.+anElasticMat.getBetaMinus().getImaginary()*l1));
            y12=K12*Math.atan(-anElasticMat.getBetaMinus().getImaginary()*l);
            y21=K12*Math.atan(-anElasticMat.getBetaMinus().getImaginary()*l);
            y22=K22*val-K22*Math.log(m/(betas[0]*r))-C22*Math.log((1.-anElasticMat.getBetaMinus().getImaginary()*l2)/(1.+anElasticMat.getBetaMinus().getImaginary()*l2));
                        
//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
            
        }else if(anElasticMat.getMaterialCase()==OrthotropicMaterial2D.OrthotropicMaterialCase.ZERO){
            double K11 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[1]));
            double K12 = (Math.pow(betas[3]/2., 2))/(4.*Math.PI*anElasticMat.getBetaPlus().getReal());
            double K22 = ((Math.pow(anElasticMat.getBetaPlus().getReal(), 2)-betas[3]/2.)*betas[3]/2.+(Math.pow(anElasticMat.getBetaMinus().getReal(), 2)+betas[3]/2.)*Math.pow(anElasticMat.getBetaPlus().getReal(), 2))/
                    (8.*Math.sqrt(2.)*Math.PI*anElasticMat.getBetaPlus().getReal()*Math.sqrt(betas[0]));
            
            double C11 = (Math.pow(betas[3]/2., 2))/
                    (4.*Math.sqrt(2.)*Math.PI*Math.sqrt(betas[1]));
            double C22 = (Math.pow(betas[3]/2., 2))/
                    (4.*Math.sqrt(2.)*Math.PI*Math.sqrt(betas[0]));
            
            y11=K11*val-K11*Math.log(m/(betas[1]*r))+C11*l1;
            y12=K12*l;
            y21=K12*l;
            y22=K22*val-K22*Math.log(m/(betas[0]*r))+C22*l2;
                        
//            F.set(0, 0, c*c*y11-s*c*y21-s*c*y12+s*s*y22);
//            F.set(0, 1, s*c*y11+c*c*y21-s*s*y12-s*c*y22);
//
//            F.set(1, 0, s*c*y11-s*s*y21+c*c*y12-s*c*y22);
//            F.set(1, 1, s*s*y11+s*c*y21+s*c*y12+c*c*y22);
            
            F.set(0, 0, y11);
            F.set(0, 1, y21);

            F.set(1, 0, y12);
            F.set(1, 1, y22);
            F=Rot.times(F.times(Rot.transpose()));
        }else{
            System.err.println("Problem in get_u_fund of class "+this.getClass());
            System.exit(ndofs);
        }
        
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

    @Override
    public AbstractMatrix get_s_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_r_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
