/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gendomain;

import java.util.List;
import java.util.Vector;
/**
 *
 * @author pchr
 */
public class ConstraintEquation {
    protected int numOfNTerms=0;
    protected int numOfETerms=0;
    protected int id;
    protected double[] Constraintvalue;
    protected List<ConstraintTerm> TermsOfConstraint = new Vector() ;
    
    // constructor
    public ConstraintEquation(){}

    public ConstraintEquation(int id){
        this.id=id;
    }
    
    public ConstraintEquation(int id, double[] Constraintvalue){
        this.id=id;
        this.Constraintvalue=new double[Constraintvalue.length];
        System.arraycopy(Constraintvalue, 0, this.Constraintvalue, 0, Constraintvalue.length);
    }
    
    
    public ConstraintEquation(int id, double Constraintvalue){
        this.id=id;
        this.Constraintvalue=new double[1];
        this.Constraintvalue[0]=Constraintvalue;
    }
    
    public void addConstraintTerm(ConstraintTerm cterm){
        this.TermsOfConstraint.add(numOfNTerms,cterm);
        ++this.numOfNTerms;
    }
    
    public double getConstraintvalue(int step){
        double val;
        if(Constraintvalue!=null){
            val=Constraintvalue[step];
        }else{
            val=0.;
        }
        return val;
    }
    
    public double getConstraintvalue(){
        double val;
        if(Constraintvalue!=null){
            val=Constraintvalue[0];
        }else{
            val=0.;
        }
        return val;
    }
    
    public int getID(){
        return this.id;
    }
}
