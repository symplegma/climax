/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gendomain;

/**
 *
 * @author pchr
 */
public class ConstraintTerm {
    protected Node theNode;
    protected int Variable; 
    // remained from jBEM: Variable is 0 for tractions, 1 for displacement, 2 for velocity
    protected int dof_sequence;
    protected double coef;
//    private static int numberOfConstraintTerms = 0;
    //private int id;
    
    // constructor
    public ConstraintTerm(){}
    
    public ConstraintTerm(Node theNode, int dof_sequence, double coef, int Variable){
        if(theNode==null)throw new NullPointerException();
        this.theNode=theNode;
        this.Variable=Variable;
        this.dof_sequence=dof_sequence;
        this.coef=coef;
//        ++numberOfConstraintTerms;
//        id=numberOfConstraintTerms;
    }
    
    public ConstraintTerm(Node theNode, int dof_sequence, int Variable){
        if(theNode==null)throw new NullPointerException();
        this.theNode=theNode;
        this.Variable=Variable;
        this.dof_sequence=dof_sequence;
        this.coef=1.0;
//        ++numberOfConstraintTerms;
//        id=numberOfConstraintTerms;
    }
    
    public ConstraintTerm(Node theNode, int dof_sequence){
        if(theNode==null)throw new NullPointerException();
        this.theNode=theNode;
        this.Variable=0;
        this.dof_sequence=dof_sequence;
        this.coef=1.0;
//        ++numberOfConstraintTerms;
//        id=numberOfConstraintTerms;
    }
    
    public ConstraintTerm(Node theNode){
        if(theNode==null)throw new NullPointerException();
        this.theNode=theNode;
        this.Variable=0;
        this.dof_sequence=1;
        this.coef=1.0;
//        ++numberOfConstraintTerms;
//        id=numberOfConstraintTerms;
    }

//    public int getID(){
//        return this.id;
//    }

    public int getNodeID(){
        return this.theNode.getID();
    }
}
