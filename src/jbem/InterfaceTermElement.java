/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class InterfaceTermElement{
    protected Domain theDomain;
    protected Element theElement;
    protected int theNode;
    protected int Variable;
    // Variable is 0 for tractions, 1 for displacement, 2 for velocity
    protected int dof_sequence;
    protected double coef;

    // constructor
    public InterfaceTermElement(){}

    public InterfaceTermElement(Domain theDomain,Element theElement, int theNode, int dof_sequence, double coef){
        this.theDomain=theDomain;
        this.theElement=theElement;
        this.theNode=theNode;
        this.Variable=0;
        this.dof_sequence=dof_sequence;
        this.coef=coef;
    }

    public InterfaceTermElement(Domain theDomain,Element theElement, int theNode, int dof_sequence){
        this.theDomain=theDomain;
        this.theElement=theElement;
        this.theNode=theNode;
        this.Variable=0;
        this.dof_sequence=dof_sequence;
        this.coef=1.0;
    }

    public int getDOF(int ndof){
        int dof;
        Node aNode=theElement.getNode(theNode);
        dof = aNode.getpEFTable(theElement, ndof)[dof_sequence-1];
        return dof;
    }

    public double getCoef(){
        return this.coef;
    }
    
    public int getNodeID(){
        return this.theNode;
    }

    public void Print(){
         String SVariable;
         switch(Variable){
             case 0: SVariable = " on traction"; break;
             case 1: SVariable = " on displacement"; break;
             case 2: SVariable = " on velocity"; break;
             default:SVariable = " ! unknown !"; break;
         }
         System.out.println("on Domain: "+this.theDomain.getID()+" on Element: "+theElement.getID()+", on Node: "+this.theNode+SVariable+", dof sequence: "+dof_sequence+", with coefficient = "+this.coef);
     }
}
