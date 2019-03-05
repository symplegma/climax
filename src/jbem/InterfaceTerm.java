/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class InterfaceTerm{
    protected Domain theDomain;
    protected Node theNode;
    protected int Variable;
    // Variable is 0 for tractions, 1 for displacement, 2 for velocity
    protected int dof_sequence;
    protected double coef;

    public InterfaceTerm(Domain theDomain,Node theNode, int Variable, int dof_sequence, double coef){
        this.theDomain=theDomain;
        this.theNode=theNode;
        this.Variable=Variable;
        this.dof_sequence=dof_sequence;
        this.coef=coef;
    }

    public InterfaceTerm(Domain theDomain,Node theNode, int Variable, int dof_sequence){
        this.theDomain=theDomain;
        this.theNode=theNode;
        this.Variable=Variable;
        this.dof_sequence=dof_sequence;
    }

    public int getDOF(){
        int dof=0;
        switch(Variable){
            case 1: dof = theNode.getuEFTable()[dof_sequence-1];break;
            case 2: dof = theNode.getvEFTable()[dof_sequence-1];break;
            default: System.exit(201); break;
        }
        return dof;
    }

    public int getNodeID(){
        return this.theNode.getID();
    }

    public Node getNode(){
        return this.theNode;
    }

    public int getDomainID(){
        return this.theDomain.getID();
    }

    public int getDOFsequence(){
        return this.dof_sequence;
    }

    public double getCoef(){
        return this.coef;
    }

    public void Print(){
        String SVariable;
        switch(Variable){
            case 0: SVariable = " on traction"; break;
            case 1: SVariable = " on displacement"; break;
            case 2: SVariable = " on velocity"; break;
            default:SVariable = " ! unknown !"; break;
        }
        System.out.println("on Domain: "+this.theDomain.getID()+" on Node: "+this.theNode.getID()+SVariable+", dof sequence: "+dof_sequence+", with coefficient = "+this.coef);
    }

}
