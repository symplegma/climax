/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class MinimizationConstraintTerm {
    protected InterfaceNode theNode;
    protected int Variable;
    // Variable is  0 for normal displ, 1 for tangent, 2 for damage, 
    //              3 for slip_pos, 4 for slip_neg, 
    //              5 for plast_pos, 6 for plast_neg,
    //              7 for normal_pos, 8 normal_neg
    //              9 for tangent_pos, 10 tangent_neg
    protected double coef;
    //private static int numberOfConstraintTerms = 0;
    //private int id;
    private InterfaceElement onElement;

    public MinimizationConstraintTerm(int Variable){
        this.Variable=Variable;
        this.coef=1.0;
    }

    public MinimizationConstraintTerm(InterfaceNode theNode, int Variable, double coef){
        this.theNode=theNode;
        this.Variable=Variable;
        this.coef=coef;
        //++numberOfConstraintTerms;
        //id=numberOfConstraintTerms;
    }

    public MinimizationConstraintTerm(InterfaceNode theNode, int Variable){
        this.theNode=theNode;
        this.Variable=Variable;
        this.coef=1.0;
        //++numberOfConstraintTerms;
        //id=numberOfConstraintTerms;
    }

    public void setOnElement(InterfaceElement onElement){this.onElement=onElement;}

    public int getINodeID(){
        return this.theNode.getID();
    }

    public InterfaceNode getINode(){return this.theNode;}

    public InterfaceElement getIElement(){return this.onElement;}

    public int getDOF(){
        int dof;
        switch(Variable){
            case 0: dof = theNode.getunEFTable(); break;
            case 1: dof = theNode.getutEFTable();break;
            case 2: 
                if(theNode!=null){dof = theNode.getzEFTable(onElement.getID());}else{dof = this.onElement.getINodeHier(1).getzEFTable(onElement.getID());}
                break;
            case 3:
                if(theNode!=null){dof = theNode.gets_posEFTable(onElement.getID());}else{dof = this.onElement.getINodeHier(1).gets_posEFTable(onElement.getID());}
                break;
            case 4:
                if(theNode!=null){dof = theNode.gets_negEFTable(onElement.getID());}else{dof = this.onElement.getINodeHier(1).gets_negEFTable(onElement.getID());}
                break;
            case 5:
                if(theNode!=null){dof = theNode.getp_posEFTable(onElement.getID());}else{dof = this.onElement.getINodeHier(1).getp_posEFTable(onElement.getID());}
                break;
            case 6:
                if(theNode!=null){dof = theNode.getp_negEFTable(onElement.getID());}else{dof = this.onElement.getINodeHier(1).getp_negEFTable(onElement.getID());}
                break;
            case 7: dof = theNode.getun_posEFTable();break;
            case 8: dof = theNode.getun_negEFTable();break;
            case 9: dof = theNode.getut_posEFTable();break;
            case 10: dof = theNode.getut_negEFTable();break;
            default: dof=0; break;
        }
        return dof;
    }

    public void setCoef(double cc){this.coef=cc;}

    public double getCoef(){return this.coef;}

    public int getVariable(){return this.Variable;}

    public void Print() {
        String SVariable;
        switch(Variable){
            case 0: SVariable = " on normal displacement"; break;
            case 1: SVariable = " on tangential displacement"; break;
            case 2: SVariable = " on damage"; break;
            case 3: SVariable = " on slip pos"; break;
            case 4: SVariable = " on slip neg"; break;
            case 5: SVariable = " on plast pos"; break;
            case 6: SVariable = " on plast neg"; break;
            case 7: SVariable = " on normal pos"; break;
            case 8: SVariable = " on normal neg"; break;
            case 9: SVariable = " on tangent pos"; break;
            case 10: SVariable = " on tangent neg"; break;
            default:SVariable = " ! unknown !"; break;
        }
        if(this.onElement!=null){
            if(this.theNode!=null)System.out.println("on InterfaceNode: "+this.theNode.getID()+SVariable+" (dof= "+getDOF()+"), with coefficient = "+this.coef+", on element with id = "+this.onElement.getID());
            if(this.theNode==null)System.out.println("on InterfaceNode: "+0+SVariable+" (dof= "+getDOF()+"), with coefficient = "+this.coef+", on interface element with id = "+this.onElement.getID());
        }else{
            if(this.theNode!=null)System.out.println("on InterfaceNode: "+this.theNode.getID()+SVariable+" (dof= "+getDOF()+"), with coefficient = "+this.coef+", on element with id = "+0);
            if(this.theNode==null)System.out.println("on InterfaceNode: "+0+SVariable+" (dof= "+getDOF()+"), with coefficient = "+this.coef+", on interface element with id = "+0);
        }
    }
}
