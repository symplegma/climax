/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 *
 * @author pchr
 */
public class InterfaceEquation {
    protected int numOfNTerms=0;
    protected int numOfETerms=0;
    protected int id;
    protected double[] Interfacevalue;
    protected List<InterfaceTerm> TermsOfNodeInterface = new Vector() ;
    protected List<InterfaceTermElement> TermsOfElemInterface = new Vector() ;

    public InterfaceEquation(int id, double[] Constraintvalue){
        this.id=id;
        this.Interfacevalue=new double[Constraintvalue.length];
        for(int i=0; i<Constraintvalue.length; i++){
            this.Interfacevalue[i]=Constraintvalue[i];
        }
    }

    public InterfaceEquation(int id){
        this.id=id;
    }

    // methods
    public void setInterfacevalue(double[] Constraintvalue){
        this.Interfacevalue=new double[Constraintvalue.length];
        for(int i=0; i<Constraintvalue.length; i++){
            this.Interfacevalue[i]=Constraintvalue[i];
        }
    }

    public void addNodeInterfaceTerm(InterfaceTerm cterm){
        this.TermsOfNodeInterface.add(numOfNTerms,cterm);
        ++this.numOfNTerms;
    }

    public void addElemInterfaceTerm(InterfaceTermElement cterm){
        this.TermsOfElemInterface.add(numOfETerms,cterm);
        ++this.numOfETerms;
    }

    public double getInterfacevalue(int step){
        double val;
        if(Interfacevalue!=null){
            val=Interfacevalue[step];
        }else{
            val=0.;
        }
        return val;
    }

    public int getID(){
        return this.id;
    }

    public List getNodeInterfaceTerm(){
        return this.TermsOfNodeInterface;
    }

    public List getElementConstraintTerm(){
        return this.TermsOfElemInterface;
    }

    public void Print(){
         System.out.println("Interface equation id : "+this.id+", num of NTerms = "+numOfNTerms+", num of ETerms = "+numOfETerms);
         System.out.println("constant value (length = "+Interfacevalue.length+") of equation Constraintvalue[0] = "+Interfacevalue[0]);
         for(Iterator<InterfaceTerm> nit=this.getNodeInterfaceTerm().iterator() ; nit.hasNext();){
             InterfaceTerm NodeConstraintTerm = nit.next();
             NodeConstraintTerm.Print();
         }
         for(Iterator<InterfaceTermElement> nit=this.getElementConstraintTerm().iterator() ; nit.hasNext();){
             InterfaceTermElement aConstraintTermElement = nit.next();
             aConstraintTermElement.Print();
         }
         System.out.println();
     }

}
