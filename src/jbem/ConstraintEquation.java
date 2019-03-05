/*******************************************************************************
* Climax.                                                                      *
* Copyright (C) 2009-2017 C.G. Panagiotopoulos [http://www.symplegma.org]      *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): C.G. Panagiotopoulos (pchr76@gmail.com)
// *****************************************************************************

package jbem;

import java.util.Iterator;
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
    protected List<ConstraintTerm> TermsOfNodeConstraint = new Vector() ;
    protected List<ConstraintTermElement> TermsOfElemConstraint = new Vector() ;
    private static int numofConstraintEquations=0;
    
    // constructor
    public ConstraintEquation(){}

    public ConstraintEquation(int id){
        this.id=id;
        numofConstraintEquations++;
    }
    
    public ConstraintEquation(int id, double[] Constraintvalue){
        this.id=id;
        numofConstraintEquations++;
        this.Constraintvalue=new double[Constraintvalue.length];
        System.arraycopy(Constraintvalue, 0, this.Constraintvalue, 0, Constraintvalue.length);
    }
    
    
    // methods
    public void addNodeConstraintTerm(ConstraintTerm cterm){
        this.TermsOfNodeConstraint.add(numOfNTerms,cterm);
        ++this.numOfNTerms;
    }
    
    public void addElemConstraintTerm(ConstraintTermElement cterm){
        this.TermsOfElemConstraint.add(numOfETerms,cterm);
        ++this.numOfETerms;
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

    public void setCVFORMATRIXENERGYONLY(double val){
        for(int i=0;i<Constraintvalue.length;i++){
            Constraintvalue[i]=val;
        }
    }
    
    public int getID(){
        return this.id;
    }
    
    public List getNodeConstraintTerm(){
        return this.TermsOfNodeConstraint;
    }
    
    public List getElementConstraintTerm(){
        return this.TermsOfElemConstraint;
    }
    
    public void equalp(Element elem1, Element elem2, int node_id, int dof){
        ConstraintTermElement aConstraintTermElement = new 
                        ConstraintTermElement(elem1,node_id,dof, 1.0);
        this.addElemConstraintTerm(aConstraintTermElement);
        aConstraintTermElement = new 
                        ConstraintTermElement(elem2,node_id,dof, -1.0);
        this.addElemConstraintTerm(aConstraintTermElement);
            
    }
    
    public void relatep(Element elem1, Element elem2, int node_id1, int node_id2, int dof, double coef1, double coef2){
        ConstraintTermElement aConstraintTermElement = new 
                        ConstraintTermElement(elem1,node_id1,dof, coef1);
        this.addElemConstraintTerm(aConstraintTermElement);
        aConstraintTermElement = new 
                        ConstraintTermElement(elem2,node_id2,dof, coef2);
        this.addElemConstraintTerm(aConstraintTermElement);
            
    }
    
    public void equalpElement(Element elem1,int hier1,int hier2, int dof){
        ConstraintTermElement aConstraintTermElement = new 
                        ConstraintTermElement(elem1,elem1.getNodeHier(hier1).getID(),dof, 1.0);
        this.addElemConstraintTerm(aConstraintTermElement);
        aConstraintTermElement = new 
                        ConstraintTermElement(elem1,elem1.getNodeHier(hier2).getID(),dof, -1.0);
        this.addElemConstraintTerm(aConstraintTermElement);
            
    }
    
    public void minimumPotential(Node theNode, int dof){
        Element elem;
        for(int i=0;i<theNode.getConnectedElementsIds().length;i++){
            elem=theNode.getConnectedElement(theNode.getConnectedElementsIds()[i]);
            double cof=1.; //if(i==0)cof=-1.;
            for(int j=1;j<=elem.getNumNodes();j++){
                ConstraintTermElement aConstraintTermElement = new 
                        ConstraintTermElement(elem,elem.getNodeHier(j).getID(),dof, cof*0.5*elem.getShapeFunctionProduct(dof, theNode.getID(), elem.getNodeHier(j).getID()));
                this.addElemConstraintTerm(aConstraintTermElement);
            }
        }
    }


    public void Print(){
         System.out.println("Costraint equation id : "+this.id+", num of NTerms = "+numOfNTerms+", num of ETerms = "+numOfETerms);
         if(Constraintvalue!=null)System.out.println("constant value (length = "+Constraintvalue.length+") of equation Constraintvalue[0] = "+Constraintvalue[0]);
         for(Iterator<ConstraintTerm> nit=this.getNodeConstraintTerm().iterator() ; nit.hasNext();){
             ConstraintTerm NodeConstraintTerm = nit.next();
             NodeConstraintTerm.Print();
         }
         for(Iterator<ConstraintTermElement> nit=this.getElementConstraintTerm().iterator() ; nit.hasNext();){
             ConstraintTermElement aConstraintTermElement = nit.next();
             aConstraintTermElement.Print();
         }
         System.out.println();
     }
    
    public int getNumOfConstraintEquations(){
        return numofConstraintEquations;
    }

}
