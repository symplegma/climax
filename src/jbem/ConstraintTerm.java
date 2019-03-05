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

/**
 *
 * @author pchr
 */
public class ConstraintTerm {
    protected Node theNode;
    protected int Variable; 
    // Variable is 0 for tractions, 1 for displacement, 2 for velocity
    protected int dof_sequence;
    protected double coef;
    private static int numberOfConstraintTerms = 0;
    private int id;
    
    // constructor
    public ConstraintTerm(){}
    
    public ConstraintTerm(Node theNode, int Variable, int dof_sequence, double coef){
        if(theNode==null)throw new NullPointerException();
        this.theNode=theNode;
        this.Variable=Variable;
        this.dof_sequence=dof_sequence;
        this.coef=coef;
        ++numberOfConstraintTerms;
        id=numberOfConstraintTerms;
    }
    
    public ConstraintTerm(Node theNode, int Variable, int dof_sequence){
        if(theNode==null)throw new NullPointerException();
        this.theNode=theNode;
        this.Variable=Variable;
        this.dof_sequence=dof_sequence;
        this.coef=1.0;
        ++numberOfConstraintTerms;
        id=numberOfConstraintTerms;
    }

    public int getID(){
        return this.id;
    }

    public int getNodeID(){
        return this.theNode.getID();
    }
    
    public int getDOF(){
        int dof;
        switch(Variable){
            case 1: dof = theNode.getuEFTable()[dof_sequence-1];break;
            case 2: dof = theNode.getvEFTable()[dof_sequence-1];break;
            default: dof=0; break;
        }
        return dof;
    }
    
    public double getCoef(){
        return this.coef;
    }

    public int getDOF_Sequence(){
        return this.dof_sequence;
    }
    
    public void Print(){
        String SVariable;
        switch(Variable){
            case 0: SVariable = " on traction"; break;
            case 1: SVariable = " on displacement"; break;
            case 2: SVariable = " on velocity"; break;
            default:SVariable = " ! unknown !"; break;
        }
        System.out.println("on Node: "+this.theNode.getID()+SVariable+", dof sequence: "+dof_sequence+", with coefficient = "+this.coef);
    }
    

}
