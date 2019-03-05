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

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import jbem.SOE.Solver;

/**
 *
 * @author pchr
 */
abstract public class Analysis{
    protected Map<Integer,Domain> theDomains = new HashMap<Integer,Domain>();
    //Domain theDomain;
    SOE theSOE;
    int steps;
    double StepLength;
    protected Map<Integer,InterfaceEquation>
              theInterfaceEquations = new HashMap<Integer,InterfaceEquation>();
    protected boolean constructLHS=true;
    
    // constructor
    public Analysis(){}
    
    // common methods
    public void setStepLength(double len){
        this.StepLength=len;
    }

    public void putDomain(Domain aDomain){
        this.theDomains.put(aDomain.getID(), aDomain);
    }

    public void putInterfaceEquation(InterfaceEquation anInterfaceEquation){
        this.theInterfaceEquations.put(anInterfaceEquation.getID(), anInterfaceEquation);
    }

    public Map<Integer,InterfaceEquation> getInterfaceEquations(){
        return this.theInterfaceEquations;
    }

    public boolean getConstructLHS(){
        return this.constructLHS;
    }

    public void setConstructLHS(boolean t){
        this.constructLHS=t;
    }

    public void clearInterfaceEquations(){
        this.theInterfaceEquations.clear();
    }

    public void PrintInterfaceEquations(){
         System.out.println("---- Interface Equations Of Analysis, size = "+theInterfaceEquations.size());
                 System.out.println();
                 for(Iterator<InterfaceEquation>
                         it=this.theInterfaceEquations.values().iterator(); it.hasNext();){
                     InterfaceEquation theConstraintEquation = it.next();
                     theConstraintEquation.Print();
                 }

        System.out.println("-------------------------------------------------------------");
    }
    
    public void setSolver(Solver theSOLVER){
        this.theSOE.setEquationSolver(theSOLVER);
    }
    
    public void checkInterfaceConstraints(){
        for(Iterator<Domain> dt=this.theDomains.values().iterator(); dt.hasNext();){
            Domain aDomain = dt.next();
            for(Iterator<Node> nt=aDomain.getNodes().values().iterator(); nt.hasNext();){
                Node theNode = nt.next();
                System.out.println("For node : "+theNode.getID()+" ic equations :"+this.getNumOfInterfaceConstraintsOnNode(theNode.getID()));
            }
        }
    }
    
    public int getNumOfInterfaceConstraintsOnNode(int nodeid){
        int n=0;
        boolean exist=false;
        for(Iterator<InterfaceEquation>it=this.theInterfaceEquations.values().iterator(); it.hasNext();){
                     InterfaceEquation theConstraintEquation = it.next();
                     exist=false;
                     for(Iterator<InterfaceTerm> nt=theConstraintEquation.getNodeInterfaceTerm().iterator(); nt.hasNext();){
                         InterfaceTerm theNCT = nt.next();
                         if(theNCT.getNodeID()==nodeid)exist=true;
                     }
                     
                     for(Iterator<InterfaceTermElement> nt=theConstraintEquation.getElementConstraintTerm().iterator(); nt.hasNext();){
                         InterfaceTermElement theNCT = nt.next();
                         if(theNCT.getNodeID()==nodeid)exist=true;
                     }
                     if(exist)n+=1;
                 }
        return n;
    }
    
    public void printBEMmatrix(){
        this.theSOE.getA().print(12, 12);
    }
    
    public void printBEMsystem(){
        this.theSOE.getA().print(12, 12);
        this.theSOE.getB().print(12, 12);
        this.theSOE.getX().print(12, 12);
    }
    
    // abstract methods
    abstract public void init();
    abstract public void run();
}
