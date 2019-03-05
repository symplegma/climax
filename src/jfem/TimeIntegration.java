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
package jfem;

import java.util.Iterator;
import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
abstract public class TimeIntegration extends Analysis{
    protected boolean printinfo=false;
    
    protected void updateInitialAcceleration(){
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateAccls(theSOE.getX());
        }
    }
    
    protected void commit2(){
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitDisps();
            theNode.commitLoadCaseDisp(this.currentLC, step);
            if(Analysis2File)theNode.printDisps(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitVelcs();
            if(Analysis2File)theNode.printVelcs(this);
        }
        if(Analysis2File)this.AnalysisFile.println();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.commitAccls();
            if(Analysis2File)theNode.printAccls(this);
        }
    }
     
     public void ClearAccumulatedDips(){
         for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.clearAccumulatedDips();
        } 
     }
     
     protected void addM_v(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM_v();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getM_v();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            }
        }
        
    }
    
    protected void addM_a(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM_a();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getM_a();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            } 
        }

    }
    
    protected void addM_u(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM_u();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getM_u();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            } 
        }
    }
    
    protected void addC_u(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getM_u();
            theSOE.addToFext(mat, aInt, coef*a0);
            mat=elem.getK_u();
            theSOE.addToFext(mat, aInt, coef*a1);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                int[] aInt ;
                aInt=elem.getFtable();
                mat=elem.getM_u();
                theSOE.addToFext(mat, aInt, coef*a0);
                mat=elem.getK_u();
                theSOE.addToFext(mat, aInt, coef*a1);
            } 
        }
    }
    
    protected void addC_upre(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getM_upre();
            theSOE.addToFext(mat, aInt, coef*a0);
            mat=elem.getK_upre();
            theSOE.addToFext(mat, aInt, coef*a1);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                int[] aInt ;
                aInt=elem.getFtable();
                mat=elem.getM_upre();
                theSOE.addToFext(mat, aInt, coef*a0);
                mat=elem.getK_upre();
                theSOE.addToFext(mat, aInt, coef*a1);
            } 
        }
    }
    
    protected void addC_v(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getM_v();
            theSOE.addToFext(mat, aInt, coef*a0);
            mat=elem.getK_v();
            theSOE.addToFext(mat, aInt, coef*a1);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                int[] aInt ;
                aInt=elem.getFtable();
                mat=elem.getM_v();
                theSOE.addToFext(mat, aInt, coef*a0);
                mat=elem.getK_v();
                theSOE.addToFext(mat, aInt, coef*a1);
            } 
        }
    }
    
    protected void addC_a(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            int[] aInt ;
            aInt=elem.getFtable();
            mat=elem.getM_a();
            theSOE.addToFext(mat, aInt, coef*a0);
            mat=elem.getK_a();
            theSOE.addToFext(mat, aInt, coef*a1);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                int[] aInt ;
                aInt=elem.getFtable();
                mat=elem.getM_a();
                theSOE.addToFext(mat, aInt, coef*a0);
                mat=elem.getK_a();
                theSOE.addToFext(mat, aInt, coef*a1);
            } 
        }
    }
    
    protected void addM_upre(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM_upre();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getM_upre();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            } 
        }
    }
    
    protected void addK_u(double coef){
        for(Iterator<Element> it=theDomain.getElements().values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getK_u();
            int[] aInt ;
            aInt=elem.getFtable();
            theSOE.addToFext(mat, aInt, coef);
        }
        if (this.theImposer.TypeOfImposer!=1){
            for(Iterator<ConstraintElement> it=theDomain.getConstraintElements().values().iterator(); it.hasNext();){
                ConstraintElement elem = it.next();
                AbstractMatrix mat ;
                mat=elem.getK_u();
                int[] aInt ;
                aInt=elem.getFtable();
                theSOE.addToFext(mat, aInt, coef);
            } 
        }
    }
    
    @Override
    public void formLeft() {
        theSOE.clearK();
        theSOE.clearM();
        this.formStiffnessMatrix();
        this.formMassMatrix();
        theImposer.imposeLeft(this);
    }
    
    @Override
    public void update(){
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.updateDispsStep(theSOE.getX(),1.0);
            theNode.accumulateDisps(theSOE.getX());
        }
    }
     
     protected void formVelocity(){
         theSOE.clearXv();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            AbstractMatrix mat ;
            mat=theNode.getVelcsTrial_mat();
            int[] aInt ;
            aInt=theNode.getFtable();
            theSOE.addToXv(mat, aInt,1.0);
        } 
    }
     
     protected void formAcceleration(){
         theSOE.clearXa();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            AbstractMatrix mat ;
            mat=theNode.getAcclsTrial_mat();
            int[] aInt ;
            aInt=theNode.getFtable();
            theSOE.addToXa(mat, aInt,1.0);
        } 
    }
     
     protected void formDisplacement(){
         theSOE.clearXu();
        for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
            Node theNode = it.next();
            AbstractMatrix mat ;
            //mat=theNode.getDispsConvg_mat();
            mat=theNode.getDispsAccumulated_mat();
            int[] aInt ;
            aInt=theNode.getFtable();
            theSOE.addToXu(mat, aInt,1.0);
        } 
    }
}
