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

//import java.util.Iterator;
import jmat.AbstractMatrix;
import java.util.Iterator;
import javax.swing.JOptionPane;

/**
 *
 * @author pchr
 */
public class StaticAnalysis extends Analysis{
    private int  numConstructedLHS=0;
    private AbstractMatrix tmp;
    private AbstractMatrix TemperatuteVector;
    private boolean Thermoelastic=false;
    
    // constructor
    public StaticAnalysis(Domain theDomain){
        this.putDomain(theDomain);
        this.steps=1;
    }

    public StaticAnalysis(){
        this.steps=1;
    }
    
    // methods
    
    protected void formLeft(int whichStep){
        if(Thermoelastic)TemperatuteVector=new AbstractMatrix(this.theSOE.getA().getColumnDimension(),1, 0.0);
        numConstructedLHS+=1;
        int row=0;
        int totDrow=0;
        int prevpdofs=0;
        int ndof=0;
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();
            totDrow+=theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+theDomain.getRigidRemovalDOFs();
            // fix constraint equations into A
            // -------------------------------
            row+=theDomain.getuDOFs();
            int clnm;
            double val;
            ndof=theDomain.getFundamentalSolution().get_p_DOFs();
            for(Iterator<ConstraintEquation>
                    it=theDomain.getConstraintEquations().values().iterator();
                        it.hasNext();){
                ConstraintEquation theConstraintEquation = it.next();
                for(Iterator<ConstraintTerm> nit=theConstraintEquation.getNodeConstraintTerm().iterator() ; nit.hasNext();){
                    ConstraintTerm NodeConstraintTerm = nit.next();
                    clnm=NodeConstraintTerm.getDOF()-1;
                    val=NodeConstraintTerm.getCoef();
                    this.theSOE.setA(row,clnm, val);
                }
                for(Iterator<ConstraintTermElement> nit=theConstraintEquation.getElementConstraintTerm().iterator() ; nit.hasNext();){
                    ConstraintTermElement aConstraintTermElement = nit.next();
                    clnm=aConstraintTermElement.getDOF(ndof)-1;
                    val=aConstraintTermElement.getCoef();
                    this.theSOE.setA(row,clnm, val);
                }
                ++row;
            }

            // fix boundary integral equations into A
            // --------------------------------------


            int[] nodeEFT;
            int cc=0;
            for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                ++cc;
                //System.out.println("Domain "+theDomain.getID()+", collocated node "+cc+"/"+theDomain.getNumNodes());
                Node colNode = it.next();
                nodeEFT=colNode.getuEFTable();
                AbstractMatrix mat;
                AbstractMatrix vec=null;
                int[] elemEFT;
                for(Iterator<Element>et=theDomain.theElements.values().iterator(); et.hasNext();){
//                    if(theDomain.getID()==2){
//                        System.out.println("2");
//                    }
                    Element patchElem = et.next();
                    mat=patchElem.getGu(colNode, theDomain);
                    elemEFT=patchElem.getpEFTable(theDomain);
                    if(Thermoelastic)vec=patchElem.getNornalVec(theDomain);
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()+prevpdofs-1, elemEFT[iclm]-1, -mat.get(irow, iclm));
                        }
                        if(Thermoelastic){
                            val=0.;
                            for(int iclm=0;iclm<elemEFT.length;iclm++)val+=mat.get(irow, iclm)*vec.get(iclm, 0)*((ElasticMat)theDomain.getMaterial()).getExtendedThermalCoef()*theDomain.getUniformTempChange();
                            TemperatuteVector.addVal(nodeEFT[irow]-theDomain.getpDOFs()+prevpdofs-1, 0, val);
                        }
                    }
                    
                    mat=patchElem.getH(colNode, theDomain);
//                    mat=tmp.transpose().times(patchElem.getH(colNode, theDomain).times(tmp));
                    elemEFT=patchElem.getuEFTable(theDomain);for(int irow=0; irow<nodeEFT.length; irow++){
                        int ii=0;
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            if(ii==ndof){ii=0;}
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()+prevpdofs-1, elemEFT[iclm]-1, mat.get(irow, iclm));
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()+prevpdofs-1, nodeEFT[ii]-1, -mat.get(irow, iclm));
                            ++ii;
                        }
                    }
                }
            }
            // Consider infinite or semi_ifinite region (it might need more study)
            if(theDomain.getExtension()!=Domain.Extension.FINITE){
                double DVAL=1.0;
                if(theDomain.getExtension()==Domain.Extension.SEMI_INFINITE){DVAL=0.5;}
                for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                    Node colNode = it.next();
                    nodeEFT=colNode.getuEFTable();
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        this.theSOE.getA().putVal(nodeEFT[irow]-theDomain.getpDOFs()+prevpdofs-1, nodeEFT[irow]-1, DVAL+this.theSOE.getA().get(nodeEFT[irow]-theDomain.getpDOFs()+prevpdofs-1, nodeEFT[irow]-1));
                    }
                }
            }
            
            // Fix Rigid Body Removal Equations
            // --------------------------------------
            // Part W (On the removal of rigid body motion in the solution of elastostatic problems by direct BEM, 
            // Blasquez, Mantic, Paris, Cañas, eq.(14a) && qe.(27))
            if(theDomain.getRigidRemovalDOFs()>0){
                int counter=0;
                for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                    Node colNode = it.next();
                    nodeEFT=colNode.getuEFTable();
                    if(theDomain.getFixedDofs()!=null && theDomain.getForWMatrixRBM()==false){
                        for(int i=0;i<nodeEFT.length;i++){
                            for(int j=0;j<theDomain.getFixedDofs().length;j++){
                                if(nodeEFT[i]==theDomain.getFixedDofs()[j]){
                                    this.theSOE.setA(nodeEFT[i]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+counter, 1.);
                                    counter+=1;
                                }
                            }
                        }
                    }else{
                        for(int i=0;i<theDomain.getFixedModes().length;i++){
                            if(theDomain.getFixedModes()[i]==0)this.theSOE.setA(nodeEFT[0]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+i, 1.);
                            if(theDomain.getFixedModes()[i]==1)this.theSOE.setA(nodeEFT[1]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+i, 1.);
                            if(theDomain.getFixedModes()[i]==2){
                                this.theSOE.setA(nodeEFT[0]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+i, -colNode.getCoordinates()[1]);
                                this.theSOE.setA(nodeEFT[1]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+i, colNode.getCoordinates()[0]);
                            }
                        }
                        //this.theSOE.setA(nodeEFT[0]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs(), 1.);
                        //this.theSOE.setA(nodeEFT[1]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+1, 1.);
                        //this.theSOE.setA(nodeEFT[0]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+2, -colNode.getCoordinates()[1]);
                        //this.theSOE.setA(nodeEFT[1]-theDomain.getpDOFs()-1, theDomain.getpDOFs()+theDomain.getuDOFs()+2, colNode.getCoordinates()[0]);
                    }
                    
                }
            }
            // Part V.tranpose (On the removal of rigid body motion in the solution of elastostatic problems by direct BEM, 
            // Blasquez, Mantic, Paris, Cañas, eq.(14b) && qe.(27))
            if(theDomain.getRigidRemovalDOFs()>0){
                int counter=0;
                for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                    Node colNode = it.next();
                    nodeEFT=colNode.getuEFTable();
                    if(theDomain.getFixedDofs()!=null){
                        for(int i=0;i<nodeEFT.length;i++){
                            for(int j=0;j<theDomain.getFixedDofs().length;j++){
                                if(nodeEFT[i]==theDomain.getFixedDofs()[j]){
                                    this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+counter, nodeEFT[i]-1, 1.);
                                    counter+=1;
                                }
                            }
                        }
                    }else{
                        double mat;
                        int[] elemEFT;
                        for(Iterator<Element>et=theDomain.theElements.values().iterator(); et.hasNext();){
                            Element patchElem = et.next();
                            mat=patchElem.getM(colNode, theDomain);
                            elemEFT=patchElem.getuEFTable(theDomain);
                            for(int irow=0; irow<nodeEFT.length; irow++){
                                for (int iclm=0; iclm<elemEFT.length; iclm++){
                                    for(int i=0;i<theDomain.getFixedModes().length;i++){
                                        if(theDomain.getFixedModes()[i]==0)this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+i, nodeEFT[0]-1, mat);
                                        if(theDomain.getFixedModes()[i]==1)this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+i, nodeEFT[1]-1, mat);
                                        if(theDomain.getFixedModes()[i]==2){
                                            this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+i, nodeEFT[0]-1, -colNode.getCoordinates()[1]*mat);
                                            this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+i, nodeEFT[1]-1, colNode.getCoordinates()[0]*mat);
                                        }
                                    }
                                    //this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations(), nodeEFT[0]-1, mat);
                                    //this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+1, nodeEFT[1]-1, mat);
                                    //this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+2, nodeEFT[0]-1, -colNode.getCoordinates()[1]*mat);
                                    //this.theSOE.addA(prevpdofs+theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations()+2, nodeEFT[1]-1, colNode.getCoordinates()[0]*mat);
                                }
                            }
                        }
                    }
                    
                }
            }
            // rigid body motions have been "probably" constrained
            
            
            prevpdofs+=theDomain.getNumberOfConstraintEquations()-theDomain.getpDOFs()-theDomain.getRigidRemovalDOFs();
        }
        
        // fix interface equations into A
        // --------------------------------------
        for(Iterator<InterfaceEquation>
                    it=this.theInterfaceEquations.values().iterator();
                        it.hasNext();){
            InterfaceEquation theInterfaceEquation = it.next();
            for(Iterator<InterfaceTerm> nit=theInterfaceEquation.getNodeInterfaceTerm().iterator() ; nit.hasNext();){
                InterfaceTerm NodeInterfaceTerm = nit.next();
                int clnm=NodeInterfaceTerm.getDOF()-1;
                double val=NodeInterfaceTerm.getCoef();
                this.theSOE.setA(totDrow,clnm, val);
            }
            for(Iterator<InterfaceTermElement> nit=theInterfaceEquation.getElementConstraintTerm().iterator() ; nit.hasNext();){
                InterfaceTermElement aInterfaceTermElement = nit.next();
                int clnm=aInterfaceTermElement.getDOF(ndof)-1;
                double val=aInterfaceTermElement.getCoef();
                this.theSOE.setA(totDrow,clnm, val);
            }
            ++totDrow;
        }
    }
    
    protected void formRight(int whichStep){
        int row=0;
        int prevDomainID=0;
        int totDrow=0;
        double val;
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();
            totDrow+=theDomain.getuDOFs()+theDomain.getNumberOfConstraintEquations();

            row+=theDomain.getuDOFs();
            //if(prevDomainID!=0){row+=this.theDomains.get(prevDomainID).getNumberOfConstraintEquations();}
            
            for(Iterator<ConstraintEquation>
                    it=theDomain.getConstraintEquations().values().iterator();
                        it.hasNext();){
                ConstraintEquation theConstraintEquation = it.next();
                val=theConstraintEquation.getConstraintvalue(whichStep);
                this.theSOE.setB(row, val);
                ++row;
            }
            prevDomainID=theDomain.getID();
        }
        row=0;
        for(Iterator<InterfaceEquation>
                    it=this.theInterfaceEquations.values().iterator();
                        it.hasNext();){
            InterfaceEquation theInterfaceEquation = it.next();
            val=theInterfaceEquation.getInterfacevalue(whichStep);
            this.theSOE.setB(row+totDrow, val);
            ++row;
        }
        if(this.Thermoelastic){
            for(int i=0;i<this.theSOE.getA().getColumnDimension();i++){
                this.theSOE.addB(i, this.TemperatuteVector.get(i, 0));
            }
        }
    }
    
    @Override
    public void init() {
        int sumDOFs=0;
        int BIEs=0;
        int CEs=0;
        int RRs=0;
        int IEs=this.theInterfaceEquations.size();
        int numDomains=this.theDomains.size();
        if(numDomains>1 && IEs>0){
            sumDOFs=0; BIEs=0; CEs=0; RRs=0;
            for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
                Domain theDomain = dit.next();
                sumDOFs+=theDomain.getSumDOFs();
                BIEs+=theDomain.getuDOFs();
                CEs+=theDomain.getNumberOfConstraintEquations();
                RRs+=theDomain.getRigidRemovalDOFs();
                //if( (BIEs+CEs+IEs+RRs) != sumDOFs)theDomain.checkConstraints();
                //System.out.println("Domain with id: "+theDomain.getID()+" dofs: "+theDomain.getBeginDOFs()+" udofs: "+theDomain.getuDOFs()+" pdofs: "+theDomain.getpDOFs());
            }
            if( ((BIEs+CEs+IEs+RRs) != sumDOFs)){
                for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
                    dit.next().checkConstraints();
                }
                this.checkInterfaceConstraints();
                JOptionPane.showMessageDialog(null,"BIEs="+BIEs+", CEs="+CEs+", IEs="+IEs+", RRs="+RRs+", sumDOFs="+sumDOFs+'\n'+
                        "there are "+(sumDOFs-BIEs-CEs-IEs-RRs)+" missed.",
                        "error: program will exit",JOptionPane.ERROR_MESSAGE);
                System.err.println("error (302) in:"+this.getClass().toString()+", method:"+"init");
                //System.exit(302);
            }
        }else{
            for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
                sumDOFs=0; BIEs=0; CEs=0;
                Domain theDomain = dit.next();
                sumDOFs+=theDomain.getSumDOFs();
                BIEs+=theDomain.getuDOFs();
                CEs+=theDomain.getNumberOfConstraintEquations();
                RRs+=theDomain.getRigidRemovalDOFs();
                if( ((BIEs+CEs+IEs+RRs) != sumDOFs)&&(CEs!=0)) theDomain.checkConstraints();
                if( ((BIEs+CEs+IEs+RRs) != sumDOFs)&&(CEs!=0)) {
                    JOptionPane.showMessageDialog(null,"Domain with id = "+theDomain.getID()+'\n'+"BIEs="+BIEs+", CEs="+CEs+", IEs="+IEs+", RRs="+RRs+", sumDOFs="+sumDOFs+'\n'+
                            "there are "+(sumDOFs-BIEs-CEs-IEs-RRs)+" missed.",
                            "error: program will exit",JOptionPane.ERROR_MESSAGE);
                    System.err.println("error (301) in:"+this.getClass().toString()+", method:"+"init");
//                    System.exit(301);
                }
                if((BIEs+CEs+IEs+RRs) != sumDOFs) {
                    JOptionPane.showMessageDialog(null,"Domain with id = "+theDomain.getID()+'\n'+"BIEs="+BIEs+", CEs="+CEs+", IEs="+IEs+", RRs="+RRs+", sumDOFs="+sumDOFs+'\n'+
                            "there are "+(sumDOFs-BIEs-CEs-IEs-RRs)+" missed.",
                            "WARNING: ¡but without constraint equations!",JOptionPane.WARNING_MESSAGE);
                }
            }
        }
        this.theSOE=new SOE();
        this.theSOE.init(sumDOFs);
    }

    @Override
    public void run() {
        //System.out.println();System.out.println();
        //System.out.println("Left Hand Side matrix now being formulated.");
        long start = System.currentTimeMillis();
        formLeft(0);
        long elapsedTimeMillis = System.currentTimeMillis()-start;
        
        // Get elapsed time in seconds
        float elapsedTimeSec = elapsedTimeMillis/1000F;

        // Get elapsed time in minutes
        //float elapsedTimeMin = elapsedTimeMillis/(60*1000F);

        // Get elapsed time in hours
        //float elapsedTimeHour = elapsedTimeMillis/(60*60*1000F);

        // Get elapsed time in days
        //float elapsedTimeDay = elapsedTimeMillis/(24*60*60*1000F);

        
        start = System.currentTimeMillis();
        formRight(0);
        //this.theSOE.getA().print(12,12);
        //this.theSOE.getB().print(12,12);
        elapsedTimeMillis = System.currentTimeMillis()-start;
        elapsedTimeSec = elapsedTimeMillis/1000F;
        //System.out.println("Right Hand Side matrix has been formulated. "+elapsedTimeSec);
        //System.out.println("Algebraic system of equations are now being solved.");
        start = System.currentTimeMillis();
        theSOE.solve();
        elapsedTimeMillis = System.currentTimeMillis()-start;
        elapsedTimeSec = elapsedTimeMillis/1000F;
        //System.out.println("Algebraic system of equations have been solved. "+elapsedTimeSec);
        //System.out.println();System.out.println();
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();

            theDomain.updateNodes(theSOE.getX(),0);
            theDomain.updateRigidMultipliers(theSOE.getX(),0);
            /*System.out.println("Domain with id= "+theDomain.getID());
            System.out.println("==========================");
            System.out.println("           nodes          ");
            System.out.println("==========================");
            theDomain.printNodes(0);
            System.out.println("==========================");
            System.out.println("        elements          ");
            System.out.println("==========================");
            theDomain.printElements(0);
            if(theDomain.getRigidRemovalDOFs()>0){
                System.out.println("==========================");
                System.out.println("      rigid multipliers   ");
                System.out.println("==========================");
                for(int i=0;i<theDomain.getRigidRemovalDOFs();i++){
                    System.out.println("rigid multiplier : "+(i+1)+" has the value : "+theSOE.getX().get(theDomain.getSumDOFs()-theDomain.getRigidRemovalDOFs()+i, 0));
                }
            }*/
        }
    }

    public void run(int WhichStep, int Where2Save) {
        if(this.constructLHS){
            formLeft(0);
//            theSOE.getA().print(40, 40);
        }
        formRight(0);
        
//        theSOE.getB().print(20, 16);
        theSOE.solve();
        //theSOE.getX().print(20, 16);
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();
            if(Where2Save>=0){
                //theSOE.getA().print(20, 16);
                theDomain.updateNodes(theSOE.getX(),WhichStep,Where2Save);
                theDomain.updateRigidMultipliers(theSOE.getX(),WhichStep);
            }else{
                theDomain.updateNodesT(theSOE.getX());
            }
        }
    }
    
    public void setSolution(double coef,int fromWhichStep, int toWhichStep, int Where2Save){
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();
            theDomain.setResponse(coef, fromWhichStep, toWhichStep, Where2Save);
        }
    }
    
    public void setZeroSolution(int toWhichStep, int Where2Save){
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();
            theDomain.setZeroResponse(toWhichStep, Where2Save);
        }
    }
    
    public int getNumConstructedLHS(){return this.numConstructedLHS;}
    
    public void setThermoElastic(boolean v){this.Thermoelastic=v;}
    
    public boolean getThermoElastic(){return this.Thermoelastic;}
}
