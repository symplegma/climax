/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;
import jmat.AbstractMatrix;
import java.util.Iterator;
import javax.swing.JOptionPane;
/**
 *
 * @author pchr
 */
public class TransientAnalysis2 extends TransientAnalysis{
    AbstractMatrix temp;
    AbstractMatrix res;
    
    // constructor
    public TransientAnalysis2(Domain theDomain, double TimeStep, int TotalSteps){
        this.putDomain(theDomain);
        this.steps=TotalSteps;
        this.StepLength=TimeStep;
        this.temp = new AbstractMatrix(theDomain.getSumDOFs(),theDomain.getSumDOFs());
        
    }
    
    // methods
    protected void formLeft(){
        this.theSOE.zeroA();
        int row=0;
        int nd=0;
        int md=0;
        int clnm;
        double val;
        int ndof;
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();

            row+=2*theDomain.getuDOFs(); //row*=2;
            nd+=theDomain.getpDOFs();
            md+=theDomain.getuDOFs();
            ndof=theDomain.getFundamentalSolution().get_p_DOFs();

            // fix constraint equations into A
            // -------------------------------

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
            //System.out.println("constraint equations incorporated into left hand side matrix.");
            // constraint equations incorporated into A
            // fix boundary integral equations into A
            // --------------------------------------
            int[] nodeEFT;
            int numNodes=theDomain.getNumNodes();
            int numCol=0;
            for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                numCol+=1;
                Node colNode = it.next();
                nodeEFT=colNode.getuEFTable();
                AbstractMatrix mat;
                int[] elemEFT;
                for(Iterator<Element>et=theDomain.theElements.values().iterator(); et.hasNext();){
                    Element patchElem = et.next();
                    mat=patchElem.getGu(colNode, theDomain);
                    elemEFT=patchElem.getpEFTable(theDomain);
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, -mat.get(irow, iclm));
                        }
                    }
                    mat=patchElem.getGv(colNode, theDomain);
                    elemEFT=patchElem.getpEFTable(theDomain);
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1, -mat.get(irow, iclm));
                             }
                    }

                    mat=patchElem.getHres(colNode, theDomain);
                    elemEFT=patchElem.getuEFTable(theDomain);
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        int ii=0;
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            if(ii==ndof){ii=0;}
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1, nodeEFT[ii]-1, -mat.get(irow, iclm));

                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1+md, nodeEFT[ii]-1+md, -mat.get(irow, iclm));

                            ++ii;
                        }
                    }
                    mat=patchElem.getH(colNode, theDomain);
                    elemEFT=patchElem.getuEFTable(theDomain);
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        int ii=0;
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            if(ii==ndof){ii=0;}
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, mat.get(irow, iclm));

                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1+md, mat.get(irow, iclm));

                            ++ii;
                        }
                    }
                    mat=patchElem.getHdif(colNode, theDomain);
                    elemEFT=patchElem.getuEFTable(theDomain);
                    for(int irow=0; irow<nodeEFT.length; irow++){
                        int ii=0;
                        for (int iclm=0; iclm<elemEFT.length; iclm++){
                            if(ii==ndof){ii=0;}
                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, mat.get(irow, iclm));

                            this.theSOE.addA(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1+md, mat.get(irow, iclm));

                            ++ii;
                        }
                    }
                }
                //System.out.print(numCol);System.out.print(" from ");System.out.print(numNodes);
                //System.out.println(" nodes have been collocated.");
            }
        }
        
    }
    
    protected void formRight(){
        theSOE.zeroB();
        int row=0;
        double val;
        int whichStep;
        int m;
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();

            row+=2*theDomain.getuDOFs();

            m=FundamentalSolution.theFSdata.getMstep();
            //for(int i=1; i<m; i++){
            for(int i=0; i<m; i++){
                this.temp.init();
                int[] nodeEFT;
                int md=theDomain.getuDOFs();
                int ndof=theDomain.getFundamentalSolution().get_p_DOFs();

                FundamentalSolution.theFSdata.setNstep(i+1);
                //FundamentalSolution.theFSdata.setNstep(i);
                FundamentalSolution.theFSdata.setTimePos(0);
                for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                    Node colNode = it.next();
                    nodeEFT=colNode.getuEFTable();
                    AbstractMatrix mat;
                    int[] elemEFT;
                    for(Iterator<Element>et=theDomain.theElements.values().iterator(); et.hasNext();){
                        Element patchElem = et.next();
                        mat=patchElem.getGu(colNode, theDomain);
                        elemEFT=patchElem.getpEFTable(theDomain);
                        for(int irow=0; irow<nodeEFT.length; irow++){
                            for (int iclm=0; iclm<elemEFT.length; iclm++){
                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, -mat.get(irow, iclm));
                            }
                        }

                        mat=patchElem.getGv(colNode, theDomain);
                        elemEFT=patchElem.getpEFTable(theDomain);
                        for(int irow=0; irow<nodeEFT.length; irow++){
                            for (int iclm=0; iclm<elemEFT.length; iclm++){
                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1, -mat.get(irow, iclm));
                                 }
                        }

                        mat=patchElem.getH(colNode, theDomain);
                        elemEFT=patchElem.getuEFTable(theDomain);
                        for(int irow=0; irow<nodeEFT.length; irow++){
                            int ii=0;
                            for (int iclm=0; iclm<elemEFT.length; iclm++){
                                if(ii==ndof){ii=0;}
                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, mat.get(irow, iclm));

                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1+md, mat.get(irow, iclm));

                                ++ii;
                            }
                        }
                    }
                }
                //********************************************
                //System.out.print("m1 ");System.out.print("n= ");System.out.println(FundamentalSolution.theFSdata.getNstep());
                //temp.print(12, 6);
                //********************************************
                FundamentalSolution.theFSdata.setNstep(i);
                //FundamentalSolution.theFSdata.setNstep(i-1);
                FundamentalSolution.theFSdata.setTimePos(1);
                for(Iterator<Node>it=theDomain.theNodes.values().iterator(); it.hasNext();){
                    Node colNode = it.next();
                    nodeEFT=colNode.getuEFTable();
                    AbstractMatrix mat;
                    int[] elemEFT;
                    for(Iterator<Element>et=theDomain.theElements.values().iterator(); et.hasNext();){
                        Element patchElem = et.next();
                        mat=patchElem.getGu(colNode, theDomain);
                        elemEFT=patchElem.getpEFTable(theDomain);
                        for(int irow=0; irow<nodeEFT.length; irow++){
                            for (int iclm=0; iclm<elemEFT.length; iclm++){
                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, -mat.get(irow, iclm));
                            }
                        }

                        mat=patchElem.getGv(colNode, theDomain);
                        elemEFT=patchElem.getpEFTable(theDomain);
                        for(int irow=0; irow<nodeEFT.length; irow++){
                            for (int iclm=0; iclm<elemEFT.length; iclm++){
                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1, -mat.get(irow, iclm));
                            }
                        }

                        mat=patchElem.getH(colNode, theDomain);
                        elemEFT=patchElem.getuEFTable(theDomain);
                        for(int irow=0; irow<nodeEFT.length; irow++){
                            int ii=0;
                            for (int iclm=0; iclm<elemEFT.length; iclm++){
                                if(ii==ndof){ii=0;}
                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1, elemEFT[iclm]-1, mat.get(irow, iclm));

                                this.temp.addVal(nodeEFT[irow]-theDomain.getpDOFs()-1+md, elemEFT[iclm]-1+md, mat.get(irow, iclm));

                                ++ii;
                            }
                        }
                    }
                }
                //********************************************
                //System.out.print("m2 ");System.out.print("n= ");System.out.println(FundamentalSolution.theFSdata.getNstep());
                //temp.print(12, 6);
                //********************************************
                this.res=theDomain.getResponse(i);
                //********************************************
                //System.out.println("res");
                //res.print(12, 6);
                //********************************************
                this.res=temp.times(res);

                for(int j=0; j<theDomain.getSumDOFs(); j++){
                    this.theSOE.addB(j, -res.get(j, 0));
                }

            }


            whichStep=FundamentalSolution.theFSdata.getMstep();
            for(Iterator<ConstraintEquation>
                    it=theDomain.getConstraintEquations().values().iterator();
                        it.hasNext();){
                ConstraintEquation theConstraintEquation = it.next();
                val=theConstraintEquation.getConstraintvalue(whichStep);
                this.theSOE.setB(row, val);
                ++row;
            }
        }
        
    }

    @Override
    public void init() {
        int sumDOFs=0;
        int BIEs=0;
        int CEs=0;
        int IEs=this.theInterfaceEquations.size();
        int numDomains=this.theDomains.size();
        
        if(numDomains>1 && IEs>0){
            // to be constucted
        }else{
            for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
                sumDOFs=0; BIEs=0; CEs=0;
                Domain theDomain = dit.next();
                sumDOFs+=theDomain.getSumDOFs();
                BIEs+=2*theDomain.getuDOFs();
                CEs+=theDomain.getNumberOfConstraintEquations();
                if( (BIEs+CEs+IEs) != sumDOFs){
                    JOptionPane.showMessageDialog(null,"Domain with id = "+theDomain.getID()+'\n'+"BIEs="+BIEs+", CEs="+CEs+", IEs="+IEs+", sumDOFs="+sumDOFs+'\n'+
                            "there are "+(sumDOFs-BIEs-CEs-IEs)+" missed.",
                            "error: program will exit",JOptionPane.ERROR_MESSAGE);
                    System.exit(301);
                }
            }
        }
        
        
        /*for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();

            sumDOFs+=theDomain.getSumDOFs();
            BIEs=+2*theDomain.getuDOFs();
            CEs+=theDomain.getNumberOfConstraintEquations();
            theDomain.updateNodes(theSOE.getX(),0);
        }
        
        if( (BIEs+CEs) != sumDOFs){
            JOptionPane.showMessageDialog(null, "BIEs="+BIEs+", CEs="+CEs+", sumDOFs="+sumDOFs,
                    "error: program will exit",JOptionPane.ERROR_MESSAGE);
            System.exit(303);
        }*/
        
        this.theSOE=new SOE();
        this.theSOE.init(sumDOFs,this.steps);
        
    }


    @Override
    public void run() {
        long lngCurrent=0;
        long lngPrevious=0;
        FundamentalSolution.theFSdata.setTimeStep(this.StepLength);
        for(int i=1; i<this.steps; i++){
            this.currentStep=i;
            FundamentalSolution.theFSdata.setMstep(i);
            lngCurrent=System.currentTimeMillis();
            if(i==1){System.out.print(0.+": ");}else{System.out.print((lngCurrent-lngPrevious)/1000.+": ");}
            System.out.print("m=  ");System.out.print(i);
            System.out.print(" from total of ");System.out.print(this.steps);System.out.print(" steps. ");
            System.out.println();
            
            FundamentalSolution.theFSdata.setNstep(i);
            FundamentalSolution.theFSdata.setTimePos(1);
            if(i==1){
                formLeft();
                //System.out.println("A");
                //theSOE.getA().print(20, 12);
            }
            formRight();
            //System.out.println("B");
            //theSOE.getB().print(12, 6);
            theSOE.solve();
            //System.out.println("X");
            //theSOE.getX().print(12, 6);
            for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
                Domain theDomain = dit.next();
                theDomain.updateNodes(theSOE.getX(),i);
            }
            lngPrevious=lngCurrent;
        }
    }

}
