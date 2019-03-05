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
import javax.swing.JOptionPane;

public class HarmonicAnalysis extends Analysis{
    private int  numConstructedLHS=0;
    private double[] omegas;
    
    public HarmonicAnalysis(Domain theDomain, double[] oms){
        this.theDomains.put(theDomain.getID(), theDomain);
        this.omegas=new double[oms.length];
        for(int i=0;i<oms.length;i++)omegas[i]=oms[i];
    }
    
    public HarmonicAnalysis(Domain theDomain, double om){
        this.theDomains.put(theDomain.getID(), theDomain);
        this.omegas=new double[1]; omegas[0]=om;
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
