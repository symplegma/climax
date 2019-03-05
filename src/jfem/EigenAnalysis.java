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
/**
 *
 * @author pchr
 */
public class EigenAnalysis extends Analysis{
    
    // constructor
    public EigenAnalysis(){
        this(1);
    }
    
    public EigenAnalysis(int typeofBCImposer){
        theDomain = new Domain();
        this.theType=2;
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            case 2: theImposer = new ClassicLargeMass(); break;
            case 3: theImposer = new GeneralizedLargeMass(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        
    }
    
    public EigenAnalysis(Domain theDomain){
        this(theDomain,1);
    }
    
    public EigenAnalysis(Domain theDomain, int typeofBCImposer){
        this.theDomain = theDomain;
        this.theType=2;
        switch(typeofBCImposer){
            case 1: theImposer = new StaticImposer(); break;
            case 2: theImposer = new ClassicLargeMass(); break;
            case 3: theImposer = new GeneralizedLargeMass(); break;
            default:System.out.println("unkwon type of BCImposer");
            System.exit(1); break;
        }
        if(Analysis2File)this.printGenerals();
        
    }
    
    // methods
    // yet the (int) LC does not used
    public int analyse(int LC){
        int order=this.theDomain.getNdofs();
        theSOE = new SOE(order,order);
        theDomain.clearDomain();
        theDomain.setEigenLoadCases(order);
        this.formLeft();
        this.run();
        this.update();
        //this.theSOE.printSOEeigen();
        return 0;
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
    public void formRight() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void run() {
        this.theSOE.EigenRun();
    }

    @Override
    // for EigenAnalysis update incorporates commit
    public void update(){
//        int order=this.theDomain.getNdofs();
        int order=theSOE.getEigIndexLength();
        double[] eig_vals = new double[order];
        for(int i=0;i<order;i++){
            eig_vals[i]=theSOE.getEigenValues(i);
        }
//        eig_vals=theSOE.getEigenValues();
        theDomain.setEigenValues(eig_vals);
        if(Analysis2File)this.AnalysisFile.println();
        if(Analysis2File)AnalysisFile.println("number of eigval"+"     "+"eigenvalues     "+"     "+"eigenfrequencies"+"     "+"eigenperiods    ");
        if(Analysis2File)AnalysisFile.println("----------------"+"     "+"----------------"+"     "+"----------------"+"     "+"----------------");
        for(int i=0; i<eig_vals.length; i++){
             //System.out.println(Math.sqrt(eig_vals[i]));
             if(Analysis2File){
                AnalysisFile.format("%16d",i+1);
                AnalysisFile.print("     ");
                AnalysisFile.format("%16.8f",eig_vals[i]); AnalysisFile.print("     ");
                AnalysisFile.format("%16.8f",Math.sqrt(eig_vals[i])); AnalysisFile.print("     ");
                AnalysisFile.format("%16.8f",2.0*Math.PI/Math.sqrt(eig_vals[i])); AnalysisFile.println();
             }
         }
        
        for(int i=0; i<order; i++){
            if(Analysis2File){
                this.AnalysisFile.println();
                this.AnalysisFile.println("EigenVector      : "+(i+1));
                this.AnalysisFile.println("______________________");
                this.AnalysisFile.println();
            }
            for(Iterator<Node> it=theDomain.getNodes().values().iterator(); it.hasNext();){
                Node theNode = it.next();
                theNode.updateDisps(theSOE.getEigenVector(i));
                theNode.commitDisps();
                theNode.commitLoadCaseDisp(-i-1,1);
                if(Analysis2File)theNode.printDisps(this);
            }
            
            // print for trackers
            for(Iterator<Tracker> it=this.theTrackers.values().iterator(); it.hasNext();){
                Tracker theTracker = it.next();
                theTracker.write(i+1);
            }
        }
        
    }

    @Override
    public void commit() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void printGenerals() {
        String aString="Analysis with id: "+this.id;
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        aString="Type of Algorithm: "+"Eigen-values Analysis";
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        switch(this.theImposer.getType()){
            case 1: aString="Type of Boundary Conditions Imposer: "+"Penalty"; break;
            default: aString="Type of Boundary Conditions Imposer: "+"Unknown"; break;
        }
        if(Analysis2File)this.AnalysisFile.println(aString);
        
        aString="Domain of analysis id: "+this.theDomain.getID();
        if(Analysis2File)this.AnalysisFile.println(aString);
        aString="Number of degrees of freedom: "+this.theDomain.getNdofs();
        if(Analysis2File)this.AnalysisFile.println(aString);
        if(Analysis2File)this.AnalysisFile.println("===================================================");
    }
}
