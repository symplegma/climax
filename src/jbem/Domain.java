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

import geom.Line2D;
import delaunay.Pnt;
import delaunay.Triangulation;
import geom.Point;
import geom.Shape;
import geom.Triangle;
import jmat.AbstractMatrix;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import mathman.DoubleFunction;

/**
 *
 * @author pchr
 */
public class Domain {
    protected int id;
    protected Map<Integer,ResultPoint> theResultPoints = new TreeMap<Integer,ResultPoint>();
    protected Map<Integer,Node> theNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,Element> theElements = new TreeMap<Integer,Element>();
    protected Map<Integer,ConstraintEquation> 
              theConstraintEquations = new TreeMap<Integer,ConstraintEquation>();
    protected Map<Integer,Node> theConnectedNodes = new TreeMap<Integer,Node>();
    protected Material theMaterial;
    protected FundamentalSolution theFundSol;
    protected int sumDOFS;
    protected int uDOFS;
    protected int pDOFS;
    private Extension extension=Extension.FINITE; // 0 finite, 1 infinite, 2 semi-infinite
    private AbstractMatrix res;
    private AbstractMatrix Lmat;
    private int beginDOF=0;
    private boolean isDissipativeSurface=false;
    private double[] eigenvalues;
    private AbstractMatrix eigenvectors;
    private int RigidBodyRemoval=0;
    private double[][] rigidMultipliers;
    private int fixeddofs[];
    private int fixedmodes[]; // 0 for horizontal, 1 for vertical, 2 for rotational
    private boolean fixM; // THIS MEANS THAT FOR THE V matrix the rigid body
                          // are going to be used even if point-wise "constraints" are given
    protected Map<Integer,Combo> theCombos = new TreeMap<Integer,Combo>();
    protected Map<Integer,State> theStates = new TreeMap<Integer,State>();
    private int responseID=0;
    private boolean auxiliaryFields=false;
    private boolean currentFieldAuxiliary=false;
    private double tol=1.e-12;
    
    private double uniformTempChange=0.;
    private static int numberOfDomains = 0;
    
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Z_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Z_coordinate=Double.POSITIVE_INFINITY;
    protected Map<Integer,Shape> theShapes = new TreeMap<Integer,Shape>();
    private boolean isTriangulized=false;
    
    // constructor
//    public Domain(){}
    
    public Domain(){
        id = ++numberOfDomains;
    }
    
    public void clsNumberOfDomains(){numberOfDomains=0;}
    
    public Domain(int id){
        this.id=id;
    }
    
    // methods
    public void setExtension(Extension ext){
        this.extension=ext;
    }

    public int getID(){
        return this.id;
    }
    
    public void addMaterial(Material aMaterial){
        this.theMaterial=aMaterial;
        if(aMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : aMaterial.getClass().toString().equals(ViscousMaterial.class.toString()))setAuxiliaryVariables();
    }
    
    public void putPoint(ResultPoint aPoint){
        aPoint.setInDomain(this);
        this.theResultPoints.put(aPoint.getID(), aPoint);
    }
    
    public void putNode(Node aNode){
        if(aNode.getCoordinates()[0]>maximum_X_coordinate)maximum_X_coordinate=aNode.getCoordinates()[0];
        if(aNode.getCoordinates()[0]<minimum_X_coordinate)minimum_X_coordinate=aNode.getCoordinates()[0];
        if(aNode.getCoordinates().length>1){
            if(aNode.getCoordinates()[1]>maximum_Y_coordinate)maximum_Y_coordinate=aNode.getCoordinates()[1];
            if(aNode.getCoordinates()[1]<minimum_Y_coordinate)minimum_Y_coordinate=aNode.getCoordinates()[1];
        }
        if(aNode.getCoordinates().length>2){
            if(aNode.getCoordinates()[2]>maximum_Z_coordinate)maximum_Z_coordinate=aNode.getCoordinates()[2];
            if(aNode.getCoordinates()[2]<minimum_Z_coordinate)minimum_Z_coordinate=aNode.getCoordinates()[2];
        }
        this.theNodes.put(aNode.getID(), aNode);
    }
    
    public void putCombo(Combo aCombo){
        aCombo.setResponseID(responseID++);
        this.theCombos.put(aCombo.getID(), aCombo);
    }
    
    public void putState(State aState){
        aState.setResponseID(responseID++);
        this.theStates.put(aState.getID(), aState);
    }
    
    public void putElement(Element anElement){
        this.theElements.put(anElement.getID(), anElement);
    }
    
    public void putConstraintEquation(ConstraintEquation anConstraintEquation){
        this.theConstraintEquations.put(anConstraintEquation.getID(), anConstraintEquation);
    }
    
    public Extension getExtension(){
        return this.extension;
    }
    
    public Node getNode(int id){
        return this.theNodes.get(id);
    }
    
    public ResultPoint getResultPoint(int id){
        return this.theResultPoints.get(id);
    }
    
    public Combo getCombo(int id){
        return this.theCombos.get(id);
    }
    
    public State getState(int id){
        return this.theStates.get(id);
    }
    
    public Element getElement(int id){
        return this.theElements.get(id);
    }
    
    public Map<Integer,Element> getElements(){
        return this.theElements;
    }

    public Map<Integer,Node> getNodes(){
        return this.theNodes;
    }
    
    public Map<Integer,ResultPoint> getResultPoints(){
        return this.theResultPoints;
    }
    
    public Map<Integer,Combo> getCombos(){
        return this.theCombos;
    }
    
    public Map<Integer,State> getStates(){
        return this.theStates;
    }

    public Map<Integer,Node> getConnectedNodes(){
        return this.theConnectedNodes;
    }
    
    public void setSpaceIntegrator(SpaceIntegrator SI){
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.setSpaceIntegrator(SI);
        }
    }
    
    public void setSpaceIntegrator_aux(boolean b){
        this.currentFieldAuxiliary=b;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.getSpaceIntegrator().setAuxiliaryField(b);
        }
    }
    
    private SpaceIntegrator getSpaceIntegrator(){
        SpaceIntegrator si=null;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            si=theElement.getSpaceIntegrator();
            break;
        }
        return si;
    }
    
    public boolean getCurrentFieldAuxiliary(){return this.currentFieldAuxiliary;}
    
    public void setFundamentalSolution(FundamentalSolution theFundSol){
        this.theFundSol=theFundSol;
        this.theFundSol.theFSdata.setMaterial(this.theMaterial);
    }
    
    public void setDOFs(){
        if(this.getSpaceIntegrator()==null){
            System.err.println("Null SpaceIntegrator for elements of  domain: "+this.getID()+" default SpaceIntegrator is automatically defined...");
            SpaceIntegrator SI = null;
            switch(theFundSol.getSpaceDimension()){
                case 1: SI=new SpaceNodeIntegrator(); break;
                case 2: SI=new SpaceLineIntegrator(6); break;
                case 3: SI=new SpaceQuadIntegrator(6); break;
            }
            this.setSpaceIntegrator(SI);
        }
        sumDOFS=0;
        uDOFS=0;
        pDOFS=0;
        int States_Combos=this.theStates.size()+this.theCombos.size(); if(States_Combos<1)States_Combos=1;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.init_p(this.theFundSol.get_p_DOFs(),1,States_Combos);
            if(this.auxiliaryFields)theNode.init_p_aux(this.theFundSol.get_p_DOFs(),1,States_Combos);
            sumDOFS+=theNode.getp().length;
            pDOFS+=theNode.getp().length;
            theNode.init_u(this.theFundSol.get_u_DOFs(),1,States_Combos);
            if(this.auxiliaryFields)theNode.init_u_aux(this.theFundSol.get_u_DOFs(),1,States_Combos);
            sumDOFS+=theNode.getu().length;
            uDOFS+=theNode.getu().length;
            theNode.init_v(this.theFundSol.get_v_DOFs(),1,States_Combos);
            sumDOFS+=theNode.getv().length;
        }
        int nce;
        int np=beginDOF+1; int pstep=this.theFundSol.get_p_DOFs();
        int nu=beginDOF+pDOFS+1; int ustep=this.theFundSol.get_u_DOFs();
        int nv=beginDOF+pDOFS+uDOFS+1; int vstep=this.theFundSol.get_v_DOFs();
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            nce=theNode.getNumOfConnectedElements();
            theNode.setpEFTable(np, pstep*nce);
            np+=pstep*nce;
            theNode.setuEFTable(nu, ustep);
            nu+=ustep;
            theNode.setvEFTable(nv, vstep);
            nv+=vstep;
            
        }
        //this.res = new Matrix(sumDOFS,1);
    }
    
    public void setDOFs(int steps){
        if(this.getSpaceIntegrator()==null){
            System.err.println("Null SpaceIntegrator for elements of  domain: "+this.getID()+" default SpaceIntegrator is automatically defined...");
            SpaceIntegrator SI = null;
            switch(theFundSol.getSpaceDimension()){
                case 1: SI=new SpaceNodeIntegrator(); break;
                case 2: SI=new SpaceLineIntegrator(6); break;
                case 3: SI=new SpaceQuadIntegrator(6); break;
            }
            this.setSpaceIntegrator(SI);
        }
        sumDOFS=0;
        uDOFS=0;
        pDOFS=0;
        int States_Combos=this.theStates.size()+this.theCombos.size(); if(States_Combos<1)States_Combos=1;
        // in this for loop it was steps instead of steps+1 untill 30/07/2016 when I changed it.
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.init_p(this.theFundSol.get_p_DOFs(),steps+1,States_Combos);
            if(this.auxiliaryFields)theNode.init_p_aux(this.theFundSol.get_p_DOFs(),steps+1,States_Combos);
            theNode.init_ph(this.theFundSol.get_p_DOFs(),steps+1);
            sumDOFS+=theNode.getp().length;
            pDOFS+=theNode.getp().length;
            theNode.init_u(this.theFundSol.get_u_DOFs(),steps+1,States_Combos);
            if(this.auxiliaryFields)theNode.init_u_aux(this.theFundSol.get_u_DOFs(),steps+1,States_Combos);
            theNode.init_uh(this.theFundSol.get_u_DOFs(),steps+1);
            sumDOFS+=theNode.getu().length;
            uDOFS+=theNode.getu().length;
            theNode.init_v(this.theFundSol.get_v_DOFs(),steps+1,States_Combos);
            sumDOFS+=theNode.getv().length;

        }
        int nce;
        int np=beginDOF+1; int pstep=this.theFundSol.get_p_DOFs();
        int nu=beginDOF+pDOFS+1; int ustep=this.theFundSol.get_u_DOFs();
        int nv=beginDOF+pDOFS+uDOFS+1; int vstep=this.theFundSol.get_v_DOFs();
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            nce=theNode.getNumOfConnectedElements();
            theNode.setpEFTable(np, pstep*nce);
            np+=pstep*nce;
            theNode.setuEFTable(nu, ustep);
            nu+=ustep;
            theNode.setvEFTable(nv, vstep);
            nv+=vstep;
            
        }
        this.res = new AbstractMatrix(sumDOFS,1);
        this.Lmat = new AbstractMatrix(this.getNumNodes(),this.getNumNodes());
    }

    public void setBeginDOF(int beg){
        this.beginDOF=beg;
    }

    public int getBeginDOFs(){
        return this.beginDOF;
    }
    
    public int getSumDOFs(){
        return (sumDOFS+RigidBodyRemoval);
    }
    
    public int getuDOFs(){
        return this.uDOFS;
    }
    
    public int getpDOFs(){
        return this.pDOFS;
    }
    
    public int getRigidRemovalDOFs(){
        return this.RigidBodyRemoval;
    }
    
    public int getNumberOfConstraintEquations(){
        return this.theConstraintEquations.size();
    }
    
    public Map getConstraintEquations(){
        return this.theConstraintEquations;
    }
    
    public FundamentalSolution getFundamentalSolution(){
        return this.theFundSol;
    }
    
    public void updateNodes(AbstractMatrix X, int step){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            double val;
            // update p
            for(int i=0; i<aNode.getp().length; i++){
                val=X.get(aNode.getpEFTable()[i]-1, 0);
                aNode.setp(i, step, val);
            }
            // update u
            for(int i=0; i<aNode.getu().length; i++){
                val=X.get(aNode.getuEFTable()[i]-1, 0);
                aNode.setu(i, step, val);
            }
             // update v
            for(int i=0; i<aNode.getv().length; i++){
                val=X.get(aNode.getvEFTable()[i]-1, 0);
                aNode.setv(i, step, val);
            }
        }
    }
    
    public void updateNodes(AbstractMatrix X, int step, int state){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            double val;
            // update p
            for(int i=0; i<aNode.getp().length; i++){
//                val=X.get(aNode.getpEFTable()[i]-1-this.beginDOF, 0);
                val=X.get(aNode.getpEFTable()[i]-1, 0);
                if(this.auxiliaryFields==false){aNode.setp(i, step, val, state);}
                else{aNode.setp_aux(i, step, val, state);}
            }
            // update u
            for(int i=0; i<aNode.getu().length; i++){
//                val=X.get(aNode.getuEFTable()[i]-1-this.beginDOF, 0);
                val=X.get(aNode.getuEFTable()[i]-1, 0);
                if(this.auxiliaryFields==false){aNode.setu(i, step, val, state);}
                else{aNode.setu_aux(i, step, val, state);}
            }
             // update v
            for(int i=0; i<aNode.getv().length; i++){
//                val=X.get(aNode.getvEFTable()[i]-1-this.beginDOF, 0);
                val=X.get(aNode.getvEFTable()[i]-1, 0);
                aNode.setv(i, step, val, state);
            }
        }
    }
    
    public void printNodes(int step){
        printNodes(step,0);
    }
    
    public void printNodes(int step, int state){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            System.out.println("Node id: "+aNode.getID());
            System.out.print("p: ");
            for(int i=0; i<aNode.getp().length; i++){
                System.out.print(aNode.getp()[i][step][state]);
                System.out.print(" ");
            }
            System.out.println();
            System.out.print("u: ");
            for(int i=0; i<aNode.getu().length; i++){
                System.out.print(aNode.getu()[i][step][state]);
                System.out.print(" ");
            }
            System.out.println();
            System.out.print("v: ");
            for(int i=0; i<aNode.getv().length; i++){
                System.out.print(aNode.getv()[i][step][state]);
                System.out.print(" ");
            }
            System.out.println();
        }
    }
    
    public void printElements(int step){
        printElements(step, 0);
    }
    
    public void printElements(int step, int state){
        int ndf=this.theFundSol.get_p_DOFs();
        for(Iterator<Element>it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement=it.next();
            System.out.println("Element with id= "+theElement.getID());
            for(Iterator<Node>nit=theElement.getNodes().values().iterator();nit.hasNext();){
                Node theNode=nit.next();
                System.out.println("Node id: "+theNode.getID());
                System.out.print("p : ");
                for(int i=0; i<theNode.getp(theElement, ndf).length; i++){
                    System.out.print(theNode.getp(theElement, ndf,state)[i][step]);
                    System.out.print(" ");
                }
                System.out.println();
            }
        }
    }
    
    public void printNodesGeom(){
        System.out.println("Nodes of domain with id= "+this.id);
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            aNode.print();
        }
    }
    
    public void printNodeDofs(){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            System.out.println("Node id: "+aNode.getID());
            System.out.print("p :");
            for(int i=0; i<aNode.getpEFTable().length; i++){
                System.out.print(aNode.getpEFTable()[i]);
                System.out.print(" ");
            }
            System.out.println();
            System.out.print("u : ");
            for(int i=0; i<aNode.getuEFTable().length; i++){
                System.out.print(aNode.getuEFTable()[i]);
                System.out.print(" ");
            }
            System.out.println();
            System.out.print("v :");
            for(int i=0; i<aNode.getvEFTable().length; i++){
                System.out.print(aNode.getvEFTable()[i]);
                System.out.print(" ");
            }
            System.out.println();
        }
    }
    
    public int getNumNodes(){
        return this.theNodes.size();
    }

    public int getNumElements(){
        return this.theElements.size();
    }
    
    public void setResponse(double coef, int fromwhichStep, int toWhichStep, int whichState){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(this.auxiliaryFields==false){aNode.setResponse(coef, fromwhichStep, toWhichStep, whichState);}
            else{aNode.setResponse_aux(coef, fromwhichStep, toWhichStep, whichState);}
        }
    }
    
    public void setZeroResponse(int toWhichStep, int whichState){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(this.auxiliaryFields==false){aNode.setZeroResponse(toWhichStep, whichState);}
            else{aNode.setZeroResponse_aux(toWhichStep, whichState);}
        }
    }
    
    public AbstractMatrix getResponse(int step){
        this.res = new AbstractMatrix(sumDOFS,1);
        double[][][] values;
        int[] eft;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            values=aNode.getp();
            eft=aNode.getpEFTable();
            for(int i=0; i<eft.length; i++){
                this.res.putVal(eft[i]-1-this.beginDOF, 0, values[i][step][0]);
            }
            values=aNode.getu();
            eft=aNode.getuEFTable();
            for(int i=0; i<eft.length; i++){
                this.res.putVal(eft[i]-1-this.beginDOF, 0, values[i][step][0]);
            }
            values=aNode.getv();
            eft=aNode.getvEFTable();
            for(int i=0; i<eft.length; i++){
                this.res.putVal(eft[i]-1-this.beginDOF, 0, values[i][step][0]);
            }
            
        }
        return this.res;
    }
    
    public AbstractMatrix getResponse(int step, int state){
        this.res = new AbstractMatrix(sumDOFS,1);
        double[][][] values;
        int[] eft;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(this.auxiliaryFields==false){values=aNode.getp();}else{values=aNode.getp_aux();}
            eft=aNode.getpEFTable();
            for(int i=0; i<eft.length; i++){
                this.res.putVal(eft[i]-1-this.beginDOF, 0, values[i][step][state]);
            }
            if(this.auxiliaryFields==false){values=aNode.getu();}else{values=aNode.getu_aux();}
            eft=aNode.getuEFTable();
            for(int i=0; i<eft.length; i++){
                this.res.putVal(eft[i]-1-this.beginDOF, 0, values[i][step][state]);
            }
            values=aNode.getv();
            eft=aNode.getvEFTable();
            for(int i=0; i<eft.length; i++){
                this.res.putVal(eft[i]-1-this.beginDOF, 0, values[i][step][state]);
            }
            
        }
        return this.res;
    }
    
    public void printNode(int step,int id_){
        Node aNode;
        aNode=this.theNodes.get(id_);
        for(int i=0; i<aNode.getp().length; i++){
            System.out.print(aNode.getp()[i][step]);
            System.out.print(",");
        }
        for(int i=0; i<aNode.getu().length; i++){
            System.out.print(aNode.getu()[i][step]);
            System.out.print(",");
        }
        for(int i=0; i<aNode.getv().length; i++){
            System.out.print(aNode.getv()[i][step]);
            if(i<aNode.getv().length-1){System.out.print(",");}
            else{System.out.println();}
        }
    }
    
    public void printNode(int step,int id_, int state){
        Node aNode;
        aNode=this.theNodes.get(id_);
        for(int i=0; i<aNode.getp().length; i++){
            System.out.print(aNode.getp()[i][step][state]);
            System.out.print(",");
        }
        for(int i=0; i<aNode.getu().length; i++){
            System.out.print(aNode.getu()[i][step][state]);
            System.out.print(",");
        }
        for(int i=0; i<aNode.getv().length; i++){
            System.out.print(aNode.getv()[i][step][state]);
            if(i<aNode.getv().length-1){System.out.print(",");}
            else{System.out.println();}
        }
    }
    
    public void printElement(int id_, int steps, String fileName){
        PrintWriter outFile = null;
        int pdof=this.theFundSol.get_p_DOFs();
        try {
            outFile = new PrintWriter(new FileWriter(fileName));
        } catch (IOException ex) {
            Logger.getLogger(Domain.class.getName()).log(Level.SEVERE, null, ex);
        }
        Element anElement=this.theElements.get(id_);
        System.out.print(anElement.id);System.out.print(" ");
        for(int step=0; step<=steps; step++){
            outFile.print(step);outFile.print(" ");
            for(Iterator<Node> it=anElement.getNodes().values().iterator(); it.hasNext();){
                Node elemNode=it.next();
                if(step==0){System.out.print(elemNode.getID());System.out.print(" ");}
                for(int i=0; i<elemNode.getu().length; i++){
                    outFile.print(elemNode.getu()[i][step]);outFile.print(" ");
                }
                for(int i=0; i<elemNode.getv().length; i++){
                    outFile.print(elemNode.getv()[i][step]);outFile.print(" ");
                }
                for(int i=0; i<elemNode.getp(anElement, pdof).length; i++){
                    outFile.print(elemNode.getp(anElement, pdof)[i][step]);outFile.print(" ");
                }
            }
            outFile.println();
        }
        System.out.println();
        outFile.close();
    }
    
    public void printElement(int id_, int steps, String fileName, double dt){
        PrintWriter outFile = null;
        int pdof=this.theFundSol.get_p_DOFs();
        try {
            outFile = new PrintWriter(new FileWriter(fileName));
        } catch (IOException ex) {
            Logger.getLogger(Domain.class.getName()).log(Level.SEVERE, null, ex);
        }
        Element anElement=this.theElements.get(id_);
        System.out.print(anElement.id);System.out.print(" ");
        for(int step=0; step<=steps; step++){
            outFile.print(step*dt);outFile.print(" ");
            for(Iterator<Node> it=anElement.getNodes().values().iterator(); it.hasNext();){
                Node elemNode=it.next();
                if(step==0){System.out.print(elemNode.getID());System.out.print(" ");}
                for(int i=0; i<elemNode.getu().length; i++){
                    outFile.print(elemNode.getu()[i][step][0]);outFile.print(" ");
                }
                for(int i=0; i<elemNode.getv().length; i++){
                    outFile.print(elemNode.getv()[i][step][0]);outFile.print(" ");
                }
                for(int i=0; i<elemNode.getp(anElement, pdof).length; i++){
                    outFile.print(elemNode.getp(anElement, pdof)[i][step]);outFile.print(" ");
                }
            }
            outFile.println();
        }
        System.out.println();
        outFile.close();
    }
    
    public AbstractMatrix getLmat(){
        this.Lmat.init();
        int i=0;
        int j=0;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            j=0;
            for(Iterator<Node>jt=this.theNodes.values().iterator();jt.hasNext();){
                Node bNode=jt.next();
                Lmat.set(i, j, aNode.getDist(bNode));
                j++;
            }
            i++;
        }
        return this.Lmat;
    }
    
    public void generateVtk(String fileName, int steps, int DeltaStep){
        // exei graftei eleafros proxeira
        // isxyei mono gia 4komva stoixeia
        // ENTELOS EXIDIKEYMENO
        PrintWriter outFile = null;
        System.out.println();
        try {
            outFile = new PrintWriter(new FileWriter(fileName));
        } catch (IOException ex) {
            Logger.getLogger(Domain.class.getName()).log(Level.SEVERE, null, ex);
        }
        int numN=this.theNodes.size();
        int[] map = new int[numN];
        outFile.println("# vtk DataFile Version 3.0");
        outFile.println("jbem vtk output");
        outFile.println("ASCII");
        outFile.println("DATASET UNSTRUCTURED_GRID");
        outFile.println();
        outFile.println("POINTS "+numN+" float");
        int count=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node aNode=it.next();
            map[count]=aNode.getID();
            outFile.println(aNode.getCoordinates()[0]+" "+
                            aNode.getCoordinates()[1]+" "+
                            aNode.getCoordinates()[2]);
            count+=1;
        }
        outFile.println();
        numN=this.theElements.size();
        int nn=numN*5;
        outFile.println("CELLS "+numN+" "+nn);
        int n1,n2,n3,n4;
        int m1 = 0,m2 = 0,m3 = 0,m4 = 0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element aElement=it.next();
            n1=aElement.getNodeHier(1).getID();
            n2=aElement.getNodeHier(2).getID();
            n3=aElement.getNodeHier(3).getID();
            n4=aElement.getNodeHier(4).getID();
            
            for(int i=0;i<map.length;i++){
                if(n1==map[i]){
                    m1=i;
                    break;
                }
            }
            
            for(int i=0;i<map.length;i++){
                if(n2==map[i]){
                    m2=i;
                    break;
                }
            }
            
            for(int i=0;i<map.length;i++){
                if(n3==map[i]){
                    m3=i;
                    break;
                }
            }
            
            for(int i=0;i<map.length;i++){
                if(n4==map[i]){
                    m4=i;
                    break;
                }
            }
            outFile.println(4+" "+m1+" "+m2+" "+m3+" "+m4);
        }
        outFile.println();
        outFile.println("CELL_TYPES"+" "+numN);
        for(int i=1;i<=numN;i++){
            outFile.println(9);
        }
        outFile.println();
        outFile.println("POINT_DATA"+" "+this.theNodes.size());
        outFile.println();
        double[][][] d;
        //System.out.println("x y z ux uy uz");
        for(int ii=0;ii<=steps;ii+=DeltaStep){
            outFile.println("VECTORS u"+ii+" float");
            for(int i=0;i<map.length;i++){
                //System.out.println(this.getNode(map[i]).getCoordinates()[0]+" "+this.getNode(map[i]).getCoordinates()[1]+" "+this.getNode(map[i]).getCoordinates()[2]+" "+
                //        this.getNode(map[i]).getu()[0][0]+" "+this.getNode(map[i]).getu()[1][0]+" "+this.getNode(map[i]).getu()[2][0]+" ");
                d=this.getNode(map[i]).getu();
                outFile.println(d[0][ii][0]+" "+d[1][ii][0]+" "+d[2][ii][0]);
            }
            outFile.println();
            outFile.println("VECTORS ux"+ii+" float");
            for(int i=0;i<map.length;i++){
                d=this.getNode(map[i]).getu();
                outFile.println(d[0][ii][0]+" "+0.+" "+0.);
            }
            outFile.println();
            outFile.println("VECTORS uy"+ii+" float");
            for(int i=0;i<map.length;i++){
                d=this.getNode(map[i]).getu();
                outFile.println(0.+" "+d[1][ii][0]+" "+0.);
            }
            outFile.println();
            outFile.println("VECTORS uz"+ii+" float");
            for(int i=0;i<map.length;i++){
                d=this.getNode(map[i]).getu();
                outFile.println(0.+" "+0.+" "+d[2][ii][0]);
            }
            outFile.println();
        }
        
        // velocities
        /*for(int ii=0;ii<steps;ii++){
            outFile.println("VECTORS v"+ii+" float");
            for(int i=0;i<map.length;i++){
                //System.out.println(this.getNode(map[i]).getCoordinates()[0]+" "+this.getNode(map[i]).getCoordinates()[1]+" "+this.getNode(map[i]).getCoordinates()[2]+" "+
                //        this.getNode(map[i]).getu()[0][0]+" "+this.getNode(map[i]).getu()[1][0]+" "+this.getNode(map[i]).getu()[2][0]+" ");
                d=this.getNode(map[i]).getv();
                outFile.println(d[0][ii]+" "+d[1][ii]+" "+d[2][ii]);
            }
            outFile.println();
            outFile.println("VECTORS vx"+ii+" float");
            for(int i=0;i<map.length;i++){
                d=this.getNode(map[i]).getv();
                outFile.println(d[0][ii]+" "+0.+" "+0.);
            }
            outFile.println();
            outFile.println("VECTORS vy"+ii+" float");
            for(int i=0;i<map.length;i++){
                d=this.getNode(map[i]).getv();
                outFile.println(0.+" "+d[1][ii]+" "+0.);
            }
            outFile.println();
            outFile.println("VECTORS vz"+ii+" float");
            for(int i=0;i<map.length;i++){
                d=this.getNode(map[i]).getv();
                outFile.println(0.+" "+0.+" "+d[2][ii]);
            }
            outFile.println();
        }*/
        
        outFile.close();
    }
    
    public void generateVtk_Geometry(String fileName){
        // exei graftei eleafros proxeira
        // isxyei mono gia 4komva stoixeia
        // ENTELOS EXIDIKEYMENO
        PrintWriter outFile = null;
        System.out.println();
        try {
            outFile = new PrintWriter(new FileWriter(fileName));
        } catch (IOException ex) {
            Logger.getLogger(Domain.class.getName()).log(Level.SEVERE, null, ex);
        }
        int numN=this.theNodes.size();
        int[] map = new int[numN];
        outFile.println("# vtk DataFile Version 3.0");
        outFile.println("jbem vtk output");
        outFile.println("ASCII");
        outFile.println("DATASET UNSTRUCTURED_GRID");
        outFile.println();
        outFile.println("POINTS "+numN+" float");
        int count=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node aNode=it.next();
            map[count]=aNode.getID();
            outFile.println(aNode.getCoordinates()[0]+" "+
                            aNode.getCoordinates()[1]+" "+
                            aNode.getCoordinates()[2]);
            count+=1;
        }
        outFile.println();
        numN=this.theElements.size();
        int nn=numN*5;
        outFile.println("CELLS "+numN+" "+nn);
        int n1,n2,n3,n4;
        int m1 = 0,m2 = 0,m3 = 0,m4 = 0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element aElement=it.next();
            n1=aElement.getNodeHier(1).getID();
            n2=aElement.getNodeHier(2).getID();
            n3=aElement.getNodeHier(3).getID();
            n4=aElement.getNodeHier(4).getID();
            
            for(int i=0;i<map.length;i++){
                if(n1==map[i]){
                    m1=i;
                    break;
                }
            }
            
            for(int i=0;i<map.length;i++){
                if(n2==map[i]){
                    m2=i;
                    break;
                }
            }
            
            for(int i=0;i<map.length;i++){
                if(n3==map[i]){
                    m3=i;
                    break;
                }
            }
            
            for(int i=0;i<map.length;i++){
                if(n4==map[i]){
                    m4=i;
                    break;
                }
            }
            outFile.println(4+" "+m1+" "+m2+" "+m3+" "+m4);
        }
        outFile.println();
        outFile.println("CELL_TYPES"+" "+numN);
        for(int i=1;i<=numN;i++){
            outFile.println(9);
        }
        outFile.println();
        
        outFile.close();
    }
    
    public double getMinL(){
        double minL=0.;
        double anL=0.;
        int indicator=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            EPlane theElement = (EPlane) it.next();
            anL=theElement.getMinL();
            if(indicator==0){minL=anL;}
            if(anL<minL){minL=anL;}
            indicator+=1;
        }
        return minL;
    }
    
    public double getMaxL(){
        double minL=0.;
        double anL=0.;
        int indicator=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            EPlane theElement = (EPlane) it.next();
            anL=theElement.getMaxL();
            if(indicator==0){minL=anL;}
            if(anL>minL){minL=anL;}
            indicator+=1;
        }
        return minL;
    }

    public void printElementsConectivity(){
        int nnode;
        System.out.println("Domain id="+this.id+", has "+this.theElements.size()+" elements with id and nodes:");
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            System.out.print(theElement.getID()+" ");
            nnode=theElement.getNumNodes();
            for(int i=1;i<=nnode;i++){System.out.print(theElement.getNodeHier(i).getID()+" ");}
            System.out.println(theElement.getBCType());
        }
    }

    public void printElementsConectivity(javax.swing.JTextArea TextArea){
        int nnode;
        TextArea.append("Domain id="+this.id+", has "+this.theElements.size()+" elements with id and nodes:"+'\n');
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            ELine theElement = (ELine) it.next();
            //TextArea.append(theElement.getID()+" "+" of Length="+theElement.getLength()+'\n');
            nnode=theElement.getNumNodes();
            for(int i=1;i<=nnode;i++){
                TextArea.append(theElement.getNodeHier(i).getID()+" ");
            }
            TextArea.append(""+'\n');
        }
    }

    public double getMinX(){
        double mx=0.;
        double x1;
        int c=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            ++c;
            Node Node1 = it.next();
            x1=Node1.getCoordinates()[0];
            if(c==1){mx=x1;}
            if(mx>x1){mx=x1;}
        }
        return mx;
    }

    public double getMaxX(){
        double mx=0.;
        double x1;
        int c=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            ++c;
            Node Node1 = it.next();
            x1=Node1.getCoordinates()[0];
            if(c==1){mx=x1;}
            if(mx<x1){mx=x1;}
        }
        return mx;
    }
    
    public double getMinZ(){
        double mx=0.;
        double z1;
        int c=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            ++c;
            Node Node1 = it.next();
            z1=Node1.getCoordinates()[2];
            if(c==1){mx=z1;}
            if(mx>z1){mx=z1;}
        }
        return mx;
    }

    public double getMaxZ(){
        double mx=0.;
        double z1;
        int c=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            ++c;
            Node Node1 = it.next();
            z1=Node1.getCoordinates()[2];
            if(c==1){mx=z1;}
            if(mx<z1){mx=z1;}
        }
        return mx;
    }

    public double getMinY(){
        double mx=0.;
        double x1;
        int c=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            ++c;
            Node Node1 = it.next();
            x1=Node1.getCoordinates()[1];
            if(c==1){mx=x1;}
            if(mx>x1){mx=x1;}
        }
        return mx;
    }

    public double getMaxY(){
        double mx=0.;
        double x1;
        int c=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            ++c;
            Node Node1 = it.next();
            x1=Node1.getCoordinates()[1];
            if(c==1){mx=x1;}
            if(mx<x1){mx=x1;}
        }
        return mx;
    }

    public double getWork(int wstep){
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            work+=elem.getWork(this,wstep);
        }
        return work;
    }
    
    public double getWork(int wstep, int state){
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            work+=elem.getWork(this,wstep,state);
        }
        return work;
    }
    
    public double getTotalPotentialEnergy(int wstep){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork(this,wstep);
            if(!elem.isNeumann()){work+=d;}else{work-=d;}
        }
        return work;
    }
    
    public double getTotalPotentialEnergy(int wstep, int state){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork(this,wstep,state);
            if(!elem.isNeumann()){work+=d;}else{work-=d;}
        }
        return work;
    }
    
    public double getPower(int wstep){
        double power=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            power+=elem.getPower(this,wstep);
        }
        return power;
    }
    
    public double getTIPower(int wstep){
        double power=0.;
        double dt = ((TimeFundamentalSolution) this.theFundSol).getDT();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            power+=elem.getTIPower(this,wstep,dt);
        }
        return power;
    }
    
    public double getWork_hmg(){
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            work+=elem.getWork_hmg(this);
        }
        return work;
    }
    
    
    public double getWork_Gamma_C(int wstep){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork(this,wstep);
            if(elem.isContact()){work+=d;}
        }
        return work;
    }
    
    public double getWork_Gamma_C(int wstep, int state){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork(this,wstep,state);
            if(elem.isContact()){work+=d;}
        }
        return work;
    }
    
    public double getWork_Gamma_D(int wstep, int state){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork(this,wstep,state);
            if(elem.isDirichlet()){work+=d;}
        }
        return work;
    }
    
    public double getWork_Gamma_N(int wstep, int state){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork(this,wstep,state);
            if(elem.isNeumann()){work+=d;}
        }
        return work;
    }
    
    public double getTotalPotentialEnergy_hmg(){
        double work=0.;
        double d;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            d=elem.getWork_hmg(this);
            if(!elem.isNeumann()){work+=d;}else{work+=d-2*d;}
        }
        return work;
    }

    public double getGeneralisedForce(int wstep, int nodeID, int dof){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        double d;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            d=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof);
            if(!this.theElements.get(elemIDS[i]).isNeumann()){work+=d;}else{work+=d-2*d;}
        }
        return work;
    }
    
    public double getGeneralisedForce(int wstep, int nodeID, int dof, int state){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        double d;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            d=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof, state);
            if(!this.theElements.get(elemIDS[i]).isNeumann()){work+=d;}else{work+=d-2*d;}
        }
        return work;
    }
    
    public double getGeneralisedForce_Gamma_C(int wstep, int nodeID, int dof){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            if(this.theElements.get(elemIDS[i]).isContact())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof);
        }
        return work;
    }
    
    public double getGeneralisedForce_Gamma_C(int wstep, int nodeID, int dof, int state){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            if(this.theElements.get(elemIDS[i]).isContact())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof, state);
        }
        return work;
    }
    
    public double getGeneralisedForce_Gamma_D(int wstep, int nodeID, int dof, int state){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            if(this.theElements.get(elemIDS[i]).isDirichlet())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof, state);
        }
        return work;
    }
    
    public double getGeneralisedForce_Gamma_N(int wstep, int nodeID, int dof, int state){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            if(this.theElements.get(elemIDS[i]).isNeumann())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof, state);
        }
        return work;
    }
    
    public double getGeneralisedForce_hmg(int nodeID, int dof){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        double d;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            d=this.theElements.get(elemIDS[i]).getGeneralisedForce_hmg(this, nodeID, dof);
            if(!this.theElements.get(elemIDS[i]).isNeumann()){work+=d;}else{work+=d-2*d;}
            //work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, wstep, nodeID, dof);
        }
        return work;
    }


    public double getSideInequality(int step, int bound){
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            //if(elem.isExternalLoading())work+=elem.getSideInequality(this, step, bound);
            work+=elem.getSideInequality(this, step, bound);
        }
        return work;
    }
    
    public void setDissipativeSurface(){
        this.isDissipativeSurface=true;
    }

    public boolean isDissipativeSurface(){
        return this.isDissipativeSurface;
    }

    public Material getMaterial(){
        return this.theMaterial;
    }

    public void clearConstraints(){
        theConstraintEquations.clear();
    }
    
    public void clearInternalPoints(){
        this.theResultPoints.clear();
    }
    
    public void clear(){
        theResultPoints.clear();
        theNodes.clear();
        theElements.clear();
        theConstraintEquations.clear();
        theConnectedNodes.clear();
        theMaterial=null;
        theFundSol=null;
        sumDOFS=0;
        uDOFS=0;
        pDOFS=0;
        extension=Extension.FINITE;
        res=null;
        Lmat=null;
        beginDOF=0;
        isDissipativeSurface=false;
        eigenvalues=null;
        eigenvectors=null;
        RigidBodyRemoval=0;
        rigidMultipliers=null;
        fixeddofs=null;
        fixedmodes=null;
        theCombos.clear();
        theStates.clear();
        responseID=0;
        auxiliaryFields=false;
        currentFieldAuxiliary=false;
        tol=1.e-12;
        
        uniformTempChange=0.;
        maximum_X_coordinate=Double.NEGATIVE_INFINITY;
        minimum_X_coordinate=Double.POSITIVE_INFINITY;
        maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
        minimum_Y_coordinate=Double.POSITIVE_INFINITY;
        maximum_Z_coordinate=Double.NEGATIVE_INFINITY;
        minimum_Z_coordinate=Double.POSITIVE_INFINITY;
    }

    public void updateNodesT(AbstractMatrix X){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            double val;
            if(aNode.getuh()!=null && aNode.getph()!=null){
                // update p
                for(int i=0; i<aNode.getp().length; i++){
                    val=X.get(aNode.getpEFTable()[i]-1, 0);
                    aNode.setph(i, val);
                }
                // update u
                for(int i=0; i<aNode.getu().length; i++){
                    val=X.get(aNode.getuEFTable()[i]-1, 0);
                    aNode.setuh(i, val);
                }
            }
        }
    }
    
    public void PrintConstraintEquations(){
         System.out.println("---- Constraint Equations Of Domain with id = "+this.id+", number of C.E.s = "+theConstraintEquations.size()+" -----");
                 System.out.println();
                 for(Iterator<ConstraintEquation>
                         it=this.theConstraintEquations.values().iterator(); it.hasNext();){
                     ConstraintEquation theConstraintEquation = it.next();
                     theConstraintEquation.Print();
                 }

        System.out.println("-------------------------------------------------------------");
    }

    public void setEIGENSYSTEM(double[] eigenvalues, AbstractMatrix eigenvectors){
        this.eigenvalues=eigenvalues;
        this.eigenvectors=eigenvectors;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            aNode.init_u_eig(eigenvalues.length);
            double val;
            for(int j=0;j<eigenvalues.length;j++){
                for(int i=0; i<aNode.getu().length; i++){
                    val=eigenvectors.get(aNode.getuEFTable()[i]-1, j);
                    aNode.setu_eig(i,j, val);
                }
            }
            // update u_eig

        }
    }

    public double[] getEigenvalues(){
        return this.eigenvalues;
    }

    public double getArea(){
        double A=0.;
        double x1,y1,x2,y2;
        int n;
        Element theElement;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            theElement = it.next();
            n=theElement.getNumNodes();
            for(int i=1;i<n;i++){
                x1=theElement.getNodeHier(i).getCoordinates()[0];
                y1=theElement.getNodeHier(i).getCoordinates()[1];
                x2=theElement.getNodeHier(i+1).getCoordinates()[0];
                y2=theElement.getNodeHier(i+1).getCoordinates()[1];
                A+=(x1*y2-x2*y1)/2.;
            }
        }
        return A;
    }

    public double getPerimeter(){
        double A=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            ELine2 theElement = (ELine2) it.next();
            A+=theElement.getLength();
        }
        return A;
    }

    public double getArea_eig(double h, int eig){
        double A=0.;
        double x1,y1,x2,y2;
        int n;
        Element theElement;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            theElement = it.next();
            n=theElement.getNumNodes();
            for(int i=1;i<n;i++){
                x1=theElement.getNodeHier(i).getCoordinates()[0]+h*theElement.getNodeHier(i).getu_eig()[0][eig];
                y1=theElement.getNodeHier(i).getCoordinates()[1]+h*theElement.getNodeHier(i).getu_eig()[1][eig];
                x2=theElement.getNodeHier(i+1).getCoordinates()[0]+h*theElement.getNodeHier(i+1).getu_eig()[0][eig];
                y2=theElement.getNodeHier(i+1).getCoordinates()[1]+h*theElement.getNodeHier(i+1).getu_eig()[1][eig];
                A+=(x1*y2-x2*y1)/2.;
            }
        }
        return A;
    }

    public double getPerimeter_eig(double h, int eig){
        double A=0.;
        double x1,y1,x2,y2;
        int n;
        Element theElement;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            theElement = it.next();
            n=theElement.getNumNodes();
            for(int i=1;i<n;i++){
                x1=theElement.getNodeHier(i).getCoordinates()[0]+h*theElement.getNodeHier(i).getu_eig()[0][eig];
                y1=theElement.getNodeHier(i).getCoordinates()[1]+h*theElement.getNodeHier(i).getu_eig()[1][eig];
                x2=theElement.getNodeHier(i+1).getCoordinates()[0]+h*theElement.getNodeHier(i+1).getu_eig()[0][eig];
                y2=theElement.getNodeHier(i+1).getCoordinates()[1]+h*theElement.getNodeHier(i+1).getu_eig()[1][eig];
                A+=Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
            }
        }
        return A;
    }

    public int getNumEigs(){
        int n=0;
        if(this.eigenvalues!=null)n=eigenvalues.length;
        return n;
    }

    public boolean ExistNodeWithID(int theID){
        return this.theNodes.containsKey(theID);
    }
    
    public boolean ExistNodeWithCoords(double[] coords){
        boolean answer=false;
        double dist;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            dist=0;
            for(int i=0;i<coords.length;i++){
                dist+=(coords[i]-aNode.getCoordinates()[i])*(coords[i]-aNode.getCoordinates()[i]);
            }
            if(Math.sqrt(dist)<=1.e-14)answer=true;
        }
        return answer;
    }

    public boolean ExistElementWithID(int theID){
        return this.theElements.containsKey(theID);
    }
    
    public void RemoveRigidBody(int steps){
        this.RigidBodyRemoval=3*(this.theFundSol.get_u_DOFs()-1);
        if(RigidBodyRemoval==3){
            this.rigidMultipliers = new double[RigidBodyRemoval][steps];
            this.fixedmodes = new int[RigidBodyRemoval];
            fixedmodes[0]=0; // 0 for horizontal motion
            fixedmodes[1]=1; // 1 for vertical motion
            fixedmodes[2]=2; // 2 for rotational motion
        }else{
            this.rigidMultipliers = new double[RigidBodyRemoval][steps];
            this.fixedmodes = new int[RigidBodyRemoval];
            fixedmodes[0]=0; // 0 for x motion
            fixedmodes[1]=1; // 1 for y motion
            fixedmodes[2]=2; // 2 for xy or z rotation
            fixedmodes[3]=3; // 3 for z motion
            fixedmodes[4]=4; // 4 for zx or y rotation
            fixedmodes[5]=5; // 5 for yz or x rotation
        }
    }
    
    public void FixCertainDofs(int[] dofs,int steps){
        this.RigidBodyRemoval=dofs.length;
        this.fixeddofs = new int[dofs.length];
        this.rigidMultipliers = new double[dofs.length][steps];
        System.arraycopy(dofs, 0, fixeddofs, 0, dofs.length);
        fixM=false;
    }
    
    public void FixCertainModes(int[] dofs,int steps){
        this.RigidBodyRemoval=dofs.length;
        this.fixedmodes = new int[dofs.length];
        this.rigidMultipliers = new double[dofs.length][steps];
        System.arraycopy(dofs, 0, fixedmodes, 0, dofs.length);
    }
    
    public void updateRigidMultipliers(AbstractMatrix X, int step){
        for(int i=0;i<this.RigidBodyRemoval;i++){
            this.rigidMultipliers[i][step] = X.get(X.getRowDimension()-RigidBodyRemoval+i, 0);
        }
    }
    
    public double[][] getRigidMultipliers(){return this.rigidMultipliers;}
            
    public int[] getFixedDofs(){return this.fixeddofs;}
    
    public int[] getFixedModes(){return this.fixedmodes;}
    
    public int getNumOfConstraintsOnNode(int nodeid){
        int n=0;
        boolean exist=false;
        for(Iterator<ConstraintEquation>it=this.theConstraintEquations.values().iterator(); it.hasNext();){
                     ConstraintEquation theConstraintEquation = it.next();
                     exist=false;
                     for(Iterator<ConstraintTerm> nt=theConstraintEquation.getNodeConstraintTerm().iterator(); nt.hasNext();){
                         ConstraintTerm theNCT = nt.next();
                         if(theNCT.getNodeID()==nodeid)exist=true;
                     }
                     
                     for(Iterator<ConstraintTermElement> nt=theConstraintEquation.getElementConstraintTerm().iterator(); nt.hasNext();){
                         ConstraintTermElement theNCT = nt.next();
                         if(theNCT.getNodeID()==nodeid)exist=true;
                     }
                     if(exist)n+=1;
                 }
        return n;
    }
    
    public void checkConstraints(){
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
            System.out.println("For node : "+theNode.getID()+" c equations :"+this.getNumOfConstraintsOnNode(theNode.getID()));
        }
    }
    
    public void setForWMatrixRBM(){
        this.fixedmodes = new int[3];
        fixedmodes[0]=0; // 0 for horizontal motion
        fixedmodes[1]=1; // 1 for vertical motion
        fixedmodes[2]=2; // 2 for rotational motion
        this.fixM=true;
    }
    
    public void setForWMatrixRBM(int[] dofs){
        this.fixedmodes = new int[dofs.length];
        System.arraycopy(dofs, 0, fixedmodes, 0, dofs.length);
        fixedmodes[0]=0; // 0 for horizontal motion
        fixedmodes[1]=1; // 1 for vertical motion
        fixedmodes[2]=2; // 2 for rotational motion
        this.fixM=true;
    }
    
    public boolean getForWMatrixRBM(){return this.fixM;}
    
    public void calculateResponseCombos(int step){
        for(Iterator<Combo> ct=this.getCombos().values().iterator(); ct.hasNext();){
            Combo theCombo = ct.next();
            State theState;
            AbstractMatrix X = new AbstractMatrix(this.getResponse(step).getRowDimension(),this.getResponse(step).getColumnDimension()) ;
            X.init(0.0);
            for(int i=0;i<theCombo.getComboStates().length;i++){
                theState=this.theStates.get(theCombo.getComboStates()[i]);
                X=X.plus(this.getResponse(step, theState.getResponseID()).times(theCombo.getComboCoefs()[i]));
                //X.print(22, 22);
            }
            this.updateNodes(X, step,theCombo.getResponseID());
        }
    }

    public void Aux2MainVariables(int itime, double tau) {
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
            theNode.Aux2MainVariables(itime, tau, this.theMaterial);
        }
    }
    
    public void Aux2MainVariables_(int itime, double tau) {
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
            theNode.Aux2MainVariables_(itime, tau, this.theMaterial);
        }
    }

    public void printInternalPointsGeom() {
        System.out.println("Result Points of Domain with id: "+this.id);
        for(Iterator<ResultPoint>it=this.theResultPoints.values().iterator();it.hasNext();){
            ResultPoint aNode=it.next();
            System.out.print(aNode.getID()+" ");
            for(int i=0; i<aNode.getCoordinates().length; i++){
                System.out.print(aNode.getCoordinates()[i]+" ");
            }
            System.out.println(" "+aNode.isOnBoundary());
        }
    }

    public double getTIPower_TYPE(int itime, BoundaryType theBCType) {
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            switch (theBCType) {
                case DIRICHLET:
                    if(elem.isDirichlet()){work+=elem.getTIPower(this,itime);}
                    break;

                case NEUMANN:
                    if(elem.isNeumann()){work+=elem.getTIPower(this,itime);}
                    break;

                case CONTACT:
                    if(elem.isContact()){work+=elem.getTIPower(this,itime);}
                    break;

                default:
                    work+=elem.getTIPower(this,itime);
                    break;
            }
        }
        return work;
    }
    
    public double getTIPowerDifference_TYPE(int itime, BoundaryType theBCType) {
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            switch (theBCType) {
                case DIRICHLET:
                    if(elem.isDirichlet()){work+=elem.getTIPowerDifference(this,itime);}
                    break;

                case NEUMANN:
                    if(elem.isNeumann()){work+=elem.getTIPowerDifference(this,itime);}
                    break;

                case CONTACT:
                    if(elem.isContact()){work+=elem.getTIPowerDifference(this,itime);}
                    break;

                default:
                    work+=elem.getTIPowerDifference(this,itime);
                    break;
            }
        }
        return work;
    }
    
    public enum BoundaryType {
        DIRICHLET,NEUMANN,CONTACT,UNDEFINED
    }
    
    public enum Extension {
        FINITE,INFINITE,SEMI_INFINITE
    }
    
    public double getWork_TYPE(int DispStep, int TracStep, int DispState, int TracState, BoundaryType theBCType){
        double work=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            switch (theBCType) {
                case DIRICHLET:
                    if(elem.isDirichlet()){work+=elem.getWork(this,DispStep,TracStep,DispState,TracState);}
                    break;

                case NEUMANN:
                    if(elem.isNeumann()){work+=elem.getWork(this,DispStep,TracStep,DispState,TracState);}
                    break;

                case CONTACT:
                    if(elem.isContact()){work+=elem.getWork(this,DispStep,TracStep,DispState,TracState);}
                    break;

                default:
                    work+=elem.getWork(this,DispStep,TracStep,DispState,TracState);
                    break;
            }
        }
        return work;
    }
    
    public double getGeneralisedForce_TYPE(int DispStep, int TracStep, int DispState, int TracState, BoundaryType theBCType, int nodeID, int dof){
        if(dof!=1){if(dof!=2){System.err.println("generalized force for dof: "+dof); }}
        double work=0.;
        int[] elemIDS;
        elemIDS = theNodes.get(nodeID).getConnectedElementsIds();
        switch (theBCType) {
            case DIRICHLET:
                for(int i=0;i<elemIDS.length;i++){
                    if(this.theElements.get(elemIDS[i]).isDirichlet())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, DispStep,TracStep,DispState,TracState, nodeID, dof);
                }
                break;

            case NEUMANN:
                for(int i=0;i<elemIDS.length;i++){
                    if(this.theElements.get(elemIDS[i]).isNeumann())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, DispStep,TracStep,DispState,TracState, nodeID, dof);
                }
                break;

            case CONTACT:
                for(int i=0;i<elemIDS.length;i++){
                    if(this.theElements.get(elemIDS[i]).isContact())work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, DispStep,TracStep,DispState,TracState, nodeID, dof);
                }
                break;

            default:
                for(int i=0;i<elemIDS.length;i++){
                    work+=this.theElements.get(elemIDS[i]).getGeneralisedForce(this, DispStep,TracStep,DispState,TracState, nodeID, dof);
                }
                break;
        }
        
        return work;
    }
    
    public void setAuxiliaryVariables(){
        this.auxiliaryFields=true;
    }
    
    public double getInternalPointDisp(int IPid){
        return getInternalPointDisp(IPid, 0, 0, 0);
    }
    
    public double getInternalPointDisp(int IPid, int step, int wstate){
        return getInternalPointDisp(IPid, 0, step, wstate);
    }
    
    public double getInternalPointDisp(int IPid, int wdisp, int step, int wstate){
        double val=0.;
        if(this.theResultPoints.get(IPid).isOnBoundary()){
            val=getInternalPointDispOnBoundary(IPid, wdisp, step, wstate);
        }else{
            val=getInternalPointDispInterior(IPid, wdisp, step, wstate);
        }
        return val;
    }
    
    public double getInternalPointDispNorm(int IPid){
        return getInternalPointDispNorm(IPid, 0, 0);
    }
    
    public double getInternalPointDispNorm(int IPid, int step, int wstate){
        double val=0.;
        if(this.theResultPoints.get(IPid).isOnBoundary()){
            val=getInternalPointDispOnBoundary(IPid, step, wstate);
        }else{
            val=getInternalPointDispInterior(IPid, step, wstate);
        }
        return val;
    }
    
    public double getInternalPointDispInterior(int IPid, int step, int wstate){
        double val=0.;
        int ndofs=this.getFundamentalSolution().get_u_DOFs();
        double[] vals = new double[ndofs];
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            FundamentalSolution.theFSdata.setMaterial(theMaterial);
            for(int n=0;n<ndofs;n++){
                vals[n]+=elem.getInternalPointDisp(this, theResultPoints.get(IPid),n,step,wstate);
            }
        }
        for(int n=0;n<ndofs;n++){val+=vals[n]*vals[n];}
        return Math.sqrt(val);
    }
    
    public double getInternalPointDispInterior(int IPid, int wdisp, int step, int wstate){
        double val=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            FundamentalSolution.theFSdata.setMaterial(theMaterial);
            val+=elem.getInternalPointDisp(this, theResultPoints.get(IPid),wdisp,step,wstate);
        }
        return val;
    }
    
    public double getInternalPointDispOnBoundary(int IPid, int step, int wstate){
        // returns the norm
        double val=0.;
        int ndofs=this.getFundamentalSolution().get_u_DOFs();
        double[] vals = new double[ndofs];
        int onElemID = 0;
        double dist,cor1,cor2;
        int numOfElements=0;
        Element elem;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
            
            cor1=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
            if(cor1>=-1. && cor1<=1.&& dist<=tol){
                numOfElements+=1;
                onElemID=elem.getID();
            }
        }
        
        switch(numOfElements){
            case 0:
                System.err.println("No element has been found on boundary for InternalPoint with id: "+IPid+", coordinates of point x= "+this.theResultPoints.get(IPid).getCoordinates()[0]+", y= "+this.theResultPoints.get(IPid).getCoordinates()[1]);
                val=0.;
                break;
            case 1:
                // It is on element
                if(this.getFundamentalSolution().ndofs==2){
                    elem=((ELine) this.theElements.get(onElemID));
                    for(int i=1;i<=elem.getNumNodes();i++){
                        for(int n=0;n<ndofs;n++){
                            vals[n]+=((ELine) elem).getShapeFunction(i, elem.getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0])*elem.getNodeHier(i).getu()[n][step][wstate];
                        }
                    }
                }else{
                    elem=((EPlane) this.theElements.get(onElemID));
                    for(int i=1;i<=elem.getNumNodes();i++){
                        for(int n=0;n<ndofs;n++){
                            vals[n]+=((EPlane) elem).getShapeFunction(i, elem.getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0],elem.getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[1])*elem.getNodeHier(i).getu()[n][step][wstate];
                        }
                    }
                }
                
                break;
            default:
                // It coinsides with node
                int NodeID = 0;
                for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
                    Node theNode = it.next();
                    dist=theNode.getDist(this.theResultPoints.get(IPid));
                    if(dist<=tol){
                        NodeID=theNode.getID();
                    }
                }
                for(int n=0;n<ndofs;n++){
                    vals[n]=this.theNodes.get(NodeID).getu()[n][step][wstate];
                }
                break;
        }
        for(int n=0;n<ndofs;n++){val+=vals[n]*vals[n];}
        return Math.sqrt(val);
    }
    
    public double getInternalPointDispOnBoundary(int IPid, int wdisp, int step, int wstate){
        double val=0.;
        int onElemID = 0;
        double dist,cor1,cor2;
        int numOfElements=0;
        Element elem;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
            
            cor1=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
            if(cor1>=-1. && cor1<=1.&& dist<=tol){
                numOfElements+=1;
                onElemID=elem.getID();
            }
        }
        
        switch(numOfElements){
            case 0:
                System.err.println("No element has been found on boundary for InternalPoint with id: "+IPid+", coordinates of point x= "+this.theResultPoints.get(IPid).getCoordinates()[0]+", y= "+this.theResultPoints.get(IPid).getCoordinates()[1]);
                val=0.;
                break;
            case 1:
                // It is on element
                if(this.getFundamentalSolution().ndofs==2){
                    elem=((ELine) this.theElements.get(onElemID));
                    for(int i=1;i<=elem.getNumNodes();i++){
                        val+=((ELine) elem).getShapeFunction(i, elem.getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0])*elem.getNodeHier(i).getu()[wdisp][step][wstate];
                    }
                }else{
                    elem=((EPlane) this.theElements.get(onElemID));
                    for(int i=1;i<=elem.getNumNodes();i++){
                        val+=((EPlane) elem).getShapeFunction(i, elem.getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0],elem.getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[1])*elem.getNodeHier(i).getu()[wdisp][step][wstate];
                    }
                }
                
                break;
            default:
                // It coinsides with node
                int NodeID = 0;
                for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
                    Node theNode = it.next();
                    dist=theNode.getDist(this.theResultPoints.get(IPid));
                    if(dist<=tol){
                        NodeID=theNode.getID();
                    }
                }
                val=this.theNodes.get(NodeID).getu()[wdisp][step][wstate];
                break;
        }
        return val;
    }
    
    public double getInternalPointStress(int IPid, int wstress, int step, int wstate){
        // wstress  2D: 0 σxx, 1 σyy, 2 τxy
        //          3D: 0 σxx, 1 σyy, 2 σzz, 3 τxy, 4 τyz, 5 τxz
        // wstress 2D/3D: 6 VonMises
        double val=0.;
        if(this.theResultPoints.get(IPid).isOnBoundary()){
            val=getInternalPointStressOnBoundary(IPid, wstress, step, wstate);
        }else{
            val=getInternalPointStressInterior(IPid, wstress, step, wstate);
        }
        return val;
    }
    
    public double getInternalPointStressInterior(int IPid, int wstress, int step, int wstate){
        // wstress  2D: 0 σxx, 1 σyy, 2 τxy
        //          3D: 0 σxx, 1 σyy, 2 σzz, 3 τxy, 4 τyz, 5 τxz
        //              6 VonMises
        double val=0.;
        double sxx=0.0,syy=0.0,sxy=0.0,szz=0.0;
        if(wstress==6){
            if(this.theFundSol.get_u_DOFs()==2)for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
                Element elem = it.next();
                sxx+=elem.getInternalPointStress(this, theResultPoints.get(IPid),0,step, wstate);
                syy+=elem.getInternalPointStress(this, theResultPoints.get(IPid),1,step, wstate);
                sxy+=elem.getInternalPointStress(this, theResultPoints.get(IPid),2,step, wstate);
            }
            szz=((ElasticMat)this.getMaterial()).getLame_L()*(sxx+syy);
            sxx-=((ElasticMat)this.theMaterial).getExtendedThermalCoef()*this.uniformTempChange;
            syy-=((ElasticMat)this.theMaterial).getExtendedThermalCoef()*this.uniformTempChange;
            szz-=((ElasticMat)this.theMaterial).getExtendedThermalCoef()*this.uniformTempChange;
            if(((StaticElasticity2DFS) this.theFundSol).getPlaneStress())szz=0.0;
            val=Math.sqrt(0.5*(
                (sxx-syy)*(sxx-syy)
                +(syy-szz)*(syy-szz)
                +(szz-sxx)*(szz-sxx)
                +6.0*sxy*sxy)
        );
        }else{
            for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
                Element elem = it.next();
                val+=elem.getInternalPointStress(this, theResultPoints.get(IPid),wstress,step, wstate);
            }
            switch(this.theFundSol.get_u_DOFs()){
                case 2:
                    if(wstress==0 || wstress==1){
                        val-=((ElasticMat)this.theMaterial).getExtendedThermalCoef()*this.uniformTempChange;
                    }
                    break;
                case 3:
                    if(wstress==0 || wstress==1 || wstress==2){
                        val-=((ElasticMat)this.theMaterial).getExtendedThermalCoef()*this.uniformTempChange;
                    }
                    break;
                default:
                    System.err.println("error in method getInternalPointStress of class Domain");
                    System.exit(0);
            }
        }
        return val;
    }
    
    public double getInternalPointStressAux(int IPid, int wstress, int step, int wstate){
        // wstress  2D: 0 σxx, 1 σyy, 2 τxy
        double val=0.;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            val+=elem.getInternalPointStressAux(this, theResultPoints.get(IPid),wstress,step, wstate);
        }
        return val;
    }
    
    public double getInternalPointStressOnBoundary(int IPid, int wstress, int step, int wstate){
        // wstress  2D: 0 σxx, 1 σyy, 2 τxy
        //          3D: 0 σxx, 1 σyy, 2 σzz, 3 τxy, 4 τyz, 5 τxz
        //              6 VonMises
        double val=0.;
        double sxx=0.0, syy=0.0, sxy=0.0;
        int[] onElemID;
        double dist,cor;
        int numOfElements=0;
        Element elem;
        double eyl=0., tyl=0., txl=0.;
        double sxx_l=0.,syy_l=0.,sxy_l=0.;
        double s,c;
        double v=((ElasticMat) this.theMaterial).getPoissonRatio();
        //double G=((ElasticMat) this.theMaterial).getLame_M();
        double E=((ElasticMat) this.theMaterial).getElasticModulus();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            if(this.getFundamentalSolution().ndofs==2){
                dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
                cor=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
                if(cor>=-1. && cor<=1. && dist<=tol){
                    numOfElements+=1;
                }
            }else{
                if(elem.PointOnElement(theResultPoints.get(IPid))){
                    numOfElements+=1;
                }
            }
        }
        onElemID = new int[numOfElements];
        numOfElements=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            if(this.getFundamentalSolution().ndofs==2){
                dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
                cor=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
                if(cor>=-1. && cor<=1. && dist<=tol){
                    onElemID[numOfElements]=elem.getID();
                    numOfElements+=1;
                }
            }else{
                if(elem.PointOnElement(theResultPoints.get(IPid))){
                    numOfElements+=1;
                }
            }
        }
        
        
        switch(numOfElements){
            case 0:
                System.err.println("No element has been found on boundary for InternalPoint with id: "+IPid+", coordinates of point x= "+this.theResultPoints.get(IPid).getCoordinates()[0]+", y= "+this.theResultPoints.get(IPid).getCoordinates()[1]+", tolerance being used: "+tol);
                val=0.;
                sxx=0.0; syy=0.0; sxy=0.0;
                break;
            case 1:
                // It is on element
                for(int i=0;i<onElemID.length;i++){
                    eyl=(this.theElements.get(onElemID[i]).getDispLocalonNode(2,step,wstate)[1]-
                            this.theElements.get(onElemID[i]).getDispLocalonNode(1,step,wstate)[1])/
                            ((ELine) this.theElements.get(onElemID[i])).getLength();
                    
                    for(int n=1;n<=theElements.get(onElemID[i]).getNumNodes();n++){
                        txl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[0];
                        tyl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[1];
                        sxx_l=txl;
                        sxy_l=tyl;
                        syy_l=1./(1.-v)*((E/(1.+v))*eyl+v*sxx_l);
                        if( ((StaticElasticity2DFS) this.theFundSol).getPlaneStress() )syy_l=(E*eyl+v*sxx_l);
                        c=((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[0];
                        s=-((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[1];
//                        System.out.println("element: "+onElemID[i]+" eyl: "+eyl+" tyl: "+tyl+" txl: "+txl+" sxx_l: "+sxx_l+" syy_l: "+syy_l+" sxy_l: "+sxy_l);
                        switch(wstress){
                            case 0:
                                val+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                break;
                            case 1:
                                val+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                break;
                            case 2: 
                                val+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                break;
                            case 6:
                                sxx+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                syy+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                sxy+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                        }
                    }
                }
                break;
            default:
                // It coinsides with node
                int NodeID = 0;
                for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
                    Node theNode = it.next();
                    dist=theNode.getDist(this.theResultPoints.get(IPid));
                    if(dist<=tol){
                        NodeID=theNode.getID();
                    }
                }
                
                for(int i=0;i<onElemID.length;i++){
                    eyl=(this.theElements.get(onElemID[i]).getDispLocalonNode(2,step,wstate)[1]-
                            this.theElements.get(onElemID[i]).getDispLocalonNode(1,step,wstate)[1])/
                            ((ELine) this.theElements.get(onElemID[i])).getLength();
                    
                    int n=theElements.get(onElemID[i]).getHierOfNode(NodeID);
                    txl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[0];
                    tyl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[1];
                    sxx_l=txl;
                    sxy_l=tyl;
                    syy_l=1./(1.-v)*((E/(1.+v))*eyl+v*sxx_l);
                    if( ((StaticElasticity2DFS) this.theFundSol).getPlaneStress() )syy_l=(E*eyl+v*sxx_l);
                    c=((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[0];
                    s=-((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[1];
                        
                    switch(wstress){
                        case 0:
                            val+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l);
                            break;
                        case 1:
                            val+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l);
                            break;
                        case 2: 
                            val+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l);
                            break;
                        case 6: 
                            sxx+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l);
                            sxy+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l);
                            sxy+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l);
                            break;
                    }
                    
                }
                val/=onElemID.length;
                sxx/=onElemID.length;
                syy/=onElemID.length;
                sxy/=onElemID.length;
                break;
        }
        if(wstress==6){
            double szz=((ElasticMat)this.getMaterial()).getLame_L()*(sxx+syy);
            if(((StaticElasticity2DFS)this.getFundamentalSolution()).getPlaneStress())szz=0.0;
            val=Math.sqrt(0.5*(
                (sxx-syy)*(sxx-syy)
                +(syy-szz)*(syy-szz)
                +(szz-sxx)*(szz-sxx)
                +6.0*sxy*sxy));
        }
        return val;
    }
    
    public double[][] getInternalPointStressOnBoundary(int IPid, int step, int wstate){
        // wstress  2D: 0 σxx, 1 σyy, 2 τxy
        //          3D: 0 σxx, 1 σyy, 2 σzz, 3 τxy, 4 τyz, 5 τxz //  not implemented yet.
        double[][] val = new double[3][3];
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                val[i][j]=0.0;
            }
        }
        int[] onElemID;
        double dist,cor;
        int numOfElements=0;
        Element elem;
        double eyl=0., tyl=0., txl=0.;
        double sxx_l=0.,syy_l=0.,sxy_l=0.;
        double s,c;
        double v=((ElasticMat) this.theMaterial).getPoissonRatio();
        //double G=((ElasticMat) this.theMaterial).getLame_M();
        double E=((ElasticMat) this.theMaterial).getElasticModulus();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            if(this.getFundamentalSolution().ndofs==2){
                dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
                cor=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
                if(cor>=-1. && cor<=1. && dist<=tol){
                    numOfElements+=1;
                }
            }else{
                if(elem.PointOnElement(theResultPoints.get(IPid))){
                    numOfElements+=1;
                }
            }
        }
        onElemID = new int[numOfElements];
        numOfElements=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            if(this.getFundamentalSolution().ndofs==2){
                dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
                cor=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
                if(cor>=-1. && cor<=1. && dist<=tol){
                    onElemID[numOfElements]=elem.getID();
                    numOfElements+=1;
                }
            }else{
                if(elem.PointOnElement(theResultPoints.get(IPid))){
                    numOfElements+=1;
                }
            }
        }
        
        
        switch(numOfElements){
            case 0:
                System.err.println("No element has been found on boundary for InternalPoint with id: "+IPid+", coordinates of point x= "+this.theResultPoints.get(IPid).getCoordinates()[0]+", y= "+this.theResultPoints.get(IPid).getCoordinates()[1]+", tolerance being used: "+tol);
                break;
            case 1:
                // It is on element
                for(int i=0;i<onElemID.length;i++){
                    eyl=(this.theElements.get(onElemID[i]).getDispLocalonNode(2,step,wstate)[1]-
                            this.theElements.get(onElemID[i]).getDispLocalonNode(1,step,wstate)[1])/
                            ((ELine) this.theElements.get(onElemID[i])).getLength();
                    
                    for(int n=1;n<=theElements.get(onElemID[i]).getNumNodes();n++){
                        txl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[0];
                        tyl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[1];
                        sxx_l=txl;
                        sxy_l=tyl;
                        syy_l=1./(1.-v)*((E/(1.+v))*eyl+v*sxx_l);
                        if( ((StaticElasticity2DFS) this.theFundSol).getPlaneStress() )syy_l=(E*eyl+v*sxx_l);
                        c=((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[0];
                        s=-((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[1];
//                        System.out.println("element: "+onElemID[i]+" eyl: "+eyl+" tyl: "+tyl+" txl: "+txl+" sxx_l: "+sxx_l+" syy_l: "+syy_l+" sxy_l: "+sxy_l);
                        val[0][0]+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                        val[1][1]+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                        val[0][1]+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                    }
                }
                break;
            default:
                // It coinsides with node
                int NodeID = 0;
                for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
                    Node theNode = it.next();
                    dist=theNode.getDist(this.theResultPoints.get(IPid));
                    if(dist<=tol){
                        NodeID=theNode.getID();
                    }
                }
                
                for(int i=0;i<onElemID.length;i++){
                    eyl=(this.theElements.get(onElemID[i]).getDispLocalonNode(2,step,wstate)[1]-
                            this.theElements.get(onElemID[i]).getDispLocalonNode(1,step,wstate)[1])/
                            ((ELine) this.theElements.get(onElemID[i])).getLength();
                    
                    int n=theElements.get(onElemID[i]).getHierOfNode(NodeID);
                    txl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[0];
                    tyl=theElements.get(onElemID[i]).getTractionLocalonNode(n,step,wstate)[1];
                    sxx_l=txl;
                    sxy_l=tyl;
                    syy_l=1./(1.-v)*((E/(1.+v))*eyl+v*sxx_l);
                    if( ((StaticElasticity2DFS) this.theFundSol).getPlaneStress() )syy_l=(E*eyl+v*sxx_l);
                    c=((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[0];
                    s=-((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[1];
                    val[0][0]+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l);
                    val[1][1]+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l);
                    val[0][1]+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l);
                }
                val[0][0]/=onElemID.length;
                val[1][1]/=onElemID.length;
                val[1][0]/=onElemID.length;
                val[0][1]=val[1][0]; 
                break;
        }
        return val;
    }
    
    public double getInternalPointStressOnBoundaryAux(int IPid, int wstress, int step, int wstate){
        // wstress: 0 σxx, 1 σyy, 2 τxy
        double val=0.;
        int[] onElemID;
        double dist,cor;
        int numOfElements=0;
        Element elem;
        double eyl=0., tyl=0., txl=0.;
        double sxx_l=0.,syy_l=0.,sxy_l=0.;
        double s,c;
        double v=((ElasticMat) this.theMaterial).getPoissonRatio();
        //double G=((ElasticMat) this.theMaterial).getLame_M();
        double E=((ElasticMat) this.theMaterial).getElasticModulus();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
            cor=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
            if(cor>=-1. && cor<=1. && dist<=tol){
                numOfElements+=1;
            }
        }
        onElemID = new int[numOfElements];
        numOfElements=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            elem = it.next();
            dist=elem.getMinDistOfNode(this.theResultPoints.get(IPid));
            cor=elem.getCoordMinDistOfNode(this.theResultPoints.get(IPid))[0];
            if(cor>=-1. && cor<=1. && dist<=tol){
                onElemID[numOfElements]=elem.getID();
                numOfElements+=1;
            }
        }
        
        
        switch(numOfElements){
            case 0:
                System.err.println("No element has been found on boundary for InternalPoint with id: "+IPid+", coordinates of point x= "+this.theResultPoints.get(IPid).getCoordinates()[0]+", y= "+this.theResultPoints.get(IPid).getCoordinates()[1]);
                val=0.;
                break;
            case 1:
                // It is on element
                for(int i=0;i<onElemID.length;i++){
                    eyl=(this.theElements.get(onElemID[i]).getDispLocalonNodeAux(2,step,wstate)[1]-
                            this.theElements.get(onElemID[i]).getDispLocalonNodeAux(1,step,wstate)[1])/
                            ((ELine) this.theElements.get(onElemID[i])).getLength();
                    
                    for(int n=1;n<=theElements.get(onElemID[i]).getNumNodes();n++){
                        txl=theElements.get(onElemID[i]).getTractionLocalonNodeAux(n,step,wstate)[0];
                        tyl=theElements.get(onElemID[i]).getTractionLocalonNodeAux(n,step,wstate)[1];
                        sxx_l=txl;
                        sxy_l=tyl;
                        syy_l=1./(1.-v)*((E/(1.+v))*eyl+v*sxx_l);
                        if( ((StaticElasticity2DFS) this.theFundSol).getPlaneStress() )syy_l=(E*eyl+v*sxx_l);
                        c=((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[0];
                        s=-((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[1];
//                        System.out.println("element: "+onElemID[i]+" eyl: "+eyl+" tyl: "+tyl+" txl: "+txl+" sxx_l: "+sxx_l+" syy_l: "+syy_l+" sxy_l: "+sxy_l);
                        switch(wstress){
                            case 0:
                                val+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                break;
                            case 1:
                                val+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                break;
                            case 2: 
                                val+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l)*((ELine) theElements.get(onElemID[i])).getShapeFunction(n, theElements.get(onElemID[i]).getLocalCoordinates(this.theResultPoints.get(IPid).getCoordinates())[0]);
                                break;
                        }
                    }
                }
                break;
            default:
                // It coinsides with node
                int NodeID = 0;
                for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
                    Node theNode = it.next();
                    dist=theNode.getDist(this.theResultPoints.get(IPid));
                    if(dist<=tol){
                        NodeID=theNode.getID();
                    }
                }
                
                for(int i=0;i<onElemID.length;i++){
                    eyl=(this.theElements.get(onElemID[i]).getDispLocalonNodeAux(2,step,wstate)[1]-
                            this.theElements.get(onElemID[i]).getDispLocalonNodeAux(1,step,wstate)[1])/
                            ((ELine) this.theElements.get(onElemID[i])).getLength();
                    
                    int n=theElements.get(onElemID[i]).getHierOfNode(NodeID);
                    txl=theElements.get(onElemID[i]).getTractionLocalonNodeAux(n,step,wstate)[0];
                    tyl=theElements.get(onElemID[i]).getTractionLocalonNodeAux(n,step,wstate)[1];
                    sxx_l=txl;
                    sxy_l=tyl;
                    syy_l=1./(1.-v)*((E/(1.+v))*eyl+v*sxx_l);
                    if( ((StaticElasticity2DFS) this.theFundSol).getPlaneStress() )syy_l=(E*eyl+v*sxx_l);
                    c=((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[0];
                    s=-((ELine) theElements.get(onElemID[i])).getNormal(theElements.get(onElemID[i]).getNodeHier(n).getID())[1];
                        
                    switch(wstress){
                        case 0:
                            val+=(c*c*sxx_l+s*s*syy_l+2.*c*s*sxy_l);
                            break;
                        case 1:
                            val+=(s*s*sxx_l+c*c*syy_l-2.*c*s*sxy_l);
                            break;
                        case 2: 
                            val+=(-c*s*sxx_l+c*s*syy_l+(c*c-s*s)*sxy_l);
                            break;
                    }
                    
                }
                val/=onElemID.length;
                break;
        }
        
        return val;
    }
    
    public double getInternalPointStrain(int IPid, int wstrain, int step, int wstate) {
        // wstrain: 0 εxx, 1 εyy, 2 γxy
        double val=0.;
        AbstractMatrix StressVec,CMatrix;
        switch(this.theFundSol.ndofs){
            case 2:
                if( ((StaticElasticity2DFS) theFundSol).getPlaneStress() ){
                    System.err.println("Strains for plane stress case do not implmented yet"); break;
                }
                StressVec = new AbstractMatrix(3,1); StressVec.init();
                StressVec.set(0, 0, getInternalPointStress(IPid,0,step,wstate));
                StressVec.set(1, 0, getInternalPointStress(IPid,1,step,wstate));
                StressVec.set(2, 0, getInternalPointStress(IPid,2,step,wstate));
                CMatrix = new AbstractMatrix(3,3); CMatrix.init();
                CMatrix.set(0, 0, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(0, 1, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 1, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 0, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(2, 2, 1./((ElasticMat) this.theMaterial).getLame_M() );
                val=(CMatrix.times(StressVec)).get(wstrain, 0);
                break;
            case 3: 
                break;
            default: System.err.println("Strains not implemented for "+this.theFundSol.ndofs+ "dimensions");
                    break;
        }
        return val;
    }
    
    public double getInternalPointStrainOnBoundary(int IPid, int wstrain, int step, int wstate) {
        // wstrain: 0 εxx, 1 εyy, 2 γxy
        double val=0.;
        AbstractMatrix StressVec,CMatrix;
        switch(this.theFundSol.ndofs){
            case 2:
                if( ((StaticElasticity2DFS) theFundSol).getPlaneStress() ){
                    System.err.println("Strains for plane stress case do not implmented yet"); break;
                }
                StressVec = new AbstractMatrix(3,1); StressVec.init();
                StressVec.set(0, 0, getInternalPointStressOnBoundary(IPid,0,step,wstate));
                StressVec.set(1, 0, getInternalPointStressOnBoundary(IPid,1,step,wstate));
                StressVec.set(2, 0, getInternalPointStressOnBoundary(IPid,2,step,wstate));
                CMatrix = new AbstractMatrix(3,3); CMatrix.init();
                CMatrix.set(0, 0, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(0, 1, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 1, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 0, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(2, 2, 1./((ElasticMat) this.theMaterial).getLame_M() );
                val=(CMatrix.times(StressVec)).get(wstrain, 0);
                break;
            case 3: 
                break;
            default: System.err.println("Strains not implemented for "+this.theFundSol.ndofs+ "dimensions");
                    break;
        }
        return val;
    }
    
    public double getInternalPointEnergyDensity(int IPid, int step, int wstate){
        double val=0.;
        AbstractMatrix StrainVec,StressVec,CMatrix;
        switch(this.theFundSol.ndofs){
            case 2:
                if( ((StaticElasticity2DFS) theFundSol).getPlaneStress() ){
                    System.err.println("EnergyDensity for plane stress case do not implmented yet"); break;
                }
                StressVec = new AbstractMatrix(3,1); StressVec.init();
                StressVec.set(0, 0, getInternalPointStress(IPid,0,step,wstate));
                StressVec.set(1, 0, getInternalPointStress(IPid,1,step,wstate));
                StressVec.set(2, 0, getInternalPointStress(IPid,2,step,wstate));
                CMatrix = new AbstractMatrix(3,3); CMatrix.init();
                CMatrix.set(0, 0, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(0, 1, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 1, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 0, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(2, 2, 1./((ElasticMat) this.theMaterial).getLame_M() );
                StrainVec = new AbstractMatrix(3,1); StrainVec.init();
                StrainVec.set(0, 0, (CMatrix.times(StressVec)).get(0, 0));
                StrainVec.set(1, 0, (CMatrix.times(StressVec)).get(1, 0));
                StrainVec.set(2, 0, (CMatrix.times(StressVec)).get(2, 0));
                val= ((StressVec.transpose()).times(StrainVec)).get(0, 0);
                val/=2.;
                break;
            case 3: 
                break;
            default: System.err.println("EnergyDensity not implemented for "+this.theFundSol.ndofs+ "dimensions");
                    break;
        }
        return val;
    }
    
    public double getInternalPointEnergyDensityOnBoundary(int IPid, int step, int wstate){
        double val=0.;
        AbstractMatrix StrainVec,StressVec,CMatrix;
        switch(this.theFundSol.ndofs){
            case 2:
                if( ((StaticElasticity2DFS) theFundSol).getPlaneStress() ){
                    System.err.println("EnergyDensity for plane stress case do not implmented yet"); break;
                }
                StressVec = new AbstractMatrix(3,1); StressVec.init();
                StressVec.set(0, 0, getInternalPointStressOnBoundary(IPid,0,step,wstate));
                StressVec.set(1, 0, getInternalPointStressOnBoundary(IPid,1,step,wstate));
                StressVec.set(2, 0, getInternalPointStressOnBoundary(IPid,2,step,wstate));
                CMatrix = new AbstractMatrix(3,3); CMatrix.init();
                CMatrix.set(0, 0, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(0, 1, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 1, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 0, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(2, 2, 1./((ElasticMat) this.theMaterial).getLame_M() );
                StrainVec = new AbstractMatrix(3,1); StrainVec.init();
                StrainVec.set(0, 0, (CMatrix.times(StressVec)).get(0, 0));
                StrainVec.set(1, 0, (CMatrix.times(StressVec)).get(1, 0));
                StrainVec.set(2, 0, (CMatrix.times(StressVec)).get(2, 0));
                val= ((StressVec.transpose()).times(StrainVec)).get(0, 0);
                val/=2.;
                break;
            case 3: 
                break;
            default: System.err.println("EnergyDensity not implemented for "+this.theFundSol.ndofs+ "dimensions");
                    break;
        }
        return val;
    }
    
    private double getInternalPointEnergyDensity(int IPid, int stress_step, int strain_step, int wstate){
        double val=0.;
        AbstractMatrix StrainVec,StressVec,CMatrix,StressVec_onStrainStep;
        switch(this.theFundSol.ndofs){
            case 2:
                if( ((StaticElasticity2DFS) theFundSol).getPlaneStress() ){
                    System.err.println("EnergyDensity for plane stress case do not implmented yet"); break;
                }
                StressVec = new AbstractMatrix(3,1); StressVec.init();
                StressVec.set(0, 0, getInternalPointStress(IPid,0,stress_step,wstate));
                StressVec.set(1, 0, getInternalPointStress(IPid,1,stress_step,wstate));
                StressVec.set(2, 0, getInternalPointStress(IPid,2,stress_step,wstate));
                
                StressVec_onStrainStep = new AbstractMatrix(3,1); StressVec_onStrainStep.init();
                StressVec_onStrainStep.set(0, 0, getInternalPointStress(IPid,0,strain_step,wstate));
                StressVec_onStrainStep.set(1, 0, getInternalPointStress(IPid,1,strain_step,wstate));
                StressVec_onStrainStep.set(2, 0, getInternalPointStress(IPid,2,strain_step,wstate));
                
                CMatrix = new AbstractMatrix(3,3); CMatrix.init();
                CMatrix.set(0, 0, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(0, 1, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 1, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 0, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(2, 2, 1./((ElasticMat) this.theMaterial).getLame_M() );
                StrainVec = new AbstractMatrix(3,1); StrainVec.init();
                StrainVec.set(0, 0, (CMatrix.times(StressVec_onStrainStep)).get(0, 0));
                StrainVec.set(1, 0, (CMatrix.times(StressVec_onStrainStep)).get(1, 0));
                StrainVec.set(2, 0, (CMatrix.times(StressVec_onStrainStep)).get(2, 0));
                val= ((StressVec.transpose()).times(StrainVec)).get(0, 0);
                //val/=2.;
                break;
            case 3: 
                break;
            default: System.err.println("EnergyDensity not implemented for "+this.theFundSol.ndofs+ "dimensions");
                    break;
        }
        return val;
    }
    
    private double getInternalPointEnergyDensityOnBoundary(int IPid, int stress_step, int strain_step, int wstate){
        double val=0.;
        AbstractMatrix StrainVec,StressVec,CMatrix,StressVec_onStrainStep;
        switch(this.theFundSol.ndofs){
            case 2:
                if( ((StaticElasticity2DFS) theFundSol).getPlaneStress() ){
                    System.err.println("EnergyDensity for plane stress case do not implmented yet"); break;
                }
                StressVec = new AbstractMatrix(3,1); StressVec.init();
                StressVec.set(0, 0, getInternalPointStressOnBoundary(IPid,0,stress_step,wstate));
                StressVec.set(1, 0, getInternalPointStressOnBoundary(IPid,1,stress_step,wstate));
                StressVec.set(2, 0, getInternalPointStressOnBoundary(IPid,2,stress_step,wstate));
                
                StressVec_onStrainStep = new AbstractMatrix(3,1); StressVec_onStrainStep.init();
                StressVec_onStrainStep.set(0, 0, getInternalPointStressOnBoundary(IPid,0,strain_step,wstate));
                StressVec_onStrainStep.set(1, 0, getInternalPointStressOnBoundary(IPid,1,strain_step,wstate));
                StressVec_onStrainStep.set(2, 0, getInternalPointStressOnBoundary(IPid,2,strain_step,wstate));
                
                CMatrix = new AbstractMatrix(3,3); CMatrix.init();
                CMatrix.set(0, 0, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(0, 1, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 1, 1./((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(1, 0, -((ElasticMat) this.theMaterial).getPoissonRatio()/((ElasticMat) this.theMaterial).getElasticModulus() );
                CMatrix.set(2, 2, 1./((ElasticMat) this.theMaterial).getLame_M() );
                StrainVec = new AbstractMatrix(3,1); StrainVec.init();
                StrainVec.set(0, 0, (CMatrix.times(StressVec_onStrainStep)).get(0, 0));
                StrainVec.set(1, 0, (CMatrix.times(StressVec_onStrainStep)).get(1, 0));
                StrainVec.set(2, 0, (CMatrix.times(StressVec_onStrainStep)).get(2, 0));
                val= ((StressVec.transpose()).times(StrainVec)).get(0, 0);
                //val/=2.;
                break;
            case 3: 
                break;
            default: System.err.println("EnergyDensity not implemented for "+this.theFundSol.ndofs+ "dimensions");
                    break;
        }
        return val;
    }
    
    public double getInternalPointTIEnergyDensityRate(int IPid, double tau, int step, int wstate){
        double val=0.;
        if(step>0)
            val =   getInternalPointEnergyDensity(IPid, step, step, wstate)
                    -getInternalPointEnergyDensity(IPid, (step-1), step, wstate)
                    -getInternalPointEnergyDensity(IPid, step, (step-1), wstate)
                    +getInternalPointEnergyDensity(IPid, (step-1), (step-1), wstate)
                    ;
        val/=tau;
        return val;
    }
    
    public double getInternalPointTIEnergyDensityRateOnBoundary(int IPid, double tau, int step, int wstate){
        double val=0.;
        if(step>0)
            val =   getInternalPointEnergyDensityOnBoundary(IPid, step, step, wstate)
                    -getInternalPointEnergyDensityOnBoundary(IPid, (step-1), step, wstate)
                    -getInternalPointEnergyDensityOnBoundary(IPid, step, (step-1), wstate)
                    +getInternalPointEnergyDensityOnBoundary(IPid, (step-1), (step-1), wstate)
                    ;
        val/=tau;
        return val;
    }
    
    public int[][] FindCommonNodes(Domain anotherDomain, double tolerance){
        // returns a "Matrix" (actually int[][]) where in each serie has in the first column the node of thisDomain
        // and in the second column the respective (of approximataly the same coordinates) node of theAnotherDomain
        int[][] common = null;
        //System.out.println("Common Nodes of Domains:"+this.getID()+" - "+anotherDomain.getID());
        int numOfCommonNodes=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node NodeOfThisDomain = it.next();
            for(Iterator<Node> it2=anotherDomain.theNodes.values().iterator(); it2.hasNext();){
                Node NodeOfAnotherDomain = it2.next();
                if(NodeOfThisDomain.getDist(NodeOfAnotherDomain)<=tolerance)numOfCommonNodes+=1;
            }
        }
        //System.out.println("Number of common nodes: "+numOfCommonNodes);
        common = new int[numOfCommonNodes][2];
        numOfCommonNodes=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node NodeOfThisDomain = it.next();
            for(Iterator<Node> it2=anotherDomain.theNodes.values().iterator(); it2.hasNext();){
                Node NodeOfAnotherDomain = it2.next();
                if(NodeOfThisDomain.getDist(NodeOfAnotherDomain)<=tolerance){
                    common[numOfCommonNodes][0]=NodeOfThisDomain.getID();
                    common[numOfCommonNodes][1]=NodeOfAnotherDomain.getID();
                    numOfCommonNodes+=1;
                }
            }
        }
        return common;
    }
    
    public int[][] FindCommonElements(Domain anotherDomain, double tolerance){
        // returns a "Matrix" (actually int[][]) where in each serie has in the first column the element of thisDomain
        // and in the second column the respective  element of theAnotherDomain
        // in the third the node of elemen of thisDomain
        // in the fourth the respective node of elemen of theAnotherDomain
        int[][] common = null;
        int numOfCommonNodes=0;
        Node Node_1_elem_1;
        Node Node_2_elem_1;
        Node Node_1_elem_2;
        Node Node_2_elem_2;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            Node_1_elem_1 = elem.getNodeHier(1);
            Node_2_elem_1 = elem.getNodeHier(2);
            for(Iterator<Element> it2=anotherDomain.getElements().values().iterator(); it2.hasNext();){
                Element elem2 = it2.next();
                Node_1_elem_2 = elem2.getNodeHier(1);
                Node_2_elem_2 = elem2.getNodeHier(2);
                if( ((Node_1_elem_1.getDist(Node_1_elem_2)<=tolerance)&&(Node_2_elem_1.getDist(Node_2_elem_2)<=tolerance))||
                    ((Node_2_elem_1.getDist(Node_1_elem_2)<=tolerance)&&(Node_1_elem_1.getDist(Node_2_elem_2)<=tolerance)))numOfCommonNodes+=2;
            }
        }
        common = new int[numOfCommonNodes][4];
        
        
        numOfCommonNodes=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            Node_1_elem_1 = elem.getNodeHier(1);
            Node_2_elem_1 = elem.getNodeHier(2);
            for(Iterator<Element> it2=anotherDomain.getElements().values().iterator(); it2.hasNext();){
                Element elem2 = it2.next();
                Node_1_elem_2 = elem2.getNodeHier(1);
                Node_2_elem_2 = elem2.getNodeHier(2);
                if( ((Node_1_elem_1.getDist(Node_1_elem_2)<=tolerance)&&(Node_2_elem_1.getDist(Node_2_elem_2)<=tolerance))||
                    ((Node_2_elem_1.getDist(Node_1_elem_2)<=tolerance)&&(Node_1_elem_1.getDist(Node_2_elem_2)<=tolerance))){
                    common[numOfCommonNodes][0]=elem.getID();
                    common[numOfCommonNodes][1]=elem2.getID();
                    common[numOfCommonNodes+1][0]=elem.getID();
                    common[numOfCommonNodes+1][1]=elem2.getID();
                    if((Node_1_elem_1.getDist(Node_1_elem_2)<=tolerance)&&(Node_2_elem_1.getDist(Node_2_elem_2)<=tolerance)){
                        common[numOfCommonNodes][2]=Node_1_elem_1.getID();
                        common[numOfCommonNodes][3]=Node_1_elem_2.getID();
                        
                        common[numOfCommonNodes+1][2]=Node_2_elem_1.getID();
                        common[numOfCommonNodes+1][3]=Node_2_elem_2.getID();
                    }
                    else{
                        common[numOfCommonNodes][2]=Node_2_elem_1.getID();
                        common[numOfCommonNodes][3]=Node_1_elem_2.getID();
                        
                        common[numOfCommonNodes+1][2]=Node_1_elem_1.getID();
                        common[numOfCommonNodes+1][3]=Node_2_elem_2.getID();
                    }
                    numOfCommonNodes+=2;
                }
            }
        }
        return common;
    }
    
    public boolean isMain(){
        return this.getElements().values().iterator().next().isMain();
    }
    
    public void setUniformTempChange(double val){this.uniformTempChange=val;}
    
    public double getUniformTempChange(){return this.uniformTempChange;}
    
    public double maximum_X_coordinate(){return this.maximum_X_coordinate;}
    
    public double maximum_Y_coordinate(){return this.maximum_Y_coordinate;}
    
    public double maximum_Z_coordinate(){return this.maximum_Z_coordinate;}
    
    public double minimum_X_coordinate(){return this.minimum_X_coordinate;}
    
    public double minimum_Y_coordinate(){return this.minimum_Y_coordinate;}
    
    public double minimum_Z_coordinate(){return this.minimum_Z_coordinate;}
    
    public void SetDirichletBCs(Line2D theLine, DoubleFunction df){
        SetDirichletBCs(theLine,df,1,true,null, 1, 1.0);
    }
    
    public void SetDirichletBCs(Line2D theLine, DoubleFunction df, boolean setequal){
        SetDirichletBCs(theLine,df,1,setequal, null, 1, 1.0);
    }
    
    public void SetDirichletBCs(Line2D theLine, DoubleFunction df, int wdof){
        SetDirichletBCs(theLine, df, wdof, true, null, 1, 1.0);
    }
    
    public void SetDirichletBCs(Line2D theLine, DoubleFunction df, int wdof, TimeFunction tf, int steps, double tau){
        SetDirichletBCs(theLine, df, wdof, true, tf, steps, tau);
    }
    
    public void SetDirichletBCs(Line2D theLine, DoubleFunction df, int wdof, boolean setequal, TimeFunction tf, int steps,double tau){
        ConstraintEquation aConstraintEquation=new ConstraintEquation();
        ConstraintTerm aConstraintTerm;
        double[] cvals= new double[steps+1];
        for(int i=0;i<steps+1;i++){
            if(tf!=null){
                cvals[i]=tf.call(i*tau);
            }else{
                cvals[i]=1.0;
            }
        }
        if(this.theMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : this.theMaterial.getClass().toString().equals(ViscousMaterial.class.toString())){
            cvals=((ViscousMaterial)theMaterial).TransformBoundaryDisplacement(cvals, tau);
        }
        int cid=aConstraintEquation.getNumOfConstraintEquations();
        // find also nearest Node to first Point of Line2D
        int strID = 0,endID = 0;
        double minstr=Double.POSITIVE_INFINITY;
        double minend=Double.POSITIVE_INFINITY;
        geom.Point Pstart=theLine.getPstart();
        geom.Point Pend=theLine.getPend();
        //////////////////////////////////////////////////
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theLine.contains(theNode)){
                double[] tvals = new double[cvals.length];
                for(int i=0;i<steps+1;i++){
                    tvals[i]=cvals[i]*(df.run(theNode.X(), theNode.Y()));
                }
                aConstraintEquation = new ConstraintEquation(++cid,tvals);
                aConstraintTerm = new ConstraintTerm(theNode, 1, wdof, 1.0);
                aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
                this.putConstraintEquation(aConstraintEquation);
                if(theNode.getDist(Pstart)<minstr){strID=theNode.getID();minstr=theNode.getDist(Pstart);}
                if(theNode.getDist(Pend)<minend){endID=theNode.getID();minend=theNode.getDist(Pend);}
            }
        }
        //set equal tractions of adjacent Elements joined on Node
        if(setequal)for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theLine.contains(theNode)){
                if(theNode.getID()!=strID && theNode.getID()!=endID){
                    int[] elems=theNode.getConnectedElementsIds();
                    aConstraintEquation = new ConstraintEquation(++cid);
                    aConstraintEquation.equalp(this.getElement(elems[0]), this.getElement(elems[1]), theNode.getID(), wdof);
                    this.putConstraintEquation(aConstraintEquation);
                }
            }
        }
    }
    ///////// same as above but with const val
    public void SetDirichletBCs(Line2D theLine, double val){
        SetDirichletBCs(theLine,val,1);
    }
    
    public void SetDirichletBCs(Line2D theLine, double val, boolean setequal){
        SetDirichletBCs(theLine,val,1,setequal,null,1,1.0);
    }
    
    public void SetDirichletBCs(Line2D theLine, double val, int wdof){
        SetDirichletBCs(theLine, val, wdof, true,null,1,1.0);
    }
    
    public void SetDirichletBCs(Line2D theLine, double val, int wdof, TimeFunction tf, int steps, double tau){
        SetDirichletBCs(theLine, val, wdof, true, tf, steps, tau);
    }
    
    public void SetDirichletBCs(Line2D theLine, double val, int wdof, boolean setequal, TimeFunction tf, int steps,double tau){        
        ConstraintEquation aConstraintEquation=new ConstraintEquation();
        ConstraintTerm aConstraintTerm;
        double[] cvals= new double[steps+1];
        for(int i=0;i<steps+1;i++){
            if(tf!=null){
                cvals[i]=tf.call(i*tau)*val;
            }else{
                cvals[i]=val;
            }
        }
        if(this.theMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : this.theMaterial.getClass().toString().equals(ViscousMaterial.class.toString())){
            cvals=((ViscousMaterial)theMaterial).TransformBoundaryDisplacement(cvals, tau);
        }
        int cid=aConstraintEquation.getNumOfConstraintEquations();
        // find also nearest Node to first Point of Line2D
        int strID = 0,endID = 0;
        double minstr=Double.POSITIVE_INFINITY;
        double minend=Double.POSITIVE_INFINITY;
        Point Pstart=theLine.getPstart();
        Point Pend=theLine.getPend();
        //////////////////////////////////////////////////
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theLine.contains(theNode)){
                aConstraintEquation = new ConstraintEquation(++cid,cvals);
                aConstraintTerm = new ConstraintTerm(theNode, 1, wdof, 1.0);
                aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
                this.putConstraintEquation(aConstraintEquation);
                if(theNode.getDist(Pstart)<minstr){strID=theNode.getID();minstr=theNode.getDist(Pstart);}
                if(theNode.getDist(Pend)<minend){endID=theNode.getID();minend=theNode.getDist(Pend);}
            }
        }
        //set equal tractions of adjacent Elements joined on Node
        if(setequal)for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theLine.contains(theNode)){
                if(theNode.getID()!=strID && theNode.getID()!=endID){
                    int[] elems=theNode.getConnectedElementsIds();
                    aConstraintEquation = new ConstraintEquation(++cid);
                    aConstraintEquation.equalp(this.getElement(elems[0]), this.getElement(elems[1]), theNode.getID(), wdof);
                    this.putConstraintEquation(aConstraintEquation);
                }
            }
        }
    }
    ///////////////////////////
    // on the whole group of boundary nodes
    ////
    public void SetDirichletBCs(DoubleFunction df){
        SetDirichletBCs(df,1,true,null, 1, 1.0);
    }
    
    public void SetDirichletBCs(DoubleFunction df, int wdof){
        SetDirichletBCs(df, wdof, true, null, 1, 1.0);
    }
    
    public void SetDirichletBCs(DoubleFunction df, int wdof, TimeFunction tf, int steps,double tau){
        SetDirichletBCs(df, wdof, true, tf, steps,tau);
    }
    
    public void SetDirichletBCs(DoubleFunction df, int wdof, boolean setequal, TimeFunction tf, int steps,double tau){
        ConstraintEquation aConstraintEquation=new ConstraintEquation();
        ConstraintTerm aConstraintTerm;
        double[] cvals= new double[steps+1];
        for(int i=0;i<steps+1;i++){
            if(tf!=null){
                cvals[i]=tf.call(i*tau);
            }else{
                cvals[i]=1.0;
            }
        }
        if(this.theMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : this.theMaterial.getClass().toString().equals(ViscousMaterial.class.toString())){
            cvals=((ViscousMaterial)theMaterial).TransformBoundaryDisplacement(cvals, tau);
        }
        int cid=aConstraintEquation.getNumOfConstraintEquations();
        // find also nearest Node to first Point of Line2D
        int strID = 0,endID = 0;
        //////////////////////////////////////////////////
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            double[] tvals = new double[cvals.length];
            for(int i=0;i<steps+1;i++){
                tvals[i]=cvals[i]*(df.run(theNode.X(), theNode.Y()));
            }
            aConstraintEquation = new ConstraintEquation(++cid,tvals);
            aConstraintTerm = new ConstraintTerm(theNode, 1, wdof, 1.0);
            aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);
            this.putConstraintEquation(aConstraintEquation);
        }
        //set equal tractions of adjacent Elements joined on Node
        if(setequal)for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theNode.getID()!=strID && theNode.getID()!=endID){
                int[] elems=theNode.getConnectedElementsIds();
                aConstraintEquation = new ConstraintEquation(++cid);
                aConstraintEquation.equalp(this.getElement(elems[0]), this.getElement(elems[1]), theNode.getID(), wdof);
                this.putConstraintEquation(aConstraintEquation);
            }
        }
    }
    
    public void SetNeumannBCs(DoubleFunction df_x, DoubleFunction df_y){
        SetNeumannBCs(df_x, df_y, 1,null,1,1.0);
    }
    
    public void SetNeumannBCs(DoubleFunction df_x, DoubleFunction df_y, int wdof, TimeFunction tf, int steps,double tau){
        ConstraintEquation aConstraintEquation=new ConstraintEquation();
        ConstraintTermElement aConstraintTermElement;
        double[] cvals= new double[steps+1];
        for(int i=0;i<steps+1;i++){
            if(tf!=null){
                cvals[i]=tf.call(i*tau);
            }else{
                cvals[i]=1.0;
            }
        }
        if(this.theMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : this.theMaterial.getClass().toString().equals(ViscousMaterial.class.toString())){
            cvals=((ViscousMaterial)theMaterial).TransformBoundaryTraction(cvals, tau);
        }
        int cid=aConstraintEquation.getNumOfConstraintEquations();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            boolean lieson=true;
            for(Iterator<Node> nit=elem.getNodes().values().iterator();nit.hasNext();){
                Node elNode=nit.next();
                double x=elNode.X(); double y=elNode.Y();
                double[] normal = elem.getNormal(elNode.getID());
//                    for(int i=0;i<steps+1;i++){
//                        cvals[i]*=normal[0]*df_x.run(x, y)+normal[1]*df_y.run(x,y);
//                    }
                double[] tvals = new double[cvals.length];
                for(int i=0;i<steps+1;i++){
                    tvals[i]=cvals[i]*(normal[0]*df_x.run(x, y)+normal[1]*df_y.run(x,y));
                }
                aConstraintEquation = new ConstraintEquation(++cid,tvals);
                aConstraintTermElement = new 
                    ConstraintTermElement(elem,
                    elNode.getID(), 
                    wdof, 1.0);
                aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                this.putConstraintEquation(aConstraintEquation);
            }
        }
    }
    ///////////////////////////
    
    public void SetNeumannBCs(Line2D theLine, DoubleFunction df_x, DoubleFunction df_y){
        SetNeumannBCs(theLine, df_x, df_y, 1,null,1,1.0);
    }
    
    public void SetNeumannBCs(Line2D theLine, DoubleFunction df_x, DoubleFunction df_y, int wdof, TimeFunction tf, int steps,double tau){
        ConstraintEquation aConstraintEquation=new ConstraintEquation();
        ConstraintTermElement aConstraintTermElement;
        double[] cvals= new double[steps+1];
        for(int i=0;i<steps+1;i++){
            if(tf!=null){
                cvals[i]=tf.call(i*tau);
            }else{
                cvals[i]=1.0;
            }
        }
        if(this.theMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : this.theMaterial.getClass().toString().equals(ViscousMaterial.class.toString())){
            cvals=((ViscousMaterial)theMaterial).TransformBoundaryTraction(cvals, tau);
        }
        int cid=aConstraintEquation.getNumOfConstraintEquations();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            boolean lieson=true;
            for(Iterator<Node> nit=elem.getNodes().values().iterator();nit.hasNext();){
                Node elNode=nit.next();
                if(!theLine.contains(elNode))lieson=false;
            }
            if(lieson){
                for(Iterator<Node> nit=elem.getNodes().values().iterator();nit.hasNext();){
                    Node elNode=nit.next();
                    double x=elNode.X(); double y=elNode.Y();
                    double[] normal = elem.getNormal(elNode.getID());
//                    for(int i=0;i<steps+1;i++){
//                        cvals[i]*=normal[0]*df_x.run(x, y)+normal[1]*df_y.run(x,y);
//                    }
                    double[] tvals = new double[cvals.length];
                    for(int i=0;i<steps+1;i++){
                        tvals[i]=cvals[i]*(normal[0]*df_x.run(x, y)+normal[1]*df_y.run(x,y));
                    }
                    aConstraintEquation = new ConstraintEquation(++cid,tvals);
                    aConstraintTermElement = new 
                        ConstraintTermElement(elem,
                        elNode.getID(), 
                        wdof, 1.0);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                    this.putConstraintEquation(aConstraintEquation);
                }
            }
        }
    }
    
    public void SetNeumannBCs(Line2D theLine, double val, int wdof){
        SetNeumannBCs(theLine, val, wdof, null, 1, 1.0);
    }
    
    public void SetNeumannBCs(Line2D theLine, double val, int wdof, TimeFunction tf, int steps,double tau){
        ConstraintEquation aConstraintEquation=new ConstraintEquation();
        ConstraintTermElement aConstraintTermElement;
        double[] cvals= new double[steps+1];
        for(int i=0;i<steps+1;i++){
            if(tf!=null){
                cvals[i]=tf.call(i*tau)*val;
            }else{
                cvals[i]=val;
            }
        }
        if(this.theMaterial.getClass().toString() == null ? ViscousMaterial.class.toString() == null : this.theMaterial.getClass().toString().equals(ViscousMaterial.class.toString())){
            cvals=((ViscousMaterial)theMaterial).TransformBoundaryTraction(cvals, tau);
        }
        int cid=aConstraintEquation.getNumOfConstraintEquations();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            boolean lieson=true;
            for(Iterator<Node> nit=elem.getNodes().values().iterator();nit.hasNext();){
                Node elNode=nit.next();
                if(!theLine.contains(elNode))lieson=false;
            }
            if(lieson){
                for(Iterator<Node> nit=elem.getNodes().values().iterator();nit.hasNext();){
                    Node elNode=nit.next();
                    double x=elNode.X(); double y=elNode.Y();
                    aConstraintEquation = new ConstraintEquation(++cid,cvals);
                    aConstraintTermElement = new 
                        ConstraintTermElement(elem,
                        elNode.getID(), 
                        wdof, 1.0);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);
                    this.putConstraintEquation(aConstraintEquation);
                }
            }
        }
    }
    
    public void SetNeumannBCs(Line2D theLine, double val){
        SetNeumannBCs(theLine, val, 1);
    }
    
    
    // some printing methods concerning response
    public void printResponse(Line2D theLine, int wdof, int step, int state){
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theLine.contains(theNode)){
                int nid=theNode.getID(); double x=theNode.X(); double y=theNode.Y();
                System.out.print(nid+" "+x+" "+y+" "+theNode.getu()[wdof-1][step][state]);
                int[] elems=theNode.getConnectedElementsIds();
                for(int i=0;i<elems.length;i++){
                    if(theLine.contains(this.getElement(elems[i]))){
                        double val=theNode.getp(getElement(elems[i]), this.theFundSol.get_p_DOFs(), state)[wdof-1][step];
                        System.out.print(" "+val);
                    }else{
                        System.out.print(" "+"-");
                    }
                }
                System.out.println();
            }
        }
    }
    
    public void printResponse(Line2D theLine){
        printResponse(theLine, 1, 0, 0);
    }
    
    public void printResponse(Line2D theLine, int wdof){
        printResponse(theLine, wdof, 0, 0);
    }
    
    public void printResponse(Line2D theLine, int wdof, int step){
        printResponse(theLine, wdof, step, 0);
    }
    
    public void testF(DoubleFunction df){
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            System.out.println(theNode.getID()+" "+df.run(theNode.X(), theNode.Y()));
        }
    }
    
    public void Triangulate(){
        //delaunay.Triangle.moreInfo = true;
        double mx,my,m;
        mx=Math.max(Math.abs(minimum_X_coordinate()), maximum_X_coordinate());
        my=Math.max(Math.abs(minimum_X_coordinate()), maximum_X_coordinate());
        m=Math.max(mx,my);
        delaunay.Triangle tri =
            new delaunay.Triangle(   new Pnt(-m*10000.0,m*10000.0), 
                            new Pnt(m*10000.0,m*10000.0), 
                            new Pnt(0,-m*10000.0));
        //System.out.println(tri);
        Triangulation dt = new Triangulation(tri);
        boolean repeat=false;
//        System.out.println("DelaunayTriangulation created: " + dt);
//        for(Iterator<Node>it=getNodes().values().iterator();it.hasNext();){
//            Node aNode=it.next();
//            Pnt aPnt = new Pnt(aNode.getCoordinates()[0],aNode.getCoordinates()[1]);
//            aPnt.setIndex(aNode.getID());
//            try{
//                dt.delaunayPlace(aPnt);
//            }catch(java.util.NoSuchElementException e){
//                repeat=true;
//            }
//        }
            
        for(Iterator<ResultPoint>it=this.getResultPoints().values().iterator();it.hasNext();){
            ResultPoint aNode=it.next();
            Pnt aPnt = new Pnt(aNode.getCoordinates()[0],aNode.getCoordinates()[1]);
            aPnt.setIndex(aNode.getID());
//            System.out.println("Point id = "+aNode.getID());
            try{
                dt.delaunayPlace(aPnt);
            }catch(java.util.NoSuchElementException e){
                repeat=true;
            }
        }
            
        int indtr=0;
        for (delaunay.Triangle triangle : dt) {
            if (Collections.disjoint(tri, triangle)) {
                Point[] nds = new Point[3];
                int count=0;
//                System.out.println("triangle: "+triangle.getID());
                for(Iterator<Pnt>it=triangle.iterator();it.hasNext();){
                    Pnt aPnt=it.next();
//                    System.out.print(aPnt.getIndex()+" ");
//                    if(theNodes.containsKey(aPnt.getIndex())){
//                        nds[count++] = this.theNodes.get(aPnt.getIndex());
//                    }else{
//                        nds[count++] = this.theResultPoints.get(aPnt.getIndex());
//                    }
                    nds[count++] = this.theResultPoints.get(aPnt.getIndex());
                }
//                System.out.println();
                Triangle aTrg = new Triangle(++indtr,nds);
                theShapes.put(aTrg.getID(), aTrg);
            }
        }
        isTriangulized=true;
    }
    
    public boolean isTriangulized(){return isTriangulized;}
    
    public Map<Integer,Shape> getShapes(){
        return this.theShapes;
    }
    
    public ArrayList<Node> getNodes(Line2D theLine){
        ArrayList<Node> NodesOnLine = new ArrayList<Node>();
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            if(theLine.contains(theNode)){
                NodesOnLine.add(theNode);
            }
        }
        return NodesOnLine;
    }
    
    public ArrayList<Element> getElements(Line2D theLine){
        ArrayList<Element> ElementsOnLine = new ArrayList<Element>();
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            if(theLine.contains(theElement)){
                ElementsOnLine.add(theElement);
            }
        }
        return ElementsOnLine;
    }
    
    public Map getInternalPointDisps(){
        return getInternalPointDisps(0, 0, 0);
    }
    
    public Map getInternalPointDisps(int step, int wstate){
        return getInternalPointDisps(0, step, wstate);
    }
    
    public Map getInternalPointDisps(int wdisp, int step, int wstate){
        TreeMap<Integer, Double> myMap = new TreeMap<Integer, Double>();
        ResultPoint[] InteriorPoint;
        ResultPoint[] BoundaryPoint;
        int Nbound=0;
        int Ninter=0;
        for(Iterator<ResultPoint> it=this.theResultPoints.values().iterator(); it.hasNext();){
            if(it.next().isOnBoundary())Nbound++;
        }
        BoundaryPoint = new ResultPoint[Nbound];
        InteriorPoint = new ResultPoint[theResultPoints.size()-Nbound];
        Nbound=0;
        for(Iterator<ResultPoint> it=this.theResultPoints.values().iterator(); it.hasNext();){
            ResultPoint theResultPoint = it.next();
            if(theResultPoint.isOnBoundary()){
                BoundaryPoint[Nbound++]=theResultPoint;
            }else{
                InteriorPoint[Ninter++]=theResultPoint;
            }
        }
//        if(this.theResultPoints.get(IPid).isOnBoundary()){
//            val=getInternalPointDispOnBoundary(IPid, wdisp, step, wstate);
//        }else{
//            val=getInternalPointDispInterior(IPid, wdisp, step, wstate);
//        }
        return myMap;
    }
}

