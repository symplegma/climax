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
//import java.util.Map;
//import java.util.HashMap;
import jmat.AbstractMatrix;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 *
 * @author pchr
 */
public class ConstraintElement extends Element{
    private static int numberOfConstraintElements = 0;
    private static AbstractMatrix MyMatrix;
    private double[] constraintValues;
    private double[] FirstDerconstraintValues;
    private double[] SecDerconstraintValues;
    private double ak=10.e12;
    private double am=10.e12;
    private double ac=0.0;
    protected List<Integer> NodeDOFS = new Vector() ;
    protected List<Double> coefs = new Vector() ;
    
    //constructors
    public ConstraintElement(){
    }
    
    public ConstraintElement(int id, double coef, Node aNode, int wdof){
        this.id=id;
        ++numberOfConstraintElements;
        ++this.ndofs;
        this.putNode(aNode);
        int[] dofs= new int[1];
        dofs[0]=wdof;
        aNode.setNdofs_ofNode(dofs);
        //this.NodeDOFS.put(this.numNodes, wdof);
        this.NodeDOFS.add(ndofs-1,wdof);
        this.coefs.add(ndofs-1,coef);
    }
    
    public ConstraintElement(int id, double coef, Node aNode, int wdof, double[] constraintValues){
        this.id=id;
        ++numberOfConstraintElements;
        ++this.ndofs;
        this.putNode(aNode);
        int[] dofs= new int[1];
        dofs[0]=wdof;
        aNode.setNdofs_ofNode(dofs);
        //this.NodeDOFS.put(this.numNodes, wdof);
        this.NodeDOFS.add(ndofs-1,wdof);
        this.coefs.add(ndofs-1,coef);
        this.constraintValues = new double[constraintValues.length];
        for(int i=0;i<constraintValues.length;i++){
            this.constraintValues[i]=constraintValues[i];
        }
        
    }
    
    public ConstraintElement(int id, double coef, Node aNode, int wdof,
                double[] constraintValues, double[] SecDerconstraintValues){
        this.id=id;
        ++numberOfConstraintElements;
        ++this.ndofs;
        this.putNode(aNode);
        int[] dofs= new int[1];
        dofs[0]=wdof;
        aNode.setNdofs_ofNode(dofs);
        //this.NodeDOFS.put(this.numNodes, wdof);
        this.NodeDOFS.add(ndofs-1,wdof);
        this.coefs.add(ndofs-1,coef);
        this.constraintValues = new double[constraintValues.length];
        System.arraycopy(constraintValues, 0, this.constraintValues, 0, constraintValues.length);
        this.SecDerconstraintValues = new double[SecDerconstraintValues.length];
        System.arraycopy(SecDerconstraintValues, 0, this.SecDerconstraintValues, 0, SecDerconstraintValues.length);
    }

    public void printC(){
        System.out.println("Id of consraint: "+this.id+" ..... to be fixed");
    }
    
    // methods
    public void setconstraintValues(double[] constraintValues){
        this.constraintValues = new double[constraintValues.length];
        System.arraycopy(constraintValues, 0, this.constraintValues, 0, constraintValues.length);
    }
    
    public void setSecDerconstraintValues(double[] SecDerconstraintValues){
        this.SecDerconstraintValues = new double[SecDerconstraintValues.length];
        System.arraycopy(SecDerconstraintValues, 0, this.SecDerconstraintValues, 0, SecDerconstraintValues.length);
    }
    
    public void setFirstDerconstraintValues(double[] FirstDerconstraintValues){
        this.FirstDerconstraintValues = new double[FirstDerconstraintValues.length];
        System.arraycopy(FirstDerconstraintValues, 0, this.FirstDerconstraintValues, 0, FirstDerconstraintValues.length);
    }
    
    public void addTerm(double coef, Node aNode, int wdof){
        ++this.ndofs;
        this.putNode(aNode);
        int[] dofs= new int[1];
        dofs[0]=wdof;
        aNode.setNdofs_ofNode(dofs);
        //this.NodeDOFS.put(this.numNodes, wdof);
        this.NodeDOFS.add(ndofs-1,wdof);
        this.coefs.add(ndofs-1,coef);
    }
    
    public void setAk(double val){
        this.ak=val;
    }
    
    public void setAm(double val){
        this.am=val;
    }
    
    public void setAc(double val){
        this.ac=val;
    }
    
    public AbstractMatrix getB(){
        MyMatrix = new AbstractMatrix(1,ndofs);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            int gdof;
            gdof=theNode.getFtable()[dof-1];
            double cf=coefs.get(j-1);
            MyMatrix.putVal(0, j-1, cf);
        }
        return MyMatrix;
    }
    

    @Override
    public AbstractMatrix getK() {
        MyMatrix=getB().transpose().times(getB());
        MyMatrix=MyMatrix.times(ak);
        return MyMatrix;
    }
    
    @Override
    public AbstractMatrix getDK(int param, int id) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        return MyMatrix;
    }
    
    @Override
    public AbstractMatrix getDM(int param, int id) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        return MyMatrix;
    }

    @Override
    public AbstractMatrix getM() {
        MyMatrix=getB().transpose().times(getB());
        MyMatrix=MyMatrix.times(am);
        return MyMatrix;
    }

    @Override
    int getElementNdofs() {
        return ndofs;
    }

    @Override
    int[] getFtable() {
        int[] theFtable = new int[this.ndofs];
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            int gdof;
            gdof=theNode.getFtable()[dof-1];
            theFtable[j-1]=gdof;
        }
        return theFtable;
    }

    @Override
    public AbstractMatrix getF() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            disps.addVal(j-1, 0, theNode.getDispsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getK().times(disps);
        return Fint;
    }

    @Override
    void clear() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    AbstractMatrix getM_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            velcs.addVal(j-1, 0, theNode.getVelcsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getM().times(velcs);
        return Fint;
    }

    @Override
    AbstractMatrix getM_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            accls.addVal(j-1, 0, theNode.getAcclsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getM().times(accls);
        return Fint;
    }

    @Override
    AbstractMatrix getM_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            disps.addVal(j-1, 0, theNode.getDispsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getM().times(disps);
        return Fint;
    }
    
    AbstractMatrix getK_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            disps.addVal(j-1, 0, theNode.getDispsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getK().times(disps);
        return Fint;
    }
    
    @Override
    AbstractMatrix getK_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            velcs.addVal(j-1, 0, theNode.getVelcsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getK().times(velcs);
        return Fint;
    }

    @Override
    AbstractMatrix getK_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            accls.addVal(j-1, 0, theNode.getAcclsTrial()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getK().times(accls);
        return Fint;
    }
    
    public AbstractMatrix getFloading_stif(int step){
        AbstractMatrix F = new AbstractMatrix(ndofs,1);
        double val;
        if(this.constraintValues==null){val=0.;}else{
            if(step+1>this.constraintValues.length){
                val=0.;
            }else{
                val=this.constraintValues[step];
            }
        }
        
        F=getB().transpose();
        F=F.times(ak*val);
        return F;
    }
    
    public AbstractMatrix getFloading_damp(int step){
        AbstractMatrix F = new AbstractMatrix(ndofs,1);
        double val;
        if(this.FirstDerconstraintValues==null){val=0.;}else{
            if(step+1>this.FirstDerconstraintValues.length){
                val=0.;
            }else{
                val=this.FirstDerconstraintValues[step];
            }
        }
        F=getB().transpose();
        F=F.times(ac*val);
        return F;
    }
    
    public AbstractMatrix getFloading_inert(int step){
        AbstractMatrix F = new AbstractMatrix(ndofs,1);
        double val;
        if(this.SecDerconstraintValues==null){val=0.;}else{
            if(step+1>this.SecDerconstraintValues.length){
                val=0.;
            }else{
                val=this.SecDerconstraintValues[step];
            }
        }
        F=getB().transpose();
        F=F.times(am*val);
        return F;
    }

    @Override
    void commit() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean ElementPlastified() {
        return false;
    }

    @Override
    public double getuKu() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            disps.addVal(j-1, 0, theNode.getDispsTrial()[dof-1]);
        }
        return (disps.transpose().times(this.getK().times(disps))).get(0, 0);
    }
    
    @Override
    public double getvMv() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            velcs.addVal(j-1, 0, theNode.getVelcsTrial()[dof-1]);
        }
        return (velcs.transpose().times(this.getM().times(velcs))).get(0, 0);
    }

    @Override
    public Node getNodeHierarchy(int nodeID) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getT() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        MyMatrix.init(0.0);
        for(int i=0;i<ndofs;i++){
            MyMatrix.set(i, i, 1.0);
        }
        return MyMatrix;
    }
    
    public double getVolume(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
     @Override
    public AbstractMatrix getM_upre() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            disps.addVal(j-1, 0, theNode.getDispsTrial()[dof-1]-theNode.getDeltaDisps()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getM().times(disps);
        return Fint;
    }
    
    @Override
    public AbstractMatrix getK_upre() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            ++j;
            int dof=NodeDOFS.get(j-1);
            disps.addVal(j-1, 0, theNode.getDispsTrial()[dof-1]-theNode.getDeltaDisps()[dof-1]);
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getK().times(disps);
        return Fint;
    }

}
