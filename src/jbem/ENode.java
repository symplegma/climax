/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import geom.Point;

/**
 *
 * @author pchr
 */
public final class ENode extends Element{
    private static int numENodes=0;
    
    // constructor
    public ENode(int id, Node Node1){
        ++numENodes;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1); this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.numNodes=1;
        this.setElemType();
        this.ContactElement=new boolean[1];
        this.NeumannElement=new boolean[1];
        this.DirichletElement=new boolean[1];
        this.DirichletElement[0]=false;
    }
    
    public ENode(int id, Node Node1, SpaceIntegrator theSI){
        ++numENodes;
        this.id=id;
        this.numNodes=1;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theSIntegrator=theSI;
        this.setElemType();
    }

    @Override
    public Node getNodePrint(int seq) {
        int id_=theNodesHierarchy.get(1);
        return this.theNodes.get(id_);
    }

    @Override
    public double getWork(Domain theDomain, int wstep) {
        SpaceNodeIntegrator quadSI = (SpaceNodeIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWork(this,theDomain, wstep);
    }
    
    @Override
    public double getWork(Domain theDomain, int wstep, int state) {
        SpaceNodeIntegrator quadSI = (SpaceNodeIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWork(this,theDomain, wstep,state);
    }
    
    @Override
    public double getPower(Domain theDomain, int wstep) {
        SpaceNodeIntegrator quadSI = (SpaceNodeIntegrator) this.theSIntegrator;
        return quadSI.IntegratePower(this,theDomain, wstep);
    }
    
    @Override
    public double getTIPower(Domain theDomain, int wstep, double dt) {
        SpaceNodeIntegrator quadSI = (SpaceNodeIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTIPower(this,theDomain, wstep,dt);
    }
    
    @Override
    public double getWork_hmg(Domain theDomain) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getSideInequality(Domain theDomain, int step, int bound) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected void setElemType() {
        this.ElemType=0;
    }

    @Override
    public double[] getCoordMinDistOfNode(Point aNode) {
        double[] theXSI = new double[1];
        theXSI[0]=0.;
        return theXSI;
    }

    @Override
    public double getMinDistOfNode(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedForce(Domain theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getGeneralisedForce(Domain theDomain, int wstep, int nodeID, int dof, int state) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getGeneralisedForce_hmg(Domain theDomain, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkX(Domain theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkY(Domain theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getGlobalCoordinates(double[] local) {
        double[] vals = new double[this.getNodeHier(1).getCoordinates().length];
        for(int i=0;i<this.getNodeHier(1).getCoordinates().length;i++){vals[i]=this.getNodeHier(1).getCoordinates()[i];}
        return vals;
    }

    @Override
    public double[] getLocalCoordinates(double[] global) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getCoordMinDistOfNode(InterfaceNode aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWork(Domain theDomain, int DispStep, int TracStep, int DispState, int TracState) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedForce(Domain theDomain, int DispStep, int TracStep, int DispState, int TracState, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getTIPower(Domain theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public double getInternalPointDisp(Domain theDomain, ResultPoint theInternalPoint, int wdisp,int step,int wstate){
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getInternalPointStress(Domain theDomain, ResultPoint theInternalPoint, int ws1, int step, int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getTIPowerDifference(Domain theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getDispLocalonNode(int nodeHier,int step,int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double[] getTractionLocalonNode(int nodeHier,int step,int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getDispLocalonNodeAux(int nodeHier, int step, int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getTractionLocalonNodeAux(int nodeHier, int step, int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getInternalPointStressAux(Domain theDomain, ResultPoint theInternalPoint, int wstress, int step, int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getShapeFunctionProduct(int dof, int nodeID_i, int nodeID_j) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public boolean PointOnElement(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
