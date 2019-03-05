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
public final class ELine3 extends ELine{
    private static int numELine3s=0;
    
    // constructor
    public ELine3(int id, Node Node1, Node Node2, Node Node3){
        ++numELine3s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theNodes.put(Node3.getID(), Node3);this.theNodesHierarchy.put(3, Node3.getID());
        Node3.putElement(this);
        this.numNodes=3;
        this.setElemType();
    }
    
    public ELine3(int id, Node Node1, Node Node2, Node Node3, SpaceIntegrator theSI){
        ++numELine3s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theNodes.put(Node3.getID(), Node3);this.theNodesHierarchy.put(3, Node3.getID());
        Node3.putElement(this);
        this.theSIntegrator=theSI;
        this.numNodes=3;
        this.setElemType();
    }

    @Override
    public double getShapeFunction(int wSF, double xsi) {
        double N;
        switch(wSF){
            case 1: N=0.5*(1-xsi)-0.5*(1-xsi*xsi); break;
            case 2: N=0.5*(1+xsi)-0.5*(1-xsi*xsi); break;
            case 3: N=1-xsi*xsi; break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    public double getShapeFunction_xsi(int wSF, double xsi) {
        double N;
        switch(wSF){
            case 1: N=-0.5+xsi; break;
            case 2: N= 0.5+xsi; break;
            default:N= -2.*xsi ;  break;
        }
        return N;
    }

    @Override
    public Node getNodePrint(int seq) {
        int seq2=0;
        switch(seq){
            case 1:seq2=1;break;
            case 2:seq2=3;break;
            case 3:seq2=2;break;
        }
        int id_=theNodesHierarchy.get(seq2);
        return this.theNodes.get(id_);
    }

    @Override
    protected void setElemType() {
        this.ElemType=2;
    }

    @Override
    public double[] getCoordMinDistOfNode(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getMinDistOfNode(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getCoordMinDistOfNode(InterfaceNode aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
