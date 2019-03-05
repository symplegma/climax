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
public class ELine2 extends ELine{
    private static int numELine2s=0;
    
    // constructor
    public ELine2(){}
    
    public ELine2(int id, Node Node1, Node Node2){
        ++numELine2s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.numNodes=2;
        this.setElemType();
    }
    
    public ELine2(int id, Node Node1, Node Node2, SpaceIntegrator theSI){
        ++numELine2s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theSIntegrator=theSI;
        this.numNodes=2;
        this.setElemType();
    }

    @Override
    public double getShapeFunction(int wSF, double xsi) {
        double N;
        switch(wSF){
            case 1: N=0.5*(1-xsi); break;
            case 2: N=0.5*(1+xsi); break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    public double getShapeFunction_xsi(int wSF, double xsi) {
        double N;
        switch(wSF){
            case 1: N=-0.5; break;
            case 2: N=0.5; break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    public Node getNodePrint(int seq) {
        int seq2=0;
        switch(seq){
            case 1:seq2=1;break;
            case 2:seq2=2;break;
        }
        int id_=theNodesHierarchy.get(seq2);
        return this.theNodes.get(id_);
    }

    @Override
    protected void setElemType() {
        this.ElemType=1;
    }

    @Override
    public double[] getCoordMinDistOfNode(Point aNode) {
        double[] theXSI = new double[1];
        double x1,x2,x3,y1,y2,y3;
        x1=this.getNodeHier(1).getCoordinates()[0];
        y1=this.getNodeHier(1).getCoordinates()[1];
        x2=this.getNodeHier(2).getCoordinates()[0];
        y2=this.getNodeHier(2).getCoordinates()[1];
        x3=aNode.getCoordinates()[0];
        y3=aNode.getCoordinates()[1];
        theXSI[0]=((x3-x1)*(x2-x1)+(y3-y1)*(y2-y1) )/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        theXSI[0] = 2.*theXSI[0]-1.;
        return theXSI;
    }

    @Override
    public double[] getCoordMinDistOfNode(InterfaceNode aNode) {
        double[] theXSI = new double[1];
        double x1,x2,x3,y1,y2,y3;
        x1=this.getNodeHier(1).getCoordinates()[0];
        y1=this.getNodeHier(1).getCoordinates()[1];
        x2=this.getNodeHier(2).getCoordinates()[0];
        y2=this.getNodeHier(2).getCoordinates()[1];
        x3=aNode.getCoordinates()[0];
        y3=aNode.getCoordinates()[1];
        theXSI[0]=( (x3-x1)*(x2-x1)+(y3-y1)*(y2-y1) )/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        theXSI[0] = 2.*theXSI[0]-1.;
        return theXSI;
    }

    @Override
    public double getMinDistOfNode(Point aNode) {
        double[] theXSI = new double[1];
        double x1,x2,x3,y1,y2,y3,dist = 0;
        x1=this.getNodeHier(1).getCoordinates()[0];
        y1=this.getNodeHier(1).getCoordinates()[1];
        x2=this.getNodeHier(2).getCoordinates()[0];
        y2=this.getNodeHier(2).getCoordinates()[1];
        x3=aNode.getCoordinates()[0];
        y3=aNode.getCoordinates()[1];
        theXSI[0]=( (x3-x1)*(x2-x1)+(y3-y1)*(y2-y1) )/( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) );
        theXSI[0] = 2.*theXSI[0]-1.;
        x1=this.getShapeFunction(1, theXSI[0])*this.getNodeHier(1).getCoordinates()[0]+
                this.getShapeFunction(2, theXSI[0])*this.getNodeHier(2).getCoordinates()[0];
        y1=this.getShapeFunction(1, theXSI[0])*this.getNodeHier(1).getCoordinates()[1]+
                this.getShapeFunction(2, theXSI[0])*this.getNodeHier(2).getCoordinates()[1];
        dist= Math.sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
        return dist;
    }

}
