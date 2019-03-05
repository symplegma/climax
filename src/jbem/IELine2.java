/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class IELine2 extends IELine{


    // constructor
    public IELine2(InterfaceNode Node1, InterfaceNode Node2){
        ++numofIElems;
        this.id=numofIElems;
        this.theINodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putIElement(this);
        this.theINodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putIElement(this);
        this.numINodes=2;
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
    public double[] getCoordMinDistOfINode(InterfaceNode aNode) {
        double[] theXSI = new double[1];
        double x1,x2,x3,y1,y2,y3;
        x1=this.getINodeHier(1).getCoordinates()[0];
        y1=this.getINodeHier(1).getCoordinates()[1];
        x2=this.getINodeHier(2).getCoordinates()[0];
        y2=this.getINodeHier(2).getCoordinates()[1];
        x3=aNode.getCoordinates()[0];
        y3=aNode.getCoordinates()[1];
        theXSI[0]=( (x3-x1)*(x2-x1)+(y3-y1)*(y2-y1) )/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        theXSI[0] = 2.*theXSI[0]-1.;
        return theXSI;
    }
    
    @Override
    public double getLength() {
        double L=0.;
        double x1,x2,y1,y2;
        x1=this.getINodeHier(1).getIntermediateCoords()[0];
        y1=this.getINodeHier(1).getIntermediateCoords()[1];
        x2=this.getINodeHier(2).getIntermediateCoords()[0];
        y2=this.getINodeHier(2).getIntermediateCoords()[1];
        L=Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
        return L;
    }


}
