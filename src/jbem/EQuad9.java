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
public final class EQuad9 extends EQuad{
    private static int numEQuad9s=0;
    
    // constructor
    public EQuad9(int id, Node Node1, Node Node2, Node Node3, Node Node4
                        , Node Node5, Node Node6, Node Node7, Node Node8
                        , Node Node9){
        ++numEQuad9s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theNodes.put(Node3.getID(), Node3);this.theNodesHierarchy.put(3, Node3.getID());
        Node3.putElement(this);
        this.theNodes.put(Node4.getID(), Node4);this.theNodesHierarchy.put(4, Node4.getID());
        Node4.putElement(this);
        this.theNodes.put(Node5.getID(), Node5);this.theNodesHierarchy.put(5, Node5.getID());
        Node5.putElement(this);
        this.theNodes.put(Node6.getID(), Node6);this.theNodesHierarchy.put(6, Node6.getID());
        Node6.putElement(this);
        this.theNodes.put(Node7.getID(), Node7);this.theNodesHierarchy.put(7, Node7.getID());
        Node7.putElement(this);
        this.theNodes.put(Node8.getID(), Node8);this.theNodesHierarchy.put(8, Node8.getID());
        Node8.putElement(this);
        this.theNodes.put(Node9.getID(), Node9);this.theNodesHierarchy.put(9, Node9.getID());
        Node9.putElement(this);
        this.numNodes=9;
        this.setElemType();
    }
    
    public EQuad9(int id, Node Node1, Node Node2, Node Node3, Node Node4
                        , Node Node5, Node Node6, Node Node7, Node Node8
                        , Node Node9, SpaceIntegrator theSI){
        ++numEQuad9s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theNodes.put(Node3.getID(), Node3);this.theNodesHierarchy.put(3, Node3.getID());
        Node3.putElement(this);
        this.theNodes.put(Node4.getID(), Node4);this.theNodesHierarchy.put(4, Node4.getID());
        Node4.putElement(this);
        this.theNodes.put(Node5.getID(), Node5);this.theNodesHierarchy.put(5, Node5.getID());
        Node5.putElement(this);
        this.theNodes.put(Node6.getID(), Node6);this.theNodesHierarchy.put(6, Node6.getID());
        Node6.putElement(this);
        this.theNodes.put(Node7.getID(), Node7);this.theNodesHierarchy.put(7, Node7.getID());
        Node7.putElement(this);
        this.theNodes.put(Node8.getID(), Node8);this.theNodesHierarchy.put(8, Node8.getID());
        Node8.putElement(this);
        this.theNodes.put(Node9.getID(), Node9);this.theNodesHierarchy.put(9, Node9.getID());
        Node9.putElement(this);
        this.theSIntegrator=theSI;
        this.numNodes=9;
        this.setElemType();
    }
    
    /**
     * returns the wSF (wSF may be 1 to 9)
     * at the point with local coordinates
     * xsi,eta
     */
    public double getShapeFunction(int wSF, double xsi, double eta){
        double N;
        switch(wSF){
            case 1: N=0.25*(1-xsi)*(1-eta)
                    -0.5*this.getShapeFunction(5, xsi, eta)
                    -0.5*this.getShapeFunction(8, xsi, eta)
                    -0.25*this.getShapeFunction(9, xsi, eta);   
            break;
            case 2: N=0.25*(1+xsi)*(1-eta)
                    -0.5*this.getShapeFunction(5, xsi, eta)
                    -0.5*this.getShapeFunction(6, xsi, eta)
                    -0.25*this.getShapeFunction(9, xsi, eta);  
            break; 
            case 3: N=0.25*(1+xsi)*(1+eta)
                    -0.5*this.getShapeFunction(6, xsi, eta)
                    -0.5*this.getShapeFunction(7, xsi, eta)
                    -0.25*this.getShapeFunction(9, xsi, eta); 
            break;
            case 4: N=0.25*(1-xsi)*(1+eta)
                    -0.5*this.getShapeFunction(7, xsi, eta)
                    -0.5*this.getShapeFunction(8, xsi, eta)
                    -0.25*this.getShapeFunction(9, xsi, eta); 
            break;
            case 5: N=0.5*(1-xsi*xsi)*(1-eta)
                    -0.5*this.getShapeFunction(9, xsi, eta); 
            break;
            case 6: N=0.5*(1+xsi)*(1-eta*eta)
                    -0.5*this.getShapeFunction(9, xsi, eta); 
            break;
            case 7: N=0.5*(1-xsi*xsi)*(1+eta)
                    -0.5*this.getShapeFunction(9, xsi, eta); 
            break;
            case 8: N=0.5*(1-xsi)*(1-eta*eta)
                    -0.5*this.getShapeFunction(9, xsi, eta); 
            break;
            case 9: N=(1-xsi*xsi)*(1-eta*eta); 
            break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    public double getShapeFunction_xsi(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=-0.25*(1-eta)
                    -0.5*this.getShapeFunction_xsi(5, xsi, eta)
                    -0.5*this.getShapeFunction_xsi(8, xsi, eta)
                    -0.25*this.getShapeFunction_xsi(9, xsi, eta);  
            break;
            case 2: N=0.25*(1-eta)
                    -0.5*this.getShapeFunction_xsi(5, xsi, eta)
                    -0.5*this.getShapeFunction_xsi(6, xsi, eta)
                    -0.25*this.getShapeFunction_xsi(9, xsi, eta); 
            break; 
            case 3: N=0.25*(1+eta)
                    -0.5*this.getShapeFunction_xsi(6, xsi, eta)
                    -0.5*this.getShapeFunction_xsi(7, xsi, eta)
                    -0.25*this.getShapeFunction_xsi(9, xsi, eta); 
            break;
            case 4: N=-0.25*(1+eta)
                    -0.5*this.getShapeFunction_xsi(7, xsi, eta)
                    -0.5*this.getShapeFunction_xsi(8, xsi, eta)
                    -0.25*this.getShapeFunction_xsi(9, xsi, eta); 
            break;
            case 5: N=0.5*(-2.*xsi)*(1-eta)
                    -0.5*this.getShapeFunction_xsi(9, xsi, eta); 
            break;
            case 6: N=0.5*(1-eta*eta)
                    -0.5*this.getShapeFunction_xsi(9, xsi, eta); 
            break;
            case 7: N=0.5*(-2.*xsi)*(1+eta)
                    -0.5*this.getShapeFunction_xsi(9, xsi, eta); 
            break;
            case 8: N=0.5*(-1.)*(1-eta*eta)
                    -0.5*this.getShapeFunction_xsi(9, xsi, eta);  
            break;
            case 9: N=(-2.*xsi)*(1-eta*eta); 
            break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    public double getShapeFunction_eta(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=-0.25*(1-xsi)
                    -0.5*this.getShapeFunction_eta(5, xsi, eta)
                    -0.5*this.getShapeFunction_eta(8, xsi, eta)
                    -0.25*this.getShapeFunction_eta(9, xsi, eta);  
            break;
            case 2: N=-0.25*(1+xsi)
                    -0.5*this.getShapeFunction_eta(5, xsi, eta)
                    -0.5*this.getShapeFunction_eta(6, xsi, eta)
                    -0.25*this.getShapeFunction_eta(9, xsi, eta); 
            break; 
            case 3: N=0.25*(1+xsi)
                    -0.5*this.getShapeFunction_eta(6, xsi, eta)
                    -0.5*this.getShapeFunction_eta(7, xsi, eta)
                    -0.25*this.getShapeFunction_eta(9, xsi, eta); 
            break;
            case 4: N=0.25*(1-xsi)
                    -0.5*this.getShapeFunction_eta(7, xsi, eta)
                    -0.5*this.getShapeFunction_eta(8, xsi, eta)
                    -0.25*this.getShapeFunction_eta(9, xsi, eta); 
            break;
            case 5: N=0.5*(1-xsi*xsi)*(-1.)
                    -0.5*this.getShapeFunction_eta(9, xsi, eta); 
            break;
            case 6: N=0.5*(1+xsi)*(-2*eta)
                    -0.5*this.getShapeFunction_eta(9, xsi, eta);  
            break;
            case 7: N=0.5*(1-xsi*xsi)
                    -0.5*this.getShapeFunction_eta(9, xsi, eta); 
            break;
            case 8: N=0.5*(1-xsi)*(-2.*eta)
                    -0.5*this.getShapeFunction_eta(9, xsi, eta); 
            break;
            case 9: N=(1-xsi*xsi)*(-2.*eta); 
            break;
            default:N=0. ;  break;
        }
        return N;
    }

    @Override
    public Node getNodePrint(int seq) {
        int seq2=0;
        switch(seq){
            case 1:seq2=1;break;
            case 2:seq2=5;break;
            case 3:seq2=2;break;
            case 4:seq2=6;break;
            case 5:seq2=3;break;
            case 6:seq2=7;break;
            case 7:seq2=4;break;
            case 8:seq2=8;break;
            case 9:seq2=9;break;
        }
        int id_=theNodesHierarchy.get(seq2);
        return this.theNodes.get(id_);
    }

    @Override
    protected void setElemType() {
        this.ElemType=5;
    }

    @Override
    public double getDeformedMinDistOfNode(Point aNode, int step, int wstate) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean PointOnElement(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
