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
public final class EQuad4 extends EQuad{
    private static int numEQuad4s=0;
    private static Node tNode1,tNode2,tNode3,tNode4;
    
    // constructor
    public EQuad4(int id, Node Node1, Node Node2, Node Node3, Node Node4){
        ++numEQuad4s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theNodes.put(Node3.getID(), Node3);this.theNodesHierarchy.put(3, Node3.getID());
        Node3.putElement(this);
        this.theNodes.put(Node4.getID(), Node4);this.theNodesHierarchy.put(4, Node4.getID());
        Node4.putElement(this);
        this.numNodes=4;
        this.setElemType();
    }
    
    public EQuad4(int id, Node Node1, Node Node2, Node Node3, Node Node4, SpaceIntegrator theSI){
        ++numEQuad4s;
        this.id=id;
        this.theNodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putElement(this);
        this.theNodes.put(Node2.getID(), Node2);this.theNodesHierarchy.put(2, Node2.getID());
        Node2.putElement(this);
        this.theNodes.put(Node3.getID(), Node3);this.theNodesHierarchy.put(3, Node3.getID());
        Node3.putElement(this);
        this.theNodes.put(Node4.getID(), Node4);this.theNodesHierarchy.put(4, Node4.getID());
        Node4.putElement(this);
        this.theSIntegrator=theSI;
        this.numNodes=4;
        this.setElemType();
    }
    
     /**
     * returns the wSF (wSF may be 1,2,3,4)
     * at the point with local coordinates
     * xsi,eta
     */
    public double getShapeFunction(int wSF, double xsi, double eta){
        double N;
        switch(wSF){
            case 1: N=0.25*(1-xsi)*(1-eta); break;
            case 2: N=0.25*(1+xsi)*(1-eta); break;
            case 3: N=0.25*(1+xsi)*(1+eta); break;
            case 4: N=0.25*(1-xsi)*(1+eta); break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    public double getShapeFunction_xsi(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=-0.25*(1-eta); break;
            case 2: N= 0.25*(1-eta); break;
            case 3: N= 0.25*(1+eta); break;
            case 4: N=-0.25*(1+eta); break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    public double getShapeFunction_eta(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=-0.25*(1-xsi); break;
            case 2: N=-0.25*(1+xsi); break;
            case 3: N= 0.25*(1+xsi); break;
            case 4: N= 0.25*(1-xsi); break;
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
            case 3:seq2=3;break;
            case 4:seq2=4;break;
        }
        int id_=theNodesHierarchy.get(seq2);
        return this.theNodes.get(id_);
    }

    @Override
    protected void setElemType() {
        this.ElemType=3;
    }
    
    @Override
    public double getMinDistOfNode(Point aNode) {
        double Dist=0.;
        double[] normal=this.getNormal(-1., -1.);
        Dist=Math.abs(normal[0]*(this.getNodeHier(1).getCoordinates()[0]-aNode.getCoordinates()[0])+
                normal[1]*(this.getNodeHier(1).getCoordinates()[1]-aNode.getCoordinates()[1])+
                normal[2]*(this.getNodeHier(1).getCoordinates()[2]-aNode.getCoordinates()[2]))/
                (Math.sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]));
        return Dist;
    }
    
    public double getDeformedMinDistOfNode(Point aNode, int step, int state) {
        double[] coords = new double[3];
        
        // fix tNode1
        coords[0]=this.getNodeHier(1).getCoordinates()[0]+this.getNodeHier(1).getu()[0][step][state];
        coords[1]=this.getNodeHier(1).getCoordinates()[1]+this.getNodeHier(1).getu()[1][step][state];
        coords[2]=this.getNodeHier(1).getCoordinates()[2]+this.getNodeHier(1).getu()[2][step][state];
        tNode1=new Node(1,coords);
        // fix tNode2
        coords[0]=this.getNodeHier(2).getCoordinates()[0]+this.getNodeHier(2).getu()[0][step][state];
        coords[1]=this.getNodeHier(2).getCoordinates()[1]+this.getNodeHier(2).getu()[1][step][state];
        coords[2]=this.getNodeHier(2).getCoordinates()[2]+this.getNodeHier(2).getu()[2][step][state];
        tNode2=new Node(2,coords);
        // fix tNode3
        coords[0]=this.getNodeHier(3).getCoordinates()[0]+this.getNodeHier(3).getu()[0][step][state];
        coords[1]=this.getNodeHier(3).getCoordinates()[1]+this.getNodeHier(3).getu()[1][step][state];
        coords[2]=this.getNodeHier(3).getCoordinates()[2]+this.getNodeHier(3).getu()[2][step][state];
        tNode3=new Node(3,coords);
        // fix tNode4
        coords[0]=this.getNodeHier(4).getCoordinates()[0]+this.getNodeHier(4).getu()[0][step][state];
        coords[1]=this.getNodeHier(4).getCoordinates()[1]+this.getNodeHier(4).getu()[1][step][state];
        coords[2]=this.getNodeHier(4).getCoordinates()[2]+this.getNodeHier(4).getu()[2][step][state];
        tNode4=new Node(4,coords);
        EQuad4 tempQuad=new EQuad4(1,tNode1,tNode2,tNode3,tNode4);
        
        return tempQuad.getMinDistOfNode(aNode);
    }
    
    public boolean PointOnElement(Point aNode) {
        boolean bool=false;
        double[] coordsOnPlane = new double[3];
        double a,b,c,d,e,f,x,y,z,t;
        double[] n=this.getNormal(this.getNodeHier(1).getID());
        a=n[0]; b=n[1]; c=n[2];
        
        d=this.getNodeHier(1).getCoordinates()[0];
        e=this.getNodeHier(1).getCoordinates()[1];
        f=this.getNodeHier(1).getCoordinates()[2];
        
        x=aNode.getCoordinates()[0];
        y=aNode.getCoordinates()[1];
        z=aNode.getCoordinates()[0];
        
        t=(a*d-a*x+b*e-b*y+c*f-c*z)/(a*a+b*b+c*c);
        
        coordsOnPlane[0]=x+t*a;
        coordsOnPlane[1]=y+t*b;
        coordsOnPlane[2]=z+t*c;
        // now coordsOnPlane are the coordinates of the point which is the 
        // projection of the aNode onto the plane of EQuad4 ('this')
        if(PointInTriangle(coordsOnPlane,this.getNodeHier(1).getCoordinates(),this.getNodeHier(2).getCoordinates(),this.getNodeHier(4).getCoordinates()) || 
                PointInTriangle(coordsOnPlane,this.getNodeHier(2).getCoordinates(),this.getNodeHier(3).getCoordinates(),this.getNodeHier(4).getCoordinates()) )bool=true;
        if(bool && (getMinDistOfNode(aNode)<=1.e10))bool=true;
        return bool;
    }
    
    public boolean PointInTriangle(double[] point, double[] e1, double[] e2, double[] e3){
        boolean answer=false;
        if( (area(point,e1,e2)+area(point,e2,e3)+area(point,e3,e1)-area(e1,e2,e3))<=1.e-10)answer=true;
        return answer;
    }
    
    private double area(double[] p1, double[] p2, double[] p3){
        //http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
        double a=0.;
        double x1=p1[0]-p2[0];
        double x2=p1[1]-p2[1];
        double x3=p1[2]-p2[2];
        double y1=p3[0]-p2[0];
        double y2=p3[1]-p2[1];
        double y3=p3[2]-p2[2];
        a=0.5*Math.sqrt((x2*y3-x3*y2)*(x2*y3-x3*y2)+
                        (x3*y1-x1*y3)*(x3*y1-x1*y3)+
                        (x1*y2-x2*y1)*(x1*y2-x2*y1));
        return a;
    }
}
