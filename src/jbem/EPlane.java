/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import geom.Point;

import java.util.Iterator;
import javax.swing.JOptionPane;

/**
 *
 * @author pchr
 */
abstract public class EPlane extends Element{
    // constructor
    public EPlane(){
        this.ContactElement=new boolean[3];
        this.NeumannElement=new boolean[3];
        this.DirichletElement=new boolean[3];
        for(int i=0;i<3;i++){
            this.ContactElement[i]=false;
            this.NeumannElement[i]=false;
            this.DirichletElement[i]=false;
        }
    }
    
    // abstract methods
    abstract public double getShapeFunction(int wSF, double xsi, double eta);
    abstract public double getShapeFunction_xsi(int wSF, double xsi, double eta);
    abstract public double getShapeFunction_eta(int wSF, double xsi, double eta);
    
    // methods
    public double[] getNormal(double xsi, double eta){
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        double[] coords;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            coords=theNode.getCoordinates();
            
            vxsi[0]+=getShapeFunction_xsi(i,xsi,eta)*coords[0];
            vxsi[1]+=getShapeFunction_xsi(i,xsi,eta)*coords[1];
            vxsi[2]+=getShapeFunction_xsi(i,xsi,eta)*coords[2];
            
            veta[0]+=getShapeFunction_eta(i,xsi,eta)*coords[0];
            veta[1]+=getShapeFunction_eta(i,xsi,eta)*coords[1];
            veta[2]+=getShapeFunction_eta(i,xsi,eta)*coords[2];
        }
        normal[0]=vxsi[1]*veta[2]-vxsi[2]*veta[1];
        normal[1]=vxsi[2]*veta[0]-vxsi[0]*veta[2];
        normal[2]=vxsi[0]*veta[1]-vxsi[1]*veta[0];
        double n=Math.sqrt(
                normal[0]*normal[0]+
                normal[1]*normal[1]+
                normal[2]*normal[2]
                );
        normal[0]/=n;
        normal[1]/=n;
        normal[2]/=n;
        return normal;
    }
    
    public double[] getDistance(Point aNode, double xsi, double eta){
        double[] dist=new double[3];
        double[] NodeCoords = new double[3];
        double[] pointCoords = new double[3];
        pointCoords[0]=0.; pointCoords[1]=0.; pointCoords[2]=0.;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            NodeCoords=theNode.getCoordinates();
            
            pointCoords[0]+=getShapeFunction(i,xsi,eta)*NodeCoords[0];
            pointCoords[1]+=getShapeFunction(i,xsi,eta)*NodeCoords[1];
            pointCoords[2]+=getShapeFunction(i,xsi,eta)*NodeCoords[2];
            
        }
        NodeCoords=aNode.getCoordinates();
        
        dist[0]=pointCoords[0]-NodeCoords[0];
        dist[1]=pointCoords[1]-NodeCoords[1];
        dist[2]=pointCoords[2]-NodeCoords[2];
        
        return dist;
    }
    
    public double getJacobian(double xsi, double eta){
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        double[] coords;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            coords=theNode.getCoordinates();
            
            vxsi[0]+=getShapeFunction_xsi(i,xsi,eta)*coords[0];
            vxsi[1]+=getShapeFunction_xsi(i,xsi,eta)*coords[1];
            vxsi[2]+=getShapeFunction_xsi(i,xsi,eta)*coords[2];
            
            veta[0]+=getShapeFunction_eta(i,xsi,eta)*coords[0];
            veta[1]+=getShapeFunction_eta(i,xsi,eta)*coords[1];
            veta[2]+=getShapeFunction_eta(i,xsi,eta)*coords[2];
        }
        normal[0]=vxsi[1]*veta[2]-vxsi[2]*veta[1];
        normal[1]=vxsi[2]*veta[0]-vxsi[0]*veta[2];
        normal[2]=vxsi[0]*veta[1]-vxsi[1]*veta[0];
        double Jac=0;
        //double[] vec=this.getNormal(xsi, eta);
        Jac=Math.sqrt(
                normal[0]*normal[0]+
                normal[1]*normal[1]+
                normal[2]*normal[2]
                );
        if(Jac<=0.){
            JOptionPane.showMessageDialog(null, "zero or negative jacobian for element with id: "+this.id,
                    "warning",JOptionPane.ERROR_MESSAGE);
        }
        return Jac;
    }
    
    public double getArea(){
        SpaceQuadIntegrator quadSI = (SpaceQuadIntegrator) this.theSIntegrator;
        return quadSI.IntegrateArea(this);
    }
    
    public double getAreaDeformed(int step, int state){
        SpaceQuadIntegrator quadSI = (SpaceQuadIntegrator) this.theSIntegrator;
        return quadSI.IntegrateAreaDeformed(this, step, state);
    }
    
    public double getMinL(){
        double minL=0.;
        double anL=0.;
        int frsID;
        Node frsNode;
        int secID;
        Node secNode;
        int indicator=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            frsNode = it.next();
            frsID=frsNode.getID();
            for(Iterator<Node> it2=this.theNodes.values().iterator(); it2.hasNext();){
                secNode = it2.next();
                secID=secNode.getID();
                if(frsID!=secID){
                    anL=frsNode.getDist(secNode);
                    if(indicator==0){minL=anL;}
                    if(anL<minL){minL=anL;}
                    indicator+=1;
                }
            }
        }
        return minL;
    }
    
    public double getMaxL(){
        double minL=0.;
        double anL=0.;
        int frsID;
        Node frsNode;
        int secID;
        Node secNode;
        int indicator=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            frsNode = it.next();
            frsID=frsNode.getID();
            for(Iterator<Node> it2=this.theNodes.values().iterator(); it2.hasNext();){
                secNode = it2.next();
                secID=secNode.getID();
                if(frsID!=secID){
                    anL=frsNode.getDist(secNode);
                    if(indicator==0){minL=anL;}
                    if(anL>minL){minL=anL;}
                    indicator+=1;
                }
            }
        }
        return minL;
    }

    @Override
    public double getWork(Domain theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWork(Domain theDomain, int wstep, int state) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getPower(Domain theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getTIPower(Domain theDomain, int wstep, double dt) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWork_hmg(Domain theDomain) {
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
    public double getSideInequality(Domain theDomain, int step, int bound) {
        throw new UnsupportedOperationException("Not supported yet.");
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
    public double[] getGlobalCoordinates(double[] local) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getLocalCoordinates(double[] local) {
        if(this.numNodes!=4){System.err.println("error in getLocalCoordinates of EPlane"); System.exit(this.id);}
        double[] vals = new double[1];
        double x1= this.getNodeHier(1).getCoordinates()[0];
        double y1= this.getNodeHier(1).getCoordinates()[1];
        double z1= this.getNodeHier(1).getCoordinates()[2];
        double x2= this.getNodeHier(2).getCoordinates()[0];
        double y2= this.getNodeHier(2).getCoordinates()[1];
        double z2= this.getNodeHier(2).getCoordinates()[2];
        double x3= this.getNodeHier(3).getCoordinates()[0];
        double y3= this.getNodeHier(3).getCoordinates()[1];
        double z3= this.getNodeHier(3).getCoordinates()[2];
        double x4= this.getNodeHier(4).getCoordinates()[0];
        double y4= this.getNodeHier(4).getCoordinates()[1];
        double z4= this.getNodeHier(4).getCoordinates()[3];

        vals[0]=Math.sqrt( ((local[0]-x1)*(local[0]-x1)+(local[1]-y1)*(local[1]-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) );
        vals[0]= 2.*vals[0]-1.;
        return vals;
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
    
    abstract public double getDeformedMinDistOfNode(Point aNode,int step,int wstate);
    

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
    public double[] getCoordMinDistOfNode(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
