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
abstract public class ELine extends Element{
    
    // constructor
    public ELine(){
        this.ContactElement=new boolean[2];
        this.NeumannElement=new boolean[2];
        this.DirichletElement=new boolean[2];
        for(int i=0;i<2;i++){
            this.ContactElement[i]=false;
            this.NeumannElement[i]=false;
            this.DirichletElement[i]=false;
        }
    }
    
    // abstract methods
    abstract public double getShapeFunction(int wSF, double xsi);
    abstract public double getShapeFunction_xsi(int wSF, double xsi);
    
    // methods
    public double[] getNormal(double xsi){
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        double[] coords;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            coords=theNode.getCoordinates();
            
            vxsi[0]+=getShapeFunction_xsi(i,xsi)*coords[0];
            vxsi[1]+=getShapeFunction_xsi(i,xsi)*coords[1];
            vxsi[2]+=getShapeFunction_xsi(i,xsi)*coords[2];
            
            
        }
        veta[0]=0.;
        veta[1]=0.;
        veta[2]=1.;
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

    public double[] getTangent(double xsi){
        double[] vxsi = new double[3];
        for(int i=0; i<3; i++)vxsi[i]=0.;
        double[] coords;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            coords=theNode.getCoordinates();

            vxsi[0]+=getShapeFunction_xsi(i,xsi)*coords[0];
            vxsi[1]+=getShapeFunction_xsi(i,xsi)*coords[1];
            vxsi[2]+=getShapeFunction_xsi(i,xsi)*coords[2];


        }
        double n=Math.sqrt(
                vxsi[0]*vxsi[0]+
                vxsi[1]*vxsi[1]+
                vxsi[2]*vxsi[2]
                );
        vxsi[0]/=n;
        vxsi[1]/=n;
        vxsi[2]/=n;
        return vxsi;
    }

    public double[] getNormal(int nid){
        double xsi = 0;
        switch(this.getHierOfNode(nid)){
            case 1: xsi=-1.; break;
            case 2: xsi= 1.; break;
            case 3: xsi= 0.; break;
            default: System.err.println("does not exist node with id= "+nid+" on element with id= "+id); System.exit(nid); break;
        }
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        double[] coords;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            coords=theNode.getCoordinates();

            vxsi[0]+=getShapeFunction_xsi(i,xsi)*coords[0];
            vxsi[1]+=getShapeFunction_xsi(i,xsi)*coords[1];
            vxsi[2]+=getShapeFunction_xsi(i,xsi)*coords[2];


        }
        veta[0]=0.;
        veta[1]=0.;
        veta[2]=1.;
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
    
    public double[] getDistance(Point aNode, double xsi){
        double[] dist=new double[3];
        double[] NodeCoords = new double[3];
        double[] pointCoords = new double[3];
        pointCoords[0]=0.; pointCoords[1]=0.; pointCoords[2]=0.;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            NodeCoords=theNode.getCoordinates();
            
            pointCoords[0]+=getShapeFunction(i,xsi)*NodeCoords[0];
            pointCoords[1]+=getShapeFunction(i,xsi)*NodeCoords[1];
            pointCoords[2]+=getShapeFunction(i,xsi)*NodeCoords[2];
            
        }
        NodeCoords=aNode.getCoordinates();
        
        dist[0]=pointCoords[0]-NodeCoords[0];
        dist[1]=pointCoords[1]-NodeCoords[1];
        dist[2]=pointCoords[2]-NodeCoords[2];
        
        return dist;
    }
    
    public double getJacobian(double xsi){
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        double[] coords;
        Node theNode;
        for(int i=1; i<=this.numNodes; i++){
            theNode=this.getNodeHier(i);
            coords=theNode.getCoordinates();
            
            vxsi[0]+=getShapeFunction_xsi(i,xsi)*coords[0];
            vxsi[1]+=getShapeFunction_xsi(i,xsi)*coords[1];
            vxsi[2]+=getShapeFunction_xsi(i,xsi)*coords[2];
            
            
        }
        veta[0]=0.;
        veta[1]=0.;
        veta[2]=1.;
        normal[0]=vxsi[1]*veta[2]-vxsi[2]*veta[1];
        normal[1]=vxsi[2]*veta[0]-vxsi[0]*veta[2];
        normal[2]=vxsi[0]*veta[1]-vxsi[1]*veta[0];
        double Jac=0;
        //double[] vec=this.getNormal(xsi);
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
    
    public double getLength(){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateLength(this);
    }

    public double getWork(Domain theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWork(this,theDomain, wstep);
    }
    
    public double getWork(Domain theDomain, int wstep, int state){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWork(this,theDomain, wstep, state);
    }
    
    public double getWorkX(Domain theDomain, int wstep, int state){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWorkX(this,theDomain, wstep, state);
    }
    
    public double getWorkY(Domain theDomain, int wstep, int state){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWorkY(this,theDomain, wstep, state);
    }
    
    public double getPower(Domain theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegratePower(this,theDomain, wstep);
    }
    
    @Override
    public double getTIPower(Domain theDomain, int wstep, double dt) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public double getWork_hmg(Domain theDomain){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWork_hmg(this,theDomain);
    }

    public double getWorkX(Domain theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWorkX(this,theDomain, wstep);
    }

    public double getWorkY(Domain theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWorkY(this,theDomain, wstep);
    }

    public double getGeneralisedForce(Domain theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateGeneralisedForce(this,theDomain, wstep, nodeID, dof);
    }
    
    public double getGeneralisedForce(Domain theDomain, int DispStep, int TracStep, int DispState, int TracState, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateGeneralisedForce(this,theDomain, DispStep, TracStep, DispState, TracState, nodeID, dof);
    }
    
    public double getGeneralisedForce(Domain theDomain, int wstep, int nodeID, int dof, int state){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateGeneralisedForce(this,theDomain, wstep, nodeID, dof, state);
    }
    
    public double getGeneralisedForce_hmg(Domain theDomain, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateGeneralisedForce_hmg(this,theDomain, nodeID, dof);
    }

    public double getShapeFunctionProduct(int iSF, int jSF, double xsi) {
        double N;
        N=getShapeFunction(iSF, xsi)*getShapeFunction(jSF, xsi);
        return N;
    }

    public double getShapeFunctionProduct(int iSF, int jSF, int mSF, double xsi) {
        double N;
        N=getShapeFunction(iSF, xsi)*getShapeFunction(jSF, xsi)*getShapeFunction(mSF, xsi);
        return N;
    }

    public double getShapeFunctionProduct(int iSF, int jSF, int mSF, int nSF, double xsi) {
        double N;
        N=getShapeFunction(iSF, xsi)*getShapeFunction(jSF, xsi)*getShapeFunction(mSF, xsi)*getShapeFunction(nSF, xsi);
        return N;
    }

    public double getSideInequality(Domain theDomain, int step, int bound){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateSideInequality(this, theDomain, step, bound);
    }

    @Override
    public double[] getGlobalCoordinates(double[] local) {
        double[] vals = new double[this.getNodeHier(1).getCoordinates().length];
        for(int i=0;i<this.getNodeHier(1).getCoordinates().length;i++){
            vals[i]=0.0;
            for(int j=1;j<=this.numNodes;j++){
                vals[i]+=this.getNodeHier(j).getCoordinates()[i]*this.getShapeFunction(j, local[0]);
            }
        }
        return vals;
    }

    @Override
    public double[] getLocalCoordinates(double[] local) {
        if(this.numNodes>2){System.err.println("error in getLocalCoordinates of ELine"); System.exit(this.id);}
        double[] vals = new double[1];
        double x1= this.getNodeHier(1).getCoordinates()[0];
        double y1= this.getNodeHier(1).getCoordinates()[1];
        double x2= this.getNodeHier(2).getCoordinates()[0];
        double y2= this.getNodeHier(2).getCoordinates()[1];

        vals[0]=Math.sqrt( ((local[0]-x1)*(local[0]-x1)+(local[1]-y1)*(local[1]-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) );
        vals[0]= 2.*vals[0]-1.;
        return vals;
    }
    
    @Override
    public double getWork(Domain theDomain, int DispStep, int TracStep, int DispState, int TracState) {
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateWork(this,theDomain, DispStep, TracStep, DispState, TracState);
    }
    
    @Override
    public double getTIPower(Domain theDomain, int wstep) {
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTIPower(this,theDomain, wstep);
    }
    
    @Override
    public double getTIPowerDifference(Domain theDomain, int wstep) {
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTIPowerDifference(this,theDomain, wstep);
    }
    
    public double getInternalPointDisp(Domain theDomain, ResultPoint theInternalPoint, int wdisp,int step,int wstate){
        double val=0.;
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
            val+=quadSI.IntegrateInternalPointDisp(this,theNode, theDomain, theInternalPoint, wdisp, step, wstate);
        }
        return val;
    }
    
    public double getInternalPointStress(Domain theDomain, ResultPoint theInternalPoint, int ws1,int step,int wstate){
        double val=0.;
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
//            if(this.id<=10 && step>0){
//                System.out.println("lower: "+id);
//            }
//            if(this.id<=120 && this.id>=111&& step>0){
//                System.out.println("upper: "+id);
//            }
            val+=quadSI.IntegrateInternalPointStress(this,theNode, theDomain, theInternalPoint, ws1, step, wstate);
        }
        return val;
    }
    
    public double getInternalPointStressAux(Domain theDomain, ResultPoint theInternalPoint, int ws1,int step,int wstate){
        double val=0.;
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
//            if(this.id<=10 && step>0){
//                System.out.println("lower: "+id);
//            }
//            if(this.id<=120 && this.id>=111&& step>0){
//                System.out.println("upper: "+id);
//            }
            val+=quadSI.IntegrateInternalPointStressAux(this,theNode, theDomain, theInternalPoint, ws1, step, wstate);
        }
        return val;
    }
    
    @Override
    public double[] getDispLocalonNode(int nodeHier,int step,int wstate) {
        double[] dispsLocal;
        dispsLocal=new double[this.getNodeHier(nodeHier).getuEFTable().length];
        dispsLocal[0]=this.getNodeHier(nodeHier).getu()[0][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[0]
                +this.getNodeHier(nodeHier).getu()[1][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[1];
        dispsLocal[1]=-this.getNodeHier(nodeHier).getu()[0][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[1]
                +this.getNodeHier(nodeHier).getu()[1][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[0];
        return dispsLocal;
    }
    
    @Override
    public double[] getTractionLocalonNode(int nodeHier,int step,int wstate) {
        double[] dispsLocal;
        dispsLocal=new double[this.getNodeHier(nodeHier).getuEFTable().length];
        dispsLocal[0]=this.getNodeHier(nodeHier).getp(this, 2, wstate)[0][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[0]
                +this.getNodeHier(nodeHier).getp(this, 2, wstate)[1][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[1];
        dispsLocal[1]=-this.getNodeHier(nodeHier).getp(this, 2, wstate)[0][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[1]
                +this.getNodeHier(nodeHier).getp(this, 2, wstate)[1][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[0];
        return dispsLocal;
    }
    
    @Override
    public double[] getDispLocalonNodeAux(int nodeHier,int step,int wstate) {
        double[] dispsLocal;
        dispsLocal=new double[this.getNodeHier(nodeHier).getuEFTable().length];
        dispsLocal[0]=this.getNodeHier(nodeHier).getu_aux()[0][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[0]
                +this.getNodeHier(nodeHier).getu_aux()[1][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[1];
        dispsLocal[1]=-this.getNodeHier(nodeHier).getu_aux()[0][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[1]
                +this.getNodeHier(nodeHier).getu_aux()[1][step][wstate]*this.getNormal(this.getNodeHier(nodeHier).getID())[0];
        return dispsLocal;
    }
    
    @Override
    public double[] getTractionLocalonNodeAux(int nodeHier,int step,int wstate) {
        double[] dispsLocal;
        dispsLocal=new double[this.getNodeHier(nodeHier).getuEFTable().length];
        dispsLocal[0]=this.getNodeHier(nodeHier).getp_aux(this, 2, wstate)[0][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[0]
                +this.getNodeHier(nodeHier).getp_aux(this, 2, wstate)[1][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[1];
        dispsLocal[1]=-this.getNodeHier(nodeHier).getp_aux(this, 2, wstate)[0][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[1]
                +this.getNodeHier(nodeHier).getp_aux(this, 2, wstate)[1][step]*this.getNormal(this.getNodeHier(nodeHier).getID())[0];
        return dispsLocal;
    }
    
    public double getShapeFunctionProduct(int dof, int nodeID_i, int nodeID_j){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateShapeFunctionProduct(this,dof, nodeID_i, nodeID_j);
    }
    
    @Override
    public boolean PointOnElement(Point aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
