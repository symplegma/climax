package jbem;

import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
abstract public class SpaceIntegrator {
    protected static int numofGauss;
    protected int xsi_mesh;
    protected int eta_mesh;
    protected static GaussData theGaussData;
    protected boolean DomainAuxiliaryField=false;
    
    // constructor
    public SpaceIntegrator(){
    }
    
    // methods
    public void setNumofGauss(int num){
        SpaceIntegrator.numofGauss=num;
        SpaceIntegrator.theGaussData= new GaussData(num);
        SpaceIntegrator.theGaussData.computeGaussPoints();
    }
    
    public void setNumofGauss(int num, int typeOfGauss){
        SpaceIntegrator.numofGauss=num;
        SpaceIntegrator.theGaussData= new GaussData(num);
        SpaceIntegrator.theGaussData.setGaussType(typeOfGauss);
        SpaceIntegrator.theGaussData.computeGaussPoints();
    }
    
    public void setMesh(int xmesh, int emesh){
        this.xsi_mesh=xmesh;
        this.eta_mesh=emesh;
    }
    
    public void setMesh(int xmesh){
        this.xsi_mesh=xmesh;
    }
    
    public int getXsi_mesh(){
        return this.xsi_mesh;
    }
    
    public int getEta_mesh(){
        return this.eta_mesh;
    }
    
    public int getNumofGauss(){
        return SpaceIntegrator.numofGauss;
    }
    
    public void setAuxiliaryField(boolean bool){this.DomainAuxiliaryField=bool;}
    
    public boolean getAuxiliaryField(){return this.DomainAuxiliaryField;}
    
    abstract AbstractMatrix IntegrateDu(Node colNode, Domain theDomain, Element elem, Node elemNode);
    abstract AbstractMatrix IntegrateDp(Node colNode, Domain theDomain, Element elem, Node elemNode);
    abstract double IntegrateM(Node colNode, Domain theDomain, Element elem, Node elemNode);
    abstract AbstractMatrix IntegrateDv(Node colNode, Domain theDomain, Element elem, Node elemNode);
    abstract AbstractMatrix IntegrateDp_res(Node colNode, Domain theDomain, Element elem, Node elemNode);
    abstract AbstractMatrix IntegrateDp_dif(Node colNode, Domain theDomain, Element elem, Node elemNode);

}
