/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import jbem.Domain.BoundaryType;

/**
 *
 * @author pchr
 */
abstract public class Element {
    protected int id;
    protected Map<Integer,Node> theNodes = new HashMap<Integer,Node>();
    protected Map<Integer,Integer> theNodesHierarchy = new HashMap<Integer,Integer>();
    protected int numNodes;
    protected static SpaceIntegrator theSIntegrator;
    protected boolean[] NeumannElement;
    protected boolean[] DirichletElement;
    protected boolean[] ContactElement;
    protected boolean MainElement=false;    
    protected int ElemType;

    // constructor
    public Element(){}
    
    // methods
    public int getID(){
        return this.id;
    }
    
    public Map getNodes(){
        return this.theNodes;
    }
    
    /**
     * 
     * @param id_ the id of the node which returned 
     * @return the node with id=id_
     */
    public Node getNode(int id_){
        return this.theNodes.get(id_);
    }
    
    /**
     * 
     * @param seq the hierarchy of the node returned
     * @return the Node with hierarchy seq
     */
    public Node getNodeHier(int seq){
        int id_=theNodesHierarchy.get(seq);
        return this.theNodes.get(id_);
    }
    
    /**
     * 
     * @param theNode_id is the id of the node for which the hierarchy is returned
     * @return the hierarchy of the Node with id: theNode_id
     * if no such node belongs to the element a zero value returned
     */
    public int getHierOfNode(int theNode_id){
        int seq=0;
        int i=0;
        for(Iterator<Integer> it=this.theNodesHierarchy.values().iterator(); it.hasNext();){
            ++i;
            if(it.next()==theNode_id){seq=i;}
        }
        return seq;
    }
    
    public int getNumNodes(){
        return this.numNodes;
    }
    
    public SpaceIntegrator getSpaceIntegrator(){
        return theSIntegrator;
    }
    
    public void setSpaceIntegrator(SpaceIntegrator SI){
        theSIntegrator=SI;
    }
    
    public int[] getpEFTable(Domain theDomain){
        
        int ndof=theDomain.getFundamentalSolution().get_p_DOFs();
        int[] EFT = new int[ndof*this.numNodes];
        int[] eft ;//= new int[ndof];
        int j=0;
        for(Iterator<Integer> it=this.theNodesHierarchy.values().iterator(); it.hasNext();){
            int hier = it.next();
            Node theNode;
            theNode=this.theNodes.get(hier);
            eft=theNode.getpEFTable(this, ndof);
            for(int i=0; i<eft.length; i++){
                EFT[j]=eft[i];
                ++j;
            }
        }
        return EFT;
    }
    
    public int[] getuEFTable(Domain theDomain){
        
        int ndof=theDomain.getFundamentalSolution().get_u_DOFs();
        int[] EFT = new int[ndof*this.numNodes];
        int[] eft ;
        int j=0;
        for(Iterator<Integer> it=this.theNodesHierarchy.values().iterator(); it.hasNext();){
            int hier = it.next();
            Node theNode;
            theNode=this.theNodes.get(hier);
            eft=theNode.getuEFTable();
            for(int i=0; i<eft.length; i++){
                EFT[j]=eft[i];
                ++j;
            }
        }
        return EFT;
    }
    
    public AbstractMatrix getH(Node colNode, Domain theDomain) {
        AbstractMatrix H; AbstractMatrix h;
        int nrows=theDomain.theFundSol.get_u_DOFs();
        int nclms=this.numNodes*nrows;
        H=new AbstractMatrix(nrows,nclms); H.init();
        FundamentalSolution.theFSdata.setMaterial(theDomain.theMaterial);
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            int hier=this.getHierOfNode(elemNode.getID());
            h=theSIntegrator.IntegrateDp(colNode, theDomain, this, elemNode);
            for(int irow=0; irow<nrows; irow++){
                for(int iclm=0; iclm<nrows; iclm++){
                    H.set(irow, (hier-1)*nrows+iclm, h.get(irow, iclm));
                }
            }
        }
        return H;
    }
    
    public double getM(Node colNode, Domain theDomain){
        double M=0.;
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            if(elemNode.getID()==colNode.getID())M=theSIntegrator.IntegrateM(colNode, theDomain, this, elemNode);    
        }
        return M;
    }

    public AbstractMatrix getGu(Node colNode, Domain theDomain) {
        AbstractMatrix H; AbstractMatrix h;
        int nrows=theDomain.theFundSol.get_p_DOFs();
        int nclms=this.numNodes*nrows;
        H=new AbstractMatrix(nrows,nclms); H.init();
        FundamentalSolution.theFSdata.setMaterial(theDomain.theMaterial);
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            int hier=this.getHierOfNode(elemNode.getID());      
            h=theSIntegrator.IntegrateDu(colNode, theDomain, this, elemNode);
            for(int irow=0; irow<nrows; irow++){
                for(int iclm=0; iclm<nrows; iclm++){
                    H.set(irow, (hier-1)*nrows+iclm, h.get(irow, iclm));
                }
            }
        }
        return H;
    }


    public AbstractMatrix getGv(Node colNode, Domain theDomain) {
        AbstractMatrix H; AbstractMatrix h;
        int nrows=theDomain.theFundSol.get_p_DOFs();
        int nclms=this.numNodes*nrows;
        H=new AbstractMatrix(nrows,nclms); H.init();
        FundamentalSolution.theFSdata.setMaterial(theDomain.theMaterial);
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            int hier=this.getHierOfNode(elemNode.getID());      
            h=theSIntegrator.IntegrateDv(colNode, theDomain, this, elemNode);
            for(int irow=0; irow<nrows; irow++){
                for(int iclm=0; iclm<nrows; iclm++){
                    H.set(irow, (hier-1)*nrows+iclm, h.get(irow, iclm));
                }
            }
        }
        return H;
    }
    
    public AbstractMatrix getHres(Node colNode, Domain theDomain) {
        AbstractMatrix H; AbstractMatrix h;
        int nrows=theDomain.theFundSol.get_u_DOFs();
        int nclms=this.numNodes*nrows;
        H=new AbstractMatrix(nrows,nclms); H.init();
        FundamentalSolution.theFSdata.setMaterial(theDomain.theMaterial);
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            int hier=this.getHierOfNode(elemNode.getID());
            h=theSIntegrator.IntegrateDp_res(colNode, theDomain, this, elemNode);
            for(int irow=0; irow<nrows; irow++){
                for(int iclm=0; iclm<nrows; iclm++){
                    H.set(irow, (hier-1)*nrows+iclm, h.get(irow, iclm));
                }
            }
        }
        return H;
    }
    
    public AbstractMatrix getHdif(Node colNode, Domain theDomain) {
        AbstractMatrix H; AbstractMatrix h;
        int nrows=theDomain.theFundSol.get_u_DOFs();
        int nclms=this.numNodes*nrows;
        H=new AbstractMatrix(nrows,nclms); H.init();
        FundamentalSolution.theFSdata.setMaterial(theDomain.theMaterial);
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            int hier=this.getHierOfNode(elemNode.getID());
            h=theSIntegrator.IntegrateDp_dif(colNode, theDomain, this, elemNode);
            for(int irow=0; irow<nrows; irow++){
                for(int iclm=0; iclm<nrows; iclm++){
                    H.set(irow, (hier-1)*nrows+iclm, h.get(irow, iclm));
                }
            }
        }
        return H;
    }
    
    public AbstractMatrix getNornalVec(Domain theDomain) {
        AbstractMatrix H;
        int nrows=theDomain.theFundSol.get_p_DOFs();
        int nclms=this.numNodes*nrows;
        H=new AbstractMatrix(nclms,1); H.init();
        for(Iterator<Node> it=this.getNodes().values().iterator(); it.hasNext();){
            Node elemNode=it.next();
            for(int irow=0; irow<nclms/nrows; irow++){
                H.set(irow*nrows, 0, this.getNormal(elemNode.getID())[0]);
                H.set(irow*nrows+1, 0, this.getNormal(elemNode.getID())[1]);
                if(nrows==3)H.set(irow*nrows+2, 0, this.getNormal(elemNode.getID())[2]);
            }
        }
        return H;
    }

    public boolean isMain(){
        return this.MainElement;
    }

    public void setMain(boolean b){
        this.MainElement=b;
    }

    public boolean isNeumann(){
        boolean ret=true;
        for(int i=0;i<NeumannElement.length;i++){
            if(!NeumannElement[i])ret=false;
        }
        return ret;
    }
    

    public void setNeumann(){
        for(int i=0;i<NeumannElement.length;i++){
            NeumannElement[i]=true;
        }
    }
    
    
    public boolean isContact(){
        boolean ret=true;
        for(int i=0;i<ContactElement.length;i++){
            if(!ContactElement[i])ret=false;
        }
        return ret;
    }
    

    public void setContact(){
        for(int i=0;i<ContactElement.length;i++){
            ContactElement[i]=true;
        }
    }
    
    
    public boolean isDirichlet(){
        boolean ret=true;
        for(int i=0;i<DirichletElement.length;i++){
            if(!DirichletElement[i])ret=false;
        }
        return ret;
    }
    

    public void setDirichlet(){
        for(int i=0;i<DirichletElement.length;i++){
            DirichletElement[i]=true;
        }
    }

    public int getTypeofElem(){return this.ElemType;}


    abstract public Node getNodePrint(int seq);

    abstract public double getWork(Domain theDomain, int wstep);
    
    abstract public double getWork(Domain theDomain, int wstep, int state);
    
    abstract public double getWork(Domain theDomain, int DispStep,int TracStep,int DispState,int TracState);
    
    abstract public double getTIPower(Domain theDomain, int wstep);
    
    abstract public double getTIPowerDifference(Domain theDomain, int wstep);
    
    abstract public double getPower(Domain theDomain, int wstep);
    
    abstract public double getTIPower(Domain theDomain, int wstep, double dt);
    
    abstract public double getWork_hmg(Domain theDomain);

    abstract public double getWorkX(Domain theDomain, int wstep);

    abstract public double getWorkY(Domain theDomain, int wstep);

    abstract public double getGeneralisedForce(Domain theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getGeneralisedForce(Domain theDomain, int DispStep,int TracStep,int DispState,int TracState, int nodeID, int dof);
    
    abstract public double getGeneralisedForce(Domain theDomain, int wstep, int nodeID, int dof, int state);
    
    abstract public double getGeneralisedForce_hmg(Domain theDomain, int nodeID, int dof);

    abstract public double getSideInequality(Domain theDomain, int step, int bound);

    abstract protected void setElemType();
    
    abstract public boolean PointOnElement(geom.Point aNode);
    
    /** Returns the coordinates of point on element, for which point we have 
     * the minimum distance from point aNode. (description to be check)
      @param aNode    Point for which we search for the minimum distance.
   */
    abstract public double[] getCoordMinDistOfNode(geom.Point aNode);

    abstract public double[] getCoordMinDistOfNode(InterfaceNode aNode);

    abstract public double getMinDistOfNode(geom.Point aNode);

    abstract public double[] getGlobalCoordinates(double[] local);

    abstract public double[] getLocalCoordinates(double[] global);
    
    abstract public double[] getDispLocalonNode(int nodeHier,int step,int wstate);
    
    abstract public double[] getTractionLocalonNode(int nodeHier,int step,int wstate);
    
    abstract public double[] getDispLocalonNodeAux(int nodeHier,int step,int wstate);
    
    abstract public double[] getTractionLocalonNodeAux(int nodeHier,int step,int wstate);
    
    abstract public double getShapeFunctionProduct(int dof, int nodeID_i, int nodeID_j);
    
    public BoundaryType getBCType(){
        jbem.Domain.BoundaryType thetype = BoundaryType.UNDEFINED;
        if(this.isContact()) thetype = BoundaryType.CONTACT;
        if(this.isDirichlet()) thetype = BoundaryType.DIRICHLET;
        if(this.isNeumann()) thetype = BoundaryType.NEUMANN;
        return thetype;
    }
    
    abstract public double getInternalPointDisp(Domain theDomain, ResultPoint theInternalPoint, int wdisp,int step,int wstate);
    
    abstract public double getInternalPointStress(Domain theDomain, ResultPoint theInternalPoint, int wstress ,int step,int wstate);
    
    abstract public double getInternalPointStressAux(Domain theDomain, ResultPoint theInternalPoint, int wstress ,int step,int wstate);
    
    abstract public double[] getNormal(int nid);

}
