/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author pchr
 */
abstract public class InterfaceElement {
    protected static int numofIElems=0;
    protected int id;
    protected Map<Integer,InterfaceNode> theINodes = new HashMap<Integer,InterfaceNode>();
    protected Map<Integer,Integer> theNodesHierarchy = new HashMap<Integer,Integer>();
    protected int numINodes;
    protected SpaceIntegrator theSIntegrator;

    // constructor
    public InterfaceElement(){}

    // methods
    public int getID(){
        return this.id;
    }

    public int[] getBaseNodesIDs(){
        int[] vals = null;
        InterfaceNode theINode ;
        int count=0;
        for(Iterator<InterfaceNode> it=this.theINodes.values().iterator(); it.hasNext();){
            theINode = it.next();
            if(theINode.getBaseID()!=0)count+=1;
        }
        vals = new int[count];
        count=0;
        for(Iterator<InterfaceNode> it=this.theINodes.values().iterator(); it.hasNext();){
            theINode = it.next();
            if(theINode.getBaseID()!=0){vals[count]=theINode.getBaseID(); count+=1;}
        }
        return vals;
    }

    public int getNumNodes(){
        return this.numINodes;
    }

    public Map<Integer,InterfaceNode> getNodes(){
        return this.theINodes;
    }

    public InterfaceNode getINodeHier(int seq){
        int id_=theNodesHierarchy.get(seq);
        return this.theINodes.get(id_);
    }
    
    public int getHierOfINode(int theNode_id){
        int seq=0;
        int i=0;
        for(Iterator<Integer> it=this.theNodesHierarchy.values().iterator(); it.hasNext();){
            ++i;
            if(it.next()==theNode_id){seq=i;}
        }
        return seq;
    }

    public void setSpaceIntegrator(SpaceIntegrator SI){
        this.theSIntegrator=SI;
    }
    
    abstract public double[] getCoordMinDistOfINode(InterfaceNode aNode);

    abstract public double getWorkNormal(Interface theDomain, int wstep);
    
    abstract public double getWorkNormal(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkNormalAUX(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkNormalLinear(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkNormalQuad(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkNormal_split(Interface theDomain, int wstep);

    abstract public double getWorkTangent(Interface theDomain, int wstep);
    
    abstract public double getWorkTangentAUX(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkTangent(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkTangentLinear(Interface theDomain, int wstep, double tau);
    
    abstract public double getWorkTangentQuad(Interface theDomain, int wstep, double tau);

    abstract public double getGeneralisedNormalForce(Interface theDomain, int wstep, int nodeID);

    abstract public double getGeneralisedTangentForce(Interface theDomain, int wstep, int nodeID);
    
    abstract public double getGeneralisedNormalForce(Interface theDomain, int wstep, int nodeID, double tau);

    abstract public double getGeneralisedNormalForceQuad(Interface theDomain, int wstep, int nodeID, double tau);

    abstract public double getGeneralisedNormalForceLinear(Interface theDomain, int wstep, int nodeID, double tau);

    abstract public double getGeneralisedTangentForce(Interface theDomain, int wstep, int nodeID, double tau);

    abstract public double getGeneralisedTangentForceQuad(Interface theDomain, int wstep, int nodeID, double tau);

    abstract public double getGeneralisedTangentForceLinear(Interface theDomain, int wstep, int nodeID, double tau);

    abstract public double getDissipatedDamageIEnergy(Interface theDomain, int wstep);

    abstract public double getDissipatedSlipIEnergy(Interface theDomain, int wstep);
    
    abstract public double getDissipatedSlipIEnergy_min(Interface theDomain, int wstep);
    
    abstract public double getDissipatedPlastIEnergy(Interface theDomain, int wstep);
    
    abstract public double getDissipatedPlastIEnergy_min(Interface theDomain, int wstep);
    
    abstract public double getDissipatedFrictionIEnergy(Interface theDomain, int wstep);
    
    abstract public double getIntegralDamageDrivingForce(Interface theDomain, int wstep);
    
    abstract public double getIntegralSlipDrivingForce(Interface theDomain, int wstep);
    
    abstract public double getIntegralPlastDrivingForce(Interface theDomain, int wstep);

    abstract public double getDissipatedDamageIForce(Interface theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getDissipatedDamageIConstant(Interface theDomain, int wstep, int nodeID, int dof);

    abstract public double getGeneralisedNormalDrivingForce(Interface theDomain, int wstep, int nodeID,int dof);
    
    abstract public double getGeneralisedNormalDrivingForce(Interface theDomain, int wstep, int nodeID,int dof, double tau);
    
    abstract public double getGeneralisedNormalDrivingForce_split(Interface theDomain, int wstep, int nodeID,int dof);

    abstract public double getGeneralisedTangentDrivingForce(Interface theDomain, int wstep, int nodeID,int dof);
    
    abstract public double getGeneralisedTangentDrivingForce(Interface theDomain, int wstep, int nodeID,int dof, double tau);

    abstract public double getDissipatedSlipIForce(Interface theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getDissipatedFrictionIForce(Interface theDomain, int wstep, int nodeID);
    
    abstract public double getDissipatedPlastIForce(Interface theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getDissipatedSlipIConstant(Interface theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getDissipatedPlastIConstant(Interface theDomain, int wstep, int nodeID, int dof);

    abstract public double getGeneralisedSlipDrivingForce(Interface theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getGeneralisedPlastDrivingForce(Interface theDomain, int wstep, int nodeID, int dof);
    
    abstract public double getNewCrack(int wstep);
    
    abstract public double getLength();
}
