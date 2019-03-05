/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

/**
 *
 * @author pchr
 */
public class IENode extends InterfaceElement{
    
    // constructor
    public IENode(InterfaceNode Node1){
        ++numofIElems;
        this.id=numofIElems;
        this.theINodes.put(Node1.getID(), Node1);this.theNodesHierarchy.put(1, Node1.getID());
        Node1.putIElement(this);
        this.numINodes=1;
    }

    @Override
    public double[] getCoordMinDistOfINode(InterfaceNode aNode) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkNormal(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWorkNormal(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getIntegralDamageDrivingForce(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getIntegralSlipDrivingForce(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getIntegralPlastDrivingForce(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWorkNormalQuad(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWorkNormalLinear(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkTangent(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWorkTangent(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWorkTangentLinear(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getWorkTangentQuad(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedNormalForce(Interface theDomain, int wstep, int nodeID) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedTangentForce(Interface theDomain, int wstep, int nodeID) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedDamageIEnergy(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedSlipIEnergy(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public double getDissipatedSlipIEnergy_min(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedPlastIEnergy(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedDamageIForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedNormalDrivingForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedTangentDrivingForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedSlipIForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedPlastIForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedSlipDrivingForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedPlastDrivingForce(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedSlipIConstant(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedPlastIConstant(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedDamageIConstant(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedPlastIEnergy_min(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkNormal_split(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedNormalDrivingForce_split(Interface theDomain, int wstep, int nodeID, int dof) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getNewCrack(int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedFrictionIEnergy(Interface theDomain, int wstep) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getDissipatedFrictionIForce(Interface theDomain, int wstep, int nodeID) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getLength() {
        return 0.;
    }

    @Override
    public double getGeneralisedNormalForce(Interface theDomain, int wstep, int nodeID, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedTangentForce(Interface theDomain, int wstep, int nodeID, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedNormalForceQuad(Interface theDomain, int wstep, int nodeID, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedNormalForceLinear(Interface theDomain, int wstep, int nodeID, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedTangentForceQuad(Interface theDomain, int wstep, int nodeID, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedTangentForceLinear(Interface theDomain, int wstep, int nodeID, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedNormalDrivingForce(Interface theDomain, int wstep, int nodeID, int dof, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getGeneralisedTangentDrivingForce(Interface theDomain, int wstep, int nodeID, int dof, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkNormalAUX(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getWorkTangentAUX(Interface theDomain, int wstep, double tau) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
