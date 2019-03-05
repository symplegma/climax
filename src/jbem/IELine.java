/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import javax.swing.JOptionPane;

/**
 *
 * @author pchr
 */
abstract public class IELine extends InterfaceElement{

    // constructor
    public IELine(){}

    // abstract methods
    abstract public double getShapeFunction(int wSF, double xsi);
    abstract public double getShapeFunction_xsi(int wSF, double xsi);

    // methods
    public double[] getAproximatedNormal(double xsi){
        double[] normal = new double[3];
        double[] nodNormal = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;nodNormal[i]=0.;}
        InterfaceNode theNode;
        for(int i=1; i<=this.numINodes; i++){
            theNode=this.getINodeHier(i);
            nodNormal=theNode.getNormal();

            normal[0]+=getShapeFunction(i,xsi)*nodNormal[0];
            normal[1]+=getShapeFunction(i,xsi)*nodNormal[1];
            normal[2]+=getShapeFunction(i,xsi)*nodNormal[2];


        }
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

    public double getJacobian(double xsi){
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        double[] coords;
        InterfaceNode theNode;
        for(int i=1; i<=this.numINodes; i++){
            theNode=this.getINodeHier(i);
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

    public double getNormalDisplacement(int nodeHier, int step){
        double val = this.getINodeHier(nodeHier).getun()[step];
        if(this.getINodeHier(nodeHier).getTwinID()!=0){
            if(this.getINodeHier(nodeHier).isMain()){
                val = val - this.getINodeHier(nodeHier).getTwin().getun()[step];
            }else{
                val = this.getINodeHier(nodeHier).getTwin().getun()[step] - val;
            }
        }
        return val;
    }
    
    public double getNormalDisplacement(int nodeHier, int step, ViscousMaterial Material_main, ViscousMaterial Material_secondary, double tau){
        double val=this.getINodeHier(nodeHier).getun_(Material_main, tau)[step];
        if(this.getINodeHier(nodeHier).getTwinID()!=0){
            if(this.getINodeHier(nodeHier).isMain()){
                val = this.getINodeHier(nodeHier).getun_(Material_main, tau)[step]-this.getINodeHier(nodeHier).getTwin().getun_(Material_secondary, tau)[step];
            }else{
                val = this.getINodeHier(nodeHier).getTwin().getun_(Material_main, tau)[step]-this.getINodeHier(nodeHier).getun_(Material_secondary, tau)[step];
            }
        }
        return val;
    }

    public double getTangentDisplacement(int nodeHier, int step){
        double val = this.getINodeHier(nodeHier).getut()[step];
        if(this.getINodeHier(nodeHier).getTwinID()!=0){
            if(this.getINodeHier(nodeHier).isMain()){
                val = val - this.getINodeHier(nodeHier).getTwin().getut()[step];
            }else{
                val = this.getINodeHier(nodeHier).getTwin().getut()[step] - val;
            }
        }
        return val;
    }
    
    public double getTangentDisplacement(int nodeHier, int step, ViscousMaterial Material_main, ViscousMaterial Material_secondary, double tau){
        double val=this.getINodeHier(nodeHier).getut_(Material_main, tau)[step];
        if(this.getINodeHier(nodeHier).getTwinID()!=0){
            if(this.getINodeHier(nodeHier).isMain()){
                val = this.getINodeHier(nodeHier).getut_(Material_main, tau)[step]-this.getINodeHier(nodeHier).getTwin().getut_(Material_secondary, tau)[step];
            }else{
                val = this.getINodeHier(nodeHier).getTwin().getut_(Material_main, tau)[step]-this.getINodeHier(nodeHier).getut_(Material_secondary, tau)[step];
            }
        }
        return val;
    }

    public double getDamageVariable(int nodeHier, int step){
        double val = 1. ;
        if(this.getINodeHier(nodeHier).getzEFTable(id)!=0)val = this.getINodeHier(nodeHier).getz(id)[step];
        return val;
    }

    public double getSlipVariable(int nodeHier, int step){
        double val = 0. ;
        if(this.getINodeHier(nodeHier).gets_posEFTable(id)!=0)val = this.getINodeHier(nodeHier).gets(id)[step];
        return val;
    }
    
    public double getPlastVariable(int nodeHier, int step){
        double val = 0. ;
        if(this.getINodeHier(nodeHier).getp_posEFTable(id)!=0)val = this.getINodeHier(nodeHier).getp(id)[step];
        return val;
    }

    public double getWorkNormal(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNormalStoredIEnergy(this,theDomain, wstep);
    }
    
    public double getWorkNormal(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNormalStoredIEnergy(this,theDomain, wstep, tau);
    }
    
    @Override
    public double getIntegralDamageDrivingForce(Interface theDomain, int wstep) {
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDamageDrivingForce(this,theDomain, wstep);
    }
    
    @Override
    public double getIntegralSlipDrivingForce(Interface theDomain, int wstep) {
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateSlipDrivingForce(this,theDomain, wstep);
    }
    
    @Override
    public double getIntegralPlastDrivingForce(Interface theDomain, int wstep) {
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegratePlastDrivingForce(this,theDomain, wstep);
    }
    
    public double getWorkNormalAUX(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNormalStoredIEnergyAUX(this,theDomain, wstep, tau);
    }
    
    public double getWorkNormalQuad(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNormalStoredIEnergyQuad(this,theDomain, wstep, tau);
    }
    
    public double getWorkNormalLinear(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNormalStoredIEnergyLinear(this,theDomain, wstep, tau);
    }
    
    public double getWorkNormal_split(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNormalStoredIEnergy_split(this,theDomain, wstep);
    }
    

    public double getWorkTangent(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTangentialStoredIEnergy(this,theDomain, wstep);
    }
    
    public double getWorkTangent(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTangentialStoredIEnergy(this,theDomain, wstep, tau);
    }
    
    public double getWorkTangentAUX(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTangentialStoredIEnergyAUX(this,theDomain, wstep, tau);
    }
    
    public double getWorkTangentQuad(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTangentialStoredIEnergyQuad(this,theDomain, wstep, tau);
    }
    
    public double getWorkTangentLinear(Interface theDomain, int wstep, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTangentialStoredIEnergyLinear(this,theDomain, wstep, tau);
    }

    public double getGeneralisedNormalForce(Interface theDomain, int wstep, int nodeID){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedNormalForce(this, theDomain, wstep, nodeID);
    }
    
    public double getGeneralisedNormalForce(Interface theDomain, int wstep, int nodeID, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedNormalForce(this, theDomain, wstep, nodeID, tau);
    }
    
    public double getGeneralisedNormalForceQuad(Interface theDomain, int wstep, int nodeID, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedNormalForceQuad(this, theDomain, wstep, nodeID, tau);
    }
    
    public double getGeneralisedNormalForceLinear(Interface theDomain, int wstep, int nodeID, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedNormalForceLinear(this, theDomain, wstep, nodeID, tau);
    }

    public double getGeneralisedTangentForce(Interface theDomain, int wstep, int nodeID){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedTangentForce(this, theDomain, wstep, nodeID);
    }
    
    public double getGeneralisedTangentForce(Interface theDomain, int wstep, int nodeID, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedTangentForce(this, theDomain, wstep, nodeID, tau);
    }
    
    public double getGeneralisedTangentForceQuad(Interface theDomain, int wstep, int nodeID, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedTangentForceQuad(this, theDomain, wstep, nodeID, tau);
    }
    
    public double getGeneralisedTangentForceLinear(Interface theDomain, int wstep, int nodeID, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI. IntegrateIntGeneralisedTangentForceLinear(this, theDomain, wstep, nodeID, tau);
    }
    
    public double getDissipatedDamageIEnergy(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedDamageIEnergy(this,theDomain, wstep);
    }

    public double getDissipatedSlipIEnergy(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedSlipIEnergy(this,theDomain, wstep);
    }
    
    public double getDissipatedSlipIEnergy_min(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedSlipIEnergy_min(this,theDomain, wstep);
    }
    
    public double getDissipatedPlastIEnergy_min(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedPlastIEnergy_min(this,theDomain, wstep);
    }
    
    public double getDissipatedFrictionIEnergy(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedFrictionIEnergy(this,theDomain, wstep);
    }
    
    public double getDissipatedPlastIEnergy(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedPlastIEnergy(this,theDomain, wstep);
    }
    
    public double getStoredIEnergy(Interface theDomain, int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateTangentialStoredIEnergy(this,theDomain, wstep)+quadSI.IntegrateNormalStoredIEnergy(this,theDomain, wstep);
    }


    public double getDissipatedDamageIForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedDamageIForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getDissipatedDamageIConstant(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedDamageIConstant(this, theDomain, wstep, nodeID, dof);
    }

    public double getGeneralisedNormalDrivingForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedNormalDrivingForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getGeneralisedNormalDrivingForce(Interface theDomain, int wstep, int nodeID, int dof, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedNormalDrivingForce(this, theDomain, wstep, nodeID, dof, tau);
    }
    
    public double getGeneralisedNormalDrivingForce_split(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedNormalDrivingForce_split(this, theDomain, wstep, nodeID, dof);
    }

    public double getGeneralisedTangentDrivingForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedTangentDrivingForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getGeneralisedTangentDrivingForce(Interface theDomain, int wstep, int nodeID, int dof, double tau){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedTangentDrivingForce(this, theDomain, wstep, nodeID, dof, tau);
    }

    public double getDissipatedSlipIForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedSlipIForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getDissipatedFrictionIForce(Interface theDomain, int wstep, int nodeID){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedFrictionIForce(this, theDomain, wstep, nodeID);
    }
    
    public double getDissipatedPlastIForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedPlastIForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getDissipatedSlipIConstant(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedSlipIConstant(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getDissipatedPlastIConstant(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateDissipatedPlastIConstant(this, theDomain, wstep, nodeID, dof);
    }

    public double getGeneralisedSlipDrivingForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedSlipDrivingForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getGeneralisedPlastDrivingForce(Interface theDomain, int wstep, int nodeID, int dof){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateIntGeneralisedPlastDrivingForce(this, theDomain, wstep, nodeID, dof);
    }
    
    public double getLength(){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateLength(this);
    }
    
    public double getNewCrack(int wstep){
        SpaceLineIntegrator quadSI = (SpaceLineIntegrator) this.theSIntegrator;
        return quadSI.IntegrateNewCrack(this, wstep);
    }
}
