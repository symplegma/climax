/*******************************************************************************
* Climax.                                                                      *
* Copyright (C) 2009-2017 C.G. Panagiotopoulos [http://www.symplegma.org]      *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): C.G. Panagiotopoulos (pchr76@gmail.com)
// *****************************************************************************
package jfem;
import static java.lang.Math.pow;
import java.util.Iterator;
import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class EBeam2d extends Element {
    private static int numberOfEBeam2d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    // some variables needed for elastoplastic analysis
    MaterialElasticPlasticPoint theMaterialPoint;
    public boolean RotaryInertia=false;
    
    // constructors
    public EBeam2d(int id, Node Node1, Node Node2, Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfEBeam2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        int[] dofs = new int[3];
        dofs[0]=1; dofs[1]=2; dofs[2]=6;
        Node1.setNdofs_ofNode(dofs);
        Node2.setNdofs_ofNode(dofs);
        ndofs = 6;
        dof_per_node=3;
        dimension=1;
        
        double[] coords = new double[2];
        coords[0]=Node1.getCoords()[0]+Node2.getCoords()[0];
        coords[1]=Node1.getCoords()[1]+Node2.getCoords()[1];
        theMaterialPoint = new MaterialElasticPlasticPoint();
        theMaterialPoint.setMaterial(Mater);
        theMaterialPoint.initMatPoint(1);
        theMaterialPoint.setCoords(coords,1.0);
    }
    
    // methods
    public int getNumberOfEBeam2d() {
        return numberOfEBeam2d;
    }
    
    public int getElementNdofs() {
        return ndofs;
    }  
    
    public AbstractMatrix getK() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double Ep=E;
//        if(this.ElementMaterial.getType()==0){
//            Ep=((ElasticMaterial) ElementMaterial).getElasticity();
//        } else {
//            ElastoPlastic aMaterial =(ElastoPlastic) this.ElementMaterial;
//            if(this.ElementPlastified()){
//                Ep=( aMaterial.getElasticity()*(aMaterial.getKisotropic() +aMaterial.getKkinematic()) )/
//                  ( aMaterial.getElasticity()+aMaterial.getKisotropic()+aMaterial.getKkinematic() );
//            } else {
//                Ep=aMaterial.getElasticity();
//            }
//        }
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        
        double F2;
        if(As2==0.){
            F2=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F2=(12.*E*J3)/(G*As2*L*L); // Thas it for Timoshenko beam assumption
        }
        double Flex=E*J3/(L*L*L*(1+F2));
        
        MyMatrix.putVal(0, 0,Ep*A/L);
        MyMatrix.putVal(1, 0,0.);
        MyMatrix.putVal(2, 0,0.);
        MyMatrix.putVal(3, 0,-Ep*A/L);
        MyMatrix.putVal(4, 0,0.);
        MyMatrix.putVal(5, 0,0.);
        
        MyMatrix.putVal(0, 1,0.);
        MyMatrix.putVal(1, 1,12.*Flex);
        MyMatrix.putVal(2, 1,6.*L*Flex);
        MyMatrix.putVal(3, 1,0.);
        MyMatrix.putVal(4, 1,-12.*Flex);
        MyMatrix.putVal(5, 1,6.*L*Flex);
        
        MyMatrix.putVal(0, 2,0.);
        MyMatrix.putVal(1, 2,6.*L*Flex);
        MyMatrix.putVal(2, 2,(4.+F2)*L*L*Flex);
        MyMatrix.putVal(3, 2,0.);
        MyMatrix.putVal(4, 2,-6.*L*Flex);
        MyMatrix.putVal(5, 2,(2.-F2)*L*L*Flex);
        
        MyMatrix.putVal(0, 3,-Ep*A/L);
        MyMatrix.putVal(1, 3,0.);
        MyMatrix.putVal(2, 3,0.);
        MyMatrix.putVal(3, 3,Ep*A/L);
        MyMatrix.putVal(4, 3,0.);
        MyMatrix.putVal(5, 3,0.);
        
        MyMatrix.putVal(0, 4,0.);
        MyMatrix.putVal(1, 4,-12.*Flex);
        MyMatrix.putVal(2, 4,-6.*L*Flex);
        MyMatrix.putVal(3, 4,0.);
        MyMatrix.putVal(4, 4,12.*Flex);
        MyMatrix.putVal(5, 4,-6.*L*Flex);
                
        MyMatrix.putVal(0, 5,0.);
        MyMatrix.putVal(1, 5,6.*L*Flex);
        MyMatrix.putVal(2, 5,(2.-F2)*L*L*Flex);
        MyMatrix.putVal(3, 5,0.);
        MyMatrix.putVal(4, 5,-6.*L*Flex);
        MyMatrix.putVal(5, 5,(4.+F2)*L*L*Flex);

        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getM(){
        //https://www.colorado.edu/engineering/CAS/courses.d/MFEMD.d/MFEMD.Ch18.d/MFEMD.Ch18.pdf
        //https://link.springer.com/article/10.1007/s11831-014-9108-x
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double L=this.getL();
        
        double p=((ElasticMaterial) this.ElementMaterial).getDensity();
        double E=((ElasticMaterial) this.ElementMaterial).getElasticity();
        double G=((ElasticMaterial) this.ElementMaterial).getG();
        double A=this.theCrossSection.getA();
        double J3=this.theCrossSection.getJ3();
        double As2=this.theCrossSection.getAs2();
        double m=p*A*L;
        double F2;
        if(As2==0.){
            F2=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F2=(12.*E*J3)/(G*As2*L*L); // Thas it for Timoshenko beam assumption
        }
        double r=Math.sqrt(J3/A);
        double c=m/((1+F2)*(1+F2)); // it was  c=m*L/((1+F2)*(1+F2));
        double[] mcoefs = new double[10];
        mcoefs[0]=13./35.+7*F2/10.+F2*F2/3.;  mcoefs[1]=11*L/210.+(11*F2/120.+F2*F2/24.)*L;
        mcoefs[2]=9/70.+3*F2/10.+F2*F2/6.; mcoefs[3]=13*L/420.+(3*F2/40.+F2*F2/24.)*L;
        mcoefs[4]=L*L/105.+(F2/60.+F2*F2/120.)*L*L; mcoefs[5]=L*L/140+(F2/60.+F2*F2/120.)*L*L;
        mcoefs[6]=6./5.; mcoefs[7]=(1./10.-F2/2.)*L;
        mcoefs[8]=(2./15.+F2/6.+F2*F2/3.)*L*L; mcoefs[9]=(-1./30-F2/6.+F2*F2/6.)*L*L;
        double r2lq=r*r/(L*L); if(!RotaryInertia)r2lq=m*1.0e-12;
        
        
        if(this.consistentM){
            MyMatrix.putVal(0, 0,m/3);
            MyMatrix.putVal(1, 0,0.);
            MyMatrix.putVal(2, 0,0.);
            MyMatrix.putVal(3, 0,m/6);
            MyMatrix.putVal(4, 0,0.);
            MyMatrix.putVal(5, 0,0.);

            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(1,1,c*(mcoefs[0]+mcoefs[6]*r2lq));
            MyMatrix.putVal(2,1,c*(mcoefs[1]+mcoefs[7]*r2lq));
            MyMatrix.putVal(3,1,0.);
            MyMatrix.putVal(4,1,c*(mcoefs[2]-mcoefs[6]*r2lq));
            MyMatrix.putVal(5,1,c*(-mcoefs[3]+mcoefs[7]*r2lq));

            MyMatrix.putVal(0,2,0.);
            MyMatrix.putVal(1,2,c*(mcoefs[1]+mcoefs[7]*r2lq));
            MyMatrix.putVal(2,2,c*(mcoefs[4]+mcoefs[8]*r2lq));
            MyMatrix.putVal(3,2,0.);
            MyMatrix.putVal(4,2,c*(mcoefs[3]-mcoefs[7]*r2lq));
            MyMatrix.putVal(5,2,c*(-mcoefs[5]+mcoefs[9]*r2lq));

            MyMatrix.putVal(0,3,m/6);
            MyMatrix.putVal(1,3,0.);
            MyMatrix.putVal(2,3,0.);
            MyMatrix.putVal(3,3,m/3);
            MyMatrix.putVal(4,3,0.);
            MyMatrix.putVal(5,3,0.);

            MyMatrix.putVal(0,4,0.);
            MyMatrix.putVal(1,4,c*(mcoefs[2]-mcoefs[6]*r2lq));
            MyMatrix.putVal(2,4,c*(mcoefs[3]-mcoefs[7]*r2lq));
            MyMatrix.putVal(3,4,0.);
            MyMatrix.putVal(4,4,c*(mcoefs[0]+mcoefs[6]*r2lq));
            MyMatrix.putVal(5,4,c*(-mcoefs[1]-mcoefs[7]*r2lq));

            MyMatrix.putVal(0,5,0.);
            MyMatrix.putVal(1,5,c*(-mcoefs[3]+mcoefs[7]*r2lq));
            MyMatrix.putVal(2,5,c*(-mcoefs[5]+mcoefs[9]*r2lq));
            MyMatrix.putVal(3,5,0.);
            MyMatrix.putVal(4,5,c*(-mcoefs[1]-mcoefs[7]*r2lq));
            MyMatrix.putVal(5,5,c*(mcoefs[4]+mcoefs[8]*r2lq));
        }else{
            m=m/2;
            MyMatrix.putVal(0, 0,m);
            MyMatrix.putVal(1, 0,0.);
            MyMatrix.putVal(2, 0,0.);
            MyMatrix.putVal(3, 0,0.);
            MyMatrix.putVal(4, 0,0.);
            MyMatrix.putVal(5, 0,0.);

            MyMatrix.putVal(0,1,0.);
            MyMatrix.putVal(1,1,m);
            MyMatrix.putVal(2,1,0.);
            MyMatrix.putVal(3,1,0.);
            MyMatrix.putVal(4,1,0.);
            MyMatrix.putVal(5,1,0.);

            MyMatrix.putVal(0,2,0.);
            MyMatrix.putVal(1,2,0.);
            MyMatrix.putVal(2,2,r2lq);
            MyMatrix.putVal(3,2,0.);
            MyMatrix.putVal(4,2,0.);
            MyMatrix.putVal(5,2,0.);

            MyMatrix.putVal(0,3,0.);
            MyMatrix.putVal(1,3,0.);
            MyMatrix.putVal(2,3,0.);
            MyMatrix.putVal(3,3,m);
            MyMatrix.putVal(4,3,0.);
            MyMatrix.putVal(5,3,0.);

            MyMatrix.putVal(0,4,0.);
            MyMatrix.putVal(1,4,0.);
            MyMatrix.putVal(2,4,0.);
            MyMatrix.putVal(3,4,0.);
            MyMatrix.putVal(4,4,m);
            MyMatrix.putVal(5,4,0.);

            MyMatrix.putVal(0,5,0.);
            MyMatrix.putVal(1,5,0.);
            MyMatrix.putVal(2,5,0.);
            MyMatrix.putVal(3,5,0.);
            MyMatrix.putVal(4,5,0.);
            MyMatrix.putVal(5,5,r2lq);
        }
        
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    
    public AbstractMatrix getM_v(){
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsTrial()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
            velcs.addVal(i+2, 0, theNode.getVelcsTrial()[5]);
            i+=3;
        }
        // from global to local disps
//        velcs=this.getT().times(velcs);
        return getM().times(velcs);
    }
    
    public AbstractMatrix getM_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i, 0, theNode.getAcclsTrial()[0]);
            accls.addVal(i+1, 0, theNode.getAcclsTrial()[1]);
            accls.addVal(i+2, 0, theNode.getAcclsTrial()[5]);
            i+=3;
        }
        // from global to local disps
//        accls=this.getT().times(accls);
        return getM().times(accls);
    }
    
    public AbstractMatrix getM_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]);
            i+=3;
        }
        // from global to local disps
//        disps=this.getT().times(disps);
        return getM().times(disps);
    }
    
    public AbstractMatrix getK_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]);
            i+=3;
        }
        // from global to local disps
//        disps=this.getT().times(disps);
        return getK().times(disps);
    }
    
    public AbstractMatrix getK_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsTrial()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
            velcs.addVal(i+2, 0, theNode.getVelcsTrial()[5]);
            i+=3;
        }
        // from global to local disps
//        velcs=this.getT().times(velcs);
        return getK().times(velcs);
    }
    
    public AbstractMatrix getK_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i, 0, theNode.getAcclsTrial()[0]);
            accls.addVal(i+1, 0, theNode.getAcclsTrial()[1]);
            accls.addVal(i+2, 0, theNode.getAcclsTrial()[5]);
            i+=3;
        }
        // from global to local disps
//        accls=this.getT().times(accls);
        return getK().times(accls);
    }
    
    public int[] getFtable(){
            int[] theEFTable = new int[ndofs];
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theEFTable[j]=theNode.getFtable()[0];
            theEFTable[j+1]=theNode.getFtable()[1];
            theEFTable[j+2]=theNode.getFtable()[5];
            j=j+3;
        }
        return theEFTable;
    }
    
    public double getL(){
            double x1=0.;
            double x2=0.;
            double y1=0.;
            double y2=0.;
            int j=0;
            for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
                Node theNode = it.next();
                switch(j){
                    case 0:
                        x1=theNode.getCoords()[0];
                        y1=theNode.getCoords()[1];
                        break;
                    case 1:
                        x2=theNode.getCoords()[0];
                        y2=theNode.getCoords()[1];
                        break;
                }
                ++j;
            }
            return Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        }
    
    public AbstractMatrix getT() {
            TrMatrix = new AbstractMatrix(ndofs,ndofs);
            double L=this.getL();
            double x1=0.;
            double x2=0.;
            double y1=0.;
            double y2=0.;
            int j=0;
            for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
                Node theNode = it.next();
                switch(j){
                    case 0:
                        x1=theNode.getCoords()[0];
                        y1=theNode.getCoords()[1];
                        break;
                    case 1:
                        x2=theNode.getCoords()[0];
                        y2=theNode.getCoords()[1];
                        break;
                }
                ++j;
            }
            double c=(x2-x1)/L;
            double s=(y2-y1)/L;
            
            TrMatrix.putVal(0, 0,c);
            TrMatrix.putVal(1, 0,-s);
            TrMatrix.putVal(2, 0,0.);
            TrMatrix.putVal(3, 0,0.);
            TrMatrix.putVal(4, 0,0.);
            TrMatrix.putVal(5, 0,0.);
            
            TrMatrix.putVal(0, 1,s);
            TrMatrix.putVal(1, 1,c);
            TrMatrix.putVal(2, 1,0.);
            TrMatrix.putVal(3, 1,0.);
            TrMatrix.putVal(4, 1,0.);
            TrMatrix.putVal(5, 1,0.);
            
            TrMatrix.putVal(0, 2,0.);
            TrMatrix.putVal(1, 2,0.);
            TrMatrix.putVal(2, 2,1.);
            TrMatrix.putVal(3, 2,0.);
            TrMatrix.putVal(4, 2,0.);
            TrMatrix.putVal(5, 2,0.);
            
            TrMatrix.putVal(0, 3,0.);
            TrMatrix.putVal(1, 3,0.);
            TrMatrix.putVal(2, 3,0.);
            TrMatrix.putVal(3, 3,c);
            TrMatrix.putVal(4, 3,-s);
            TrMatrix.putVal(5, 3,0.);
            
            TrMatrix.putVal(0, 4,0.);
            TrMatrix.putVal(1, 4,0.);
            TrMatrix.putVal(2, 4,0.);
            TrMatrix.putVal(3, 4,s);
            TrMatrix.putVal(4, 4,c);
            TrMatrix.putVal(5, 4,0.);
            
            TrMatrix.putVal(0, 5,0.);
            TrMatrix.putVal(1, 5,0.);
            TrMatrix.putVal(2, 5,0.);
            TrMatrix.putVal(3, 5,0.);
            TrMatrix.putVal(4, 5,0.);
            TrMatrix.putVal(5, 5,1.);
            
            return TrMatrix;
    }

    @Override
    void clear() {
    }

    @Override
    void commit() {
    }

    @Override
    public boolean ElementPlastified() {
        return false;
    }

    @Override
    public double getuKu() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
//            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
//            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
//            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]);
            disps.addVal(i, 0, theNode.getDisp(1));
            disps.addVal(i+1, 0, theNode.getDisp(2));
            disps.addVal(i+2, 0, theNode.getDisp(6));
            i+=3;
        }
        // from global to local disps
        disps=this.getT().times(disps);
        return disps.transpose().times(getK().times(disps)).get(0, 0);
    }

    @Override
    public double getvMv() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
//            velcs.addVal(i, 0, theNode.getVelcsTrial()[0]);
//            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
//            velcs.addVal(i+2, 0, theNode.getVelcsTrial()[5]);
            velcs.addVal(i, 0, theNode.getVelc(1));
            velcs.addVal(i+1, 0, theNode.getVelc(2));
            velcs.addVal(i+2, 0, theNode.getVelc(6));
            i+=3;
        }
        // from global to local disps
        velcs=this.getT().times(velcs);
        return velcs.transpose().times(getM().times(velcs)).get(0, 0);
    }

//    @Override
//    public Node getNodeHierarchy(int nodeID) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getVolume(){
        return this.getL()*this.theCrossSection.getA();
    }
    
    public double[] getNormal(int nid){
        double[] normal = new double[3];
        double[] vxsi = new double[3];
        double[] veta = new double[3];
        for(int i=0; i<3; i++){normal[i]=0.;vxsi[i]=0.;veta[i]=0.;}
        vxsi[0]=this.getNodeHierarchy(2).getCoords()[0]-this.getNodeHierarchy(1).getCoords()[0];
        vxsi[1]=this.getNodeHierarchy(2).getCoords()[1]-this.getNodeHierarchy(1).getCoords()[1];
        vxsi[2]=0.0;
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
    
//    public AbstractMatrix getF() {
//        double[] vec = new double[1];
//        double[] vec2 = new double[1];
//        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
//        int i=0;
//        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
//            Node theNode = it.next();
//            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
//            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
//            disps.addVal(i+1, 0, theNode.getDispsTrial()[5]);
//            i+=3;
//        }
//        // from global to local disps
//        disps=this.getT().times(disps);
//        
//        //ex
//        vec[0]=(-disps.get(0, 0)+disps.get(2, 0))/this.getL();
//        vec2[0]=this.theMaterialPoint.giveStrainTakeStress(vec)[0];
//        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
//        Fint.putVal(0, 0, -vec2[0]*this.theCrossSection.getA());
//        Fint.putVal(2, 0, vec2[0]*this.theCrossSection.getA());
//        Fint=this.getT().transpose().times(Fint);
//        return Fint;
//        //return this.getK_u();
//    }
    
    public AbstractMatrix getF() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]);
            i+=3;
        }
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        

        Fint=this.getK().times(disps);
        return Fint;
    }
    
    @Override
    public AbstractMatrix getF(int LC, int step){
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getLoadCaseDisps(LC,step)[0]);
            disps.addVal(i+1, 0, theNode.getLoadCaseDisps(LC,step)[1]);
            disps.addVal(i+2, 0, theNode.getLoadCaseDisps(LC,step)[5]);
            i+=3;
        }
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);

        Fint=this.getK().times(disps);
        return Fint;
    }
    
    @Override
    public AbstractMatrix getM_upre() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]-theNode.getDeltaDisps()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]-theNode.getDeltaDisps()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]-theNode.getDeltaDisps()[5]);
            i+=3;
        }
        // from global to local disps
//        disps=this.getT().times(disps);
        return getM().times(disps);
    }
    
    @Override
    public AbstractMatrix getK_upre() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]-theNode.getDeltaDisps()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]-theNode.getDeltaDisps()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]-theNode.getDeltaDisps()[5]);
            i+=3;
        }
        // from global to local disps
//        disps=this.getT().times(disps);
        return getK().times(disps);
    }
    
    @Override
    public double ShapeFunction(int wsf, double xsi){
        double N;
        switch(wsf){
            case 1: N=0.5*(1.0-xsi); break;
            case 2: N=0.5*(1.0+xsi); break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    @Override
    public AbstractMatrix getFequivalent(LoadDist theLoad){
        AbstractMatrix Fequivalent = new AbstractMatrix(ndofs,1);
        double val = theLoad.getLoadValue();
        int wd = theLoad.getwdof();
        double L=this.getL();
        switch(wd){
            case 1:
                Fequivalent.set(0, 0, val*L/2.0);
                Fequivalent.set(3, 0, val*L/2.0);
                break;
            case 2:
                Fequivalent.set(1, 0, val*L/2.0);
                Fequivalent.set(4, 0, val*L/2.0);
                Fequivalent.set(2, 0, val*L*L/12.0);
                Fequivalent.set(5, 0, -val*L*L/12.0);
                break;
            default:
                System.err.println("Not defined equivalent load for element "+this.getClass().getName()+", in respect to dof "+wd);
        }
        return Fequivalent;
    }
    
    @Override
     public double[] getStress(double[] at, int LC,int step){
         double[] stress = new double[1];
         stress[0]=0.0;
         return stress;
     }
    
    @Override
    public double getStressVonMisses(double[] at, int LC,int step){
        double[] Stress = this.getStress(at, LC, step);
        return Stress[0];
    }
    
    public void setRotaryInertiaContrinution(boolean bool){this.RotaryInertia=bool;}
    
    @Override
    public AbstractMatrix getDK(int param, int pid) {
        //param stand for
        // 1: Elasticity
        // 2: Density
        // 3: Momment
        // 4: Area
        
        // 5: Height of Cross-Section
        AbstractMatrix AM=new AbstractMatrix(ndofs,ndofs,0.0);
        switch(param){
            case 1:
                if(this.ElementMaterial.getID()==pid || pid==0)AM = getDK_1();
                break;
            case 2:
                if(this.ElementMaterial.getID()==pid || pid==0)AM = getDK_2();
                break;
            case 3:
                if(this.theCrossSection.getID()==pid || pid==0)AM = getDK_3();
                break;
            case 4:
                if(this.theCrossSection.getID()==pid || pid==0)AM = getDK_4();
                break;
            case 5:
                if(this.theCrossSection.getID()==pid || pid==0)AM = getDK_5();
                break;
        }
        return AM;
    }
    
    public AbstractMatrix getDK_5() {
        return getDK_3().times(theCrossSection.DJ3_height()).plus(getDK_4().times(theCrossSection.DA_height()));
    }
    
    public AbstractMatrix getDK_1() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            MyMatrix.putVal(0,0,A/L);
            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(0,2,0);
            MyMatrix.putVal(0,3,-(A/L));
            MyMatrix.putVal(0,4,0);
            MyMatrix.putVal(0,5,0);
            
            MyMatrix.putVal(1,0,0);
            MyMatrix.putVal(1,1,(12*J3)/Power(L,3));
            MyMatrix.putVal(1,2,(6*J3)/Power(L,2));
            MyMatrix.putVal(1,3,0);
            MyMatrix.putVal(1,4,(-12*J3)/Power(L,3));
            MyMatrix.putVal(1,5,(6*J3)/Power(L,2));
            
            MyMatrix.putVal(2,0,0);
            MyMatrix.putVal(2,1,(6*J3)/Power(L,2));
            MyMatrix.putVal(2,2,(4*J3)/L);
            MyMatrix.putVal(2,3,0);
            MyMatrix.putVal(2,4,(-6*J3)/Power(L,2));
            MyMatrix.putVal(2,5,(2*J3)/L);
            
            MyMatrix.putVal(3,0,-(A/L));
            MyMatrix.putVal(3,1,0);
            MyMatrix.putVal(3,2,0);
            MyMatrix.putVal(3,3,A/L);
            MyMatrix.putVal(3,4,0);
            MyMatrix.putVal(3,5,0);
            
            MyMatrix.putVal(4,0,0);
            MyMatrix.putVal(4,1,(-12*J3)/Power(L,3));
            MyMatrix.putVal(4,2,(-6*J3)/Power(L,2));
            MyMatrix.putVal(4,3,0);
            MyMatrix.putVal(4,4,(12*J3)/Power(L,3));
            MyMatrix.putVal(4,5,(-6*J3)/Power(L,2));
            
            MyMatrix.putVal(5,0,0);
            MyMatrix.putVal(5,1,(6*J3)/Power(L,2));
            MyMatrix.putVal(5,2,(2*J3)/L);
            MyMatrix.putVal(5,3,0);
            MyMatrix.putVal(5,4,(-6*J3)/Power(L,2));
            MyMatrix.putVal(5,5,(4*J3)/L);
        }else {
            MyMatrix.putVal(0,0,A/L);
            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(0,2,0);
            MyMatrix.putVal(0,3,-(A/L));
            MyMatrix.putVal(0,4,0);
            MyMatrix.putVal(0,5,0);
            
            MyMatrix.putVal(1,0,0);
            MyMatrix.putVal(1,1,(12*A*J3*k2)/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            MyMatrix.putVal(1,2,(6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            MyMatrix.putVal(1,3,0);
            MyMatrix.putVal(1,4,(-12*A*J3*k2)/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            MyMatrix.putVal(1,5,(6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
           
            MyMatrix.putVal(2,0,0);
            MyMatrix.putVal(2,1,(6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            MyMatrix.putVal(2,2,(4*J3*(A*k2*Power(L,2) + 6*J3*(1 + nu)))/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            MyMatrix.putVal(2,3,0);
            MyMatrix.putVal(2,4,(-6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            MyMatrix.putVal(2,5,(2*J3*(A*k2*Power(L,2) - 12*J3*(1 + nu)))/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            
            MyMatrix.putVal(3,0,-(A/L));
            MyMatrix.putVal(3,1,0);
            MyMatrix.putVal(3,2,0);
            MyMatrix.putVal(3,3,A/L);
            MyMatrix.putVal(3,4,0);
            MyMatrix.putVal(3,5,0);
            
            MyMatrix.putVal(4,0,0);
            MyMatrix.putVal(4,1,(-12*A*J3*k2)/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            MyMatrix.putVal(4,2,(-6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            MyMatrix.putVal(4,3,0);
            MyMatrix.putVal(4,4,(12*A*J3*k2)/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            MyMatrix.putVal(4,5,(-6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            
            MyMatrix.putVal(5,0,0);
            MyMatrix.putVal(5,1,(6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            MyMatrix.putVal(5,2,(2*J3*(A*k2*Power(L,2) - 12*J3*(1 + nu)))/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
            MyMatrix.putVal(5,3,0);
            MyMatrix.putVal(5,4,(-6*A*J3*k2)/(A*k2*Power(L,2) + 24*J3*(1 + nu)));
            MyMatrix.putVal(5,5,(4*J3*(A*k2*Power(L,2) + 6*J3*(1 + nu)))/(A*k2*Power(L,3) + 24*J3*L*(1 + nu)));
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getDK_2() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            // zero sensitivity
        }else {
            // zero sensitivity
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getDK_3() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double Eel  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =Eel/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            MyMatrix.putVal(0,0,0);
            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(0,2,0);
            MyMatrix.putVal(0,3,0);
            MyMatrix.putVal(0,4,0);
            MyMatrix.putVal(0,5,0);
            
            MyMatrix.putVal(1,0,0);
            MyMatrix.putVal(1,1,(12*Eel)/Power(L,3));
            MyMatrix.putVal(1,2,(6*Eel)/Power(L,2));
            MyMatrix.putVal(1,3,0);
            MyMatrix.putVal(1,4,(-12*Eel)/Power(L,3));
            MyMatrix.putVal(1,5,(6*Eel)/Power(L,2));
            
            MyMatrix.putVal(2,0,0);
            MyMatrix.putVal(2,1,(6*Eel)/Power(L,2));
            MyMatrix.putVal(2,2,(4*Eel)/L);
            MyMatrix.putVal(2,3,0);
            MyMatrix.putVal(2,4,(-6*Eel)/Power(L,2));
            MyMatrix.putVal(2,5,(2*Eel)/L);
            
            MyMatrix.putVal(3,0,0);
            MyMatrix.putVal(3,1,0);
            MyMatrix.putVal(3,2,0);
            MyMatrix.putVal(3,3,0);
            MyMatrix.putVal(3,4,0);
            MyMatrix.putVal(3,5,0);
            
            MyMatrix.putVal(4,0,0);
            MyMatrix.putVal(4,1,(-12*Eel)/Power(L,3));
            MyMatrix.putVal(4,2,(-6*Eel)/Power(L,2));
            MyMatrix.putVal(4,3,0);
            MyMatrix.putVal(4,4,(12*Eel)/Power(L,3));
            MyMatrix.putVal(4,5,(-6*Eel)/Power(L,2));
            
            MyMatrix.putVal(5,0,0);
            MyMatrix.putVal(5,1,(6*Eel)/Power(L,2));
            MyMatrix.putVal(5,2,(2*Eel)/L);
            MyMatrix.putVal(5,3,0);
            MyMatrix.putVal(5,4,(-6*Eel)/Power(L,2));
            MyMatrix.putVal(5,5,(4*Eel)/L);
        }else {
            MyMatrix.putVal(0,0,0);
            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(0,2,0);
            MyMatrix.putVal(0,3,0);
            MyMatrix.putVal(0,4,0);
            MyMatrix.putVal(0,5,0);
            
            MyMatrix.putVal(1,0,0);
            MyMatrix.putVal(1,1,(12*Power(A,2)*Eel*Power(k2,2)*L)/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(1,2,(6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(1,3,0);
            MyMatrix.putVal(1,4,(-12*Power(A,2)*Eel*Power(k2,2)*L)/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(1,5,(6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            
            MyMatrix.putVal(2,0,0);
            MyMatrix.putVal(2,1,(6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(2,2,(4*Eel*(Power(A,2)*Power(k2,2)*Power(L,4) + 12*A*J3*k2*Power(L,2)*(1 + nu) + 144*Power(J3,2)*Power(1 + nu,2)))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            MyMatrix.putVal(2,3,0);
            MyMatrix.putVal(2,4,(-6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(2,5,(-2*Eel*(-(Power(A,2)*Power(k2,2)*Power(L,4)) + 24*A*J3*k2*Power(L,2)*(1 + nu) + 288*Power(J3,2)*Power(1 + nu,2)))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            
            MyMatrix.putVal(3,0,0);
            MyMatrix.putVal(3,1,0);
            MyMatrix.putVal(3,2,0);
            MyMatrix.putVal(3,3,0);
            MyMatrix.putVal(3,4,0);
            MyMatrix.putVal(3,5,0);
            
            MyMatrix.putVal(4,0,0);
            MyMatrix.putVal(4,1,(-12*Power(A,2)*Eel*Power(k2,2)*L)/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(4,2,(-6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(4,3,0);
            MyMatrix.putVal(4,4,(12*Power(A,2)*Eel*Power(k2,2)*L)/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(4,5,(-6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            
            MyMatrix.putVal(5,0,0);
            MyMatrix.putVal(5,1,(6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(5,2,(-2*Eel*(-(Power(A,2)*Power(k2,2)*Power(L,4)) + 24*A*J3*k2*Power(L,2)*(1 + nu) + 288*Power(J3,2)*Power(1 + nu,2)))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            MyMatrix.putVal(5,3,0);
            MyMatrix.putVal(5,4,(-6*Power(A,2)*Eel*Power(k2,2)*Power(L,2))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(5,5,(4*Eel*(Power(A,2)*Power(k2,2)*Power(L,4) + 12*A*J3*k2*Power(L,2)*(1 + nu) + 144*Power(J3,2)*Power(1 + nu,2)))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getDK_4() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double Eel  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =Eel/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            MyMatrix.putVal(0,0,Eel/L);
            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(0,2,0);
            MyMatrix.putVal(0,3,-(Eel/L));
            MyMatrix.putVal(0,4,0);
            MyMatrix.putVal(0,5,0);
            
            MyMatrix.putVal(1,0,0);
            MyMatrix.putVal(1,1,0);
            MyMatrix.putVal(1,2,0);
            MyMatrix.putVal(1,3,0);
            MyMatrix.putVal(1,4,0);
            MyMatrix.putVal(1,5,0);
            
            MyMatrix.putVal(2,0,0);
            MyMatrix.putVal(2,1,0);
            MyMatrix.putVal(2,2,0);
            MyMatrix.putVal(2,3,0);
            MyMatrix.putVal(2,4,0);
            MyMatrix.putVal(2,5,0);
            
            MyMatrix.putVal(3,0,-(Eel/L));
            MyMatrix.putVal(3,1,0);
            MyMatrix.putVal(3,2,0);
            MyMatrix.putVal(3,3,Eel/L);
            MyMatrix.putVal(3,4,0);
            MyMatrix.putVal(3,5,0);
            
            MyMatrix.putVal(4,0,0);
            MyMatrix.putVal(4,1,0);
            MyMatrix.putVal(4,2,0);
            MyMatrix.putVal(4,3,0);
            MyMatrix.putVal(4,4,0);
            MyMatrix.putVal(4,5,0);
            
            MyMatrix.putVal(5,0,0);
            MyMatrix.putVal(5,1,0);
            MyMatrix.putVal(5,2,0);
            MyMatrix.putVal(5,3,0);
            MyMatrix.putVal(5,4,0);
            MyMatrix.putVal(5,5,0);
        }else {
            MyMatrix.putVal(0,0,Eel/L);
            MyMatrix.putVal(0,1,0);
            MyMatrix.putVal(0,2,0);
            MyMatrix.putVal(0,3,-(Eel/L));
            MyMatrix.putVal(0,4,0);
            MyMatrix.putVal(0,5,0);
            
            MyMatrix.putVal(1,0,0);
            MyMatrix.putVal(1,1,(288*Eel*Power(J3,2)*k2*(1 + nu))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            MyMatrix.putVal(1,2,(144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(1,3,0);
            MyMatrix.putVal(1,4,(-288*Eel*Power(J3,2)*k2*(1 + nu))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            MyMatrix.putVal(1,5,(144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            
            MyMatrix.putVal(2,0,0);
            MyMatrix.putVal(2,1,(144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(2,2,(72*Eel*Power(J3,2)*k2*L*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(2,3,0);
            MyMatrix.putVal(2,4,(-144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(2,5,(72*Eel*Power(J3,2)*k2*L*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            
            MyMatrix.putVal(3,0,-(Eel/L));
            MyMatrix.putVal(3,1,0);
            MyMatrix.putVal(3,2,0);
            MyMatrix.putVal(3,3,Eel/L);
            MyMatrix.putVal(3,4,0);
            MyMatrix.putVal(3,5,0);
            
            MyMatrix.putVal(4,0,0);
            MyMatrix.putVal(4,1,(-288*Eel*Power(J3,2)*k2*(1 + nu))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            MyMatrix.putVal(4,2,(-144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(4,3,0);
            MyMatrix.putVal(4,4,(288*Eel*Power(J3,2)*k2*(1 + nu))/(L*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            MyMatrix.putVal(4,5,(-144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            
            MyMatrix.putVal(5,0,0);
            MyMatrix.putVal(5,1,(144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(5,2,(72*Eel*Power(J3,2)*k2*L*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(5,3,0);
            MyMatrix.putVal(5,4,(-144*Eel*Power(J3,2)*k2*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
            MyMatrix.putVal(5,5,(72*Eel*Power(J3,2)*k2*L*(1 + nu))/Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2));
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    @Override
    public AbstractMatrix getDM(int param, int pid) {
        //param stand for
        // 1: Elasticity
        // 2: Density
        // 3: Momment
        // 4: Area
        
        // 5: Height of Cross-Section
        AbstractMatrix AM=new AbstractMatrix(ndofs,ndofs,0.0);
        if(this.consistentM){
            switch(param){
                case 1:
                    if(this.ElementMaterial.getID()==pid || pid==0)AM = getDM_1();
                    break;
                case 2:
                    if(this.ElementMaterial.getID()==pid || pid==0)AM = getDM_2();
                    break;
                case 3:
                    if(this.theCrossSection.getID()==pid || pid==0)AM = getDM_3();
                    break;
                case 4:
                    if(this.theCrossSection.getID()==pid || pid==0)AM = getDM_4();
                    break;
                case 5:
                    if(this.theCrossSection.getID()==pid || pid==0)AM = getDM_5();
                    break;
            }
        }
        return AM;
    }
    
    public AbstractMatrix getDM_5() {
        return getDM_3().times(theCrossSection.DJ3_height()).plus(getDM_4().times(theCrossSection.DA_height()));
    }
    
    public AbstractMatrix getDM_1() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double p= ((ElasticMaterial) this.ElementMaterial).getDensity();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            // zero sensitivity
        }else {
            // zero sensitivity
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getDM_2() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double p= ((ElasticMaterial) this.ElementMaterial).getDensity();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            if(this.RotaryInertia){
                MyMatrix.putVal(0,0,(A*L)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(A*L)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);
                
                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(A*(13 + 42*Math.sqrt(J3/A))*L)/35.);
                MyMatrix.putVal(1,2,(A*(11 + 21*Math.sqrt(J3/A))*Power(L,2))/210.);
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,-(A*L*(-9 + 7*Math.sqrt(J3/A)*L))/70.);
                MyMatrix.putVal(1,5,(A*(-13 + 42*Math.sqrt(J3/A))*Power(L,2))/420.);
                
                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(A*(11 + 21*Math.sqrt(J3/A))*Power(L,2))/210.);
                MyMatrix.putVal(2,2,(A*(1 + 14*Math.sqrt(J3/A))*Power(L,3))/105.);
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,-(A*(-13 + 42*Math.sqrt(J3/A))*Power(L,2))/420.);
                MyMatrix.putVal(2,5,(A*(-1 + Math.sqrt(J3/A))*Power(L,3))/140.);
                
                MyMatrix.putVal(3,0,(A*L)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(A*L)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);
               
                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,-(A*L*(-9 + 7*Math.sqrt(J3/A)*L))/70.);
                MyMatrix.putVal(4,2,-(A*(-13 + 42*Math.sqrt(J3/A))*Power(L,2))/420.);
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(A*(13 + 42*Math.sqrt(J3/A))*L)/35.);
                MyMatrix.putVal(4,5,-(A*(11 + 21*Math.sqrt(J3/A))*Power(L,2))/210.);
               
                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,(A*(-13 + 42*Math.sqrt(J3/A))*Power(L,2))/420.);
                MyMatrix.putVal(5,2,(A*(-1 + Math.sqrt(J3/A))*Power(L,3))/140.);
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-(A*(11 + 21*Math.sqrt(J3/A))*Power(L,2))/210.);
                MyMatrix.putVal(5,5,(A*(1 + 14*Math.sqrt(J3/A))*Power(L,3))/105.);
            }else{
                MyMatrix.putVal(0,0,(A*L)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(A*L)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);

                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(13*A*L)/35.);
                MyMatrix.putVal(1,2,(11*A*Power(L,2))/210.);
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(9*A*L)/70.);
                MyMatrix.putVal(1,5,(-13*A*Power(L,2))/420.);

                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(11*A*Power(L,2))/210.);
                MyMatrix.putVal(2,2,(A*Power(L,3))/105.);
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(13*A*Power(L,2))/420.);
                MyMatrix.putVal(2,5,-(A*Power(L,3))/140.);

                MyMatrix.putVal(3,0,(A*L)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(A*L)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);

                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(9*A*L)/70.);
                MyMatrix.putVal(4,2,(13*A*Power(L,2))/420.);
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(13*A*L)/35.);
                MyMatrix.putVal(4,5,(-11*A*Power(L,2))/210.);

                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,(-13*A*Power(L,2))/420.);
                MyMatrix.putVal(5,2,-(A*Power(L,3))/140.);
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,(-11*A*Power(L,2))/210.);
                MyMatrix.putVal(5,5,(A*Power(L,3))/105.);
            }
        }else {
            if(this.RotaryInertia){
                MyMatrix.putVal(0,0,(A*L)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(A*L)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);
                
                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(A*L*(Power(A,2)*(13 + 42*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) + 588*A*J3*k2*Power(L,2)*(1 + nu) + 6720*Power(J3,2)*Power(1 + nu,2)))/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(1,2,(A*Power(L,2)*(Power(A,2)*(11 + 21*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) - 42*A*J3*(-11 + 60*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(A*L*(Power(A,2)*Power(k2,2)*Power(L,4)*(9 - 7*Math.sqrt(J3/A)*L) + 168*A*J3*k2*Power(L,2)*(3 + 5*Math.sqrt(J3/A)*L)*(1 + nu) + 6720*Power(J3,2)*Power(1 + nu,2)))/(70.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(1,5,(A*Power(L,2)*(Power(A,2)*(-13 + 42*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) - 252*A*J3*(3 + 20*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) - 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                
                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(A*Power(L,2)*(Power(A,2)*(11 + 21*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) - 42*A*J3*(-11 + 60*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(2,2,(A*Power(L,3)*(Power(A,2)*(1 + 14*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) + 42*A*J3*(1 + 10*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 504*Power(J3,2)*(1 + 40*Math.sqrt(J3/A))*Power(1 + nu,2)))/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(A*Power(L,2)*(Power(A,2)*(13 - 42*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) + 252*A*J3*(3 + 20*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(2,5,(A*(-1 + Math.sqrt(J3/A))*Power(L,3)*(Power(A,2)*Power(k2,2)*Power(L,4) + 56*A*J3*k2*Power(L,2)*(1 + nu) + 672*Power(J3,2)*Power(1 + nu,2)))/(140.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                
                MyMatrix.putVal(3,0,(A*L)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(A*L)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);
                
                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(A*L*(Power(A,2)*Power(k2,2)*Power(L,4)*(9 - 7*Math.sqrt(J3/A)*L) + 168*A*J3*k2*Power(L,2)*(3 + 5*Math.sqrt(J3/A)*L)*(1 + nu) + 6720*Power(J3,2)*Power(1 + nu,2)))/(70.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(4,2,(A*Power(L,2)*(Power(A,2)*(13 - 42*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) + 252*A*J3*(3 + 20*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(A*L*(Power(A,2)*(13 + 42*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) + 588*A*J3*k2*Power(L,2)*(1 + nu) + 6720*Power(J3,2)*Power(1 + nu,2)))/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(4,5,-(A*Power(L,2)*(Power(A,2)*(11 + 21*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) - 42*A*J3*(-11 + 60*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                
                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,(A*Power(L,2)*(Power(A,2)*(-13 + 42*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) - 252*A*J3*(3 + 20*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) - 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(5,2,(A*(-1 + Math.sqrt(J3/A))*Power(L,3)*(Power(A,2)*Power(k2,2)*Power(L,4) + 56*A*J3*k2*Power(L,2)*(1 + nu) + 672*Power(J3,2)*Power(1 + nu,2)))/(140.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-(A*Power(L,2)*(Power(A,2)*(11 + 21*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) - 42*A*J3*(-11 + 60*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(5,5,(A*Power(L,3)*(Power(A,2)*(1 + 14*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4) + 42*A*J3*(1 + 10*Math.sqrt(J3/A))*k2*Power(L,2)*(1 + nu) + 504*Power(J3,2)*(1 + 40*Math.sqrt(J3/A))*Power(1 + nu,2)))/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            }else{
                MyMatrix.putVal(0,0,(A*L)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(A*L)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);

                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(A*L*(13*Power(A,2)*Power(k2,2)*Power(L,4) + 588*A*J3*k2*Power(L,2)*(1 + nu) + 6720*Power(J3,2)*Power(1 + nu,2)))/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(1,2,(A*Power(L,2)*(11*Power(A,2)*Power(k2,2)*Power(L,4) + 462*A*J3*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(3*A*L*(3*Power(A,2)*Power(k2,2)*Power(L,4) + 168*A*J3*k2*Power(L,2)*(1 + nu) + 2240*Power(J3,2)*Power(1 + nu,2)))/(70.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(1,5,-(A*Power(L,2)*(13*Power(A,2)*Power(k2,2)*Power(L,4) + 756*A*J3*k2*Power(L,2)*(1 + nu) + 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));

                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(A*Power(L,2)*(11*Power(A,2)*Power(k2,2)*Power(L,4) + 462*A*J3*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(2,2,(A*Power(L,3)*(Power(A,2)*Power(k2,2)*Power(L,4) + 42*A*J3*k2*Power(L,2)*(1 + nu) + 504*Power(J3,2)*Power(1 + nu,2)))/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(A*Power(L,2)*(13*Power(A,2)*Power(k2,2)*Power(L,4) + 756*A*J3*k2*Power(L,2)*(1 + nu) + 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(2,5,-(A*Power(L,3)*(Power(A,2)*Power(k2,2)*Power(L,4) + 56*A*J3*k2*Power(L,2)*(1 + nu) + 672*Power(J3,2)*Power(1 + nu,2)))/(140.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));

                MyMatrix.putVal(3,0,(A*L)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(A*L)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);

                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(3*A*L*(3*Power(A,2)*Power(k2,2)*Power(L,4) + 168*A*J3*k2*Power(L,2)*(1 + nu) + 2240*Power(J3,2)*Power(1 + nu,2)))/(70.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(4,2,(A*Power(L,2)*(13*Power(A,2)*Power(k2,2)*Power(L,4) + 756*A*J3*k2*Power(L,2)*(1 + nu) + 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(A*L*(13*Power(A,2)*Power(k2,2)*Power(L,4) + 588*A*J3*k2*Power(L,2)*(1 + nu) + 6720*Power(J3,2)*Power(1 + nu,2)))/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(4,5,-(A*Power(L,2)*(11*Power(A,2)*Power(k2,2)*Power(L,4) + 462*A*J3*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));

                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,-(A*Power(L,2)*(13*Power(A,2)*Power(k2,2)*Power(L,4) + 756*A*J3*k2*Power(L,2)*(1 + nu) + 10080*Power(J3,2)*Power(1 + nu,2)))/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(5,2,-(A*Power(L,3)*(Power(A,2)*Power(k2,2)*Power(L,4) + 56*A*J3*k2*Power(L,2)*(1 + nu) + 672*Power(J3,2)*Power(1 + nu,2)))/(140.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-(A*Power(L,2)*(11*Power(A,2)*Power(k2,2)*Power(L,4) + 462*A*J3*k2*Power(L,2)*(1 + nu) + 5040*Power(J3,2)*Power(1 + nu,2)))/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
                MyMatrix.putVal(5,5,(A*Power(L,3)*(Power(A,2)*Power(k2,2)*Power(L,4) + 42*A*J3*k2*Power(L,2)*(1 + nu) + 504*Power(J3,2)*Power(1 + nu,2)))/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),2)));
            }
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getDM_3() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double p= ((ElasticMaterial) this.ElementMaterial).getDensity();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            if(this.RotaryInertia){
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,0);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);
                
                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(3*L*p)/(5.*Math.sqrt(J3/A)));
                MyMatrix.putVal(1,2,(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,-(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(1,5,(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                
                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(2,2,(Power(L,3)*p)/(15.*Math.sqrt(J3/A)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,-(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(2,5,(Power(L,3)*p)/(280.*Math.sqrt(J3/A)));
                
                MyMatrix.putVal(3,0,0);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,0);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);
                
                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,-(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(4,2,-(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(3*L*p)/(5.*Math.sqrt(J3/A)));
                MyMatrix.putVal(4,5,-(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                
                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(5,2,(Power(L,3)*p)/(280.*Math.sqrt(J3/A)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-(Power(L,2)*p)/(20.*Math.sqrt(J3/A)));
                MyMatrix.putVal(5,5,(Power(L,3)*p)/(15.*Math.sqrt(J3/A)));
            }else{
                // zero sensitivity
            }
        }else {
            if(this.RotaryInertia){
                MyMatrix.putVal(0,0,0);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,0);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);
                
                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(-3*A*J3*k2*Power(L,3)*(56*J3*(1 + nu)*(9*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + A*k2*Power(L,2)*(-7*k2*Power(L,2) + 12*Math.sqrt(J3/A)*(1 + nu)))*p)/(35.*Power(J3/A,1.5)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,2,-(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(A*k2*Power(L,3)*(-20160*Power(J3,2)*L*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,3) + 8*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,3) + 144*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,5,-(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                
                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,-(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,2,(Power(L,3)*(126*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 45360*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 241920*Power(J3,3)*Power(1 + nu,3) + Power(A,3)*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) - 6*Math.sqrt(J3/A)*(1 + nu)))*p)/(105.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,5,(Power(L,3)*(96*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 2016*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 16128*Power(J3,3)*Power(1 + nu,3) + Power(A,3)*Power(k2,2)*Power(L,4)*(k2*Power(L,2) - 16*Math.sqrt(J3/A)*(1 + nu)))*p)/(280.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                
                MyMatrix.putVal(3,0,0);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,0);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);
                
                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(A*k2*Power(L,3)*(-20160*Power(J3,2)*L*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,3) + 8*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,3) + 144*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,2,(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(-3*A*J3*k2*Power(L,3)*(56*J3*(1 + nu)*(9*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + A*k2*Power(L,2)*(-7*k2*Power(L,2) + 12*Math.sqrt(J3/A)*(1 + nu)))*p)/(35.*Power(J3/A,1.5)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,5,(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                
                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,-(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,2,(Power(L,3)*(96*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 2016*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 16128*Power(J3,3)*Power(1 + nu,3) + Power(A,3)*Power(k2,2)*Power(L,4)*(k2*Power(L,2) - 16*Math.sqrt(J3/A)*(1 + nu)))*p)/(280.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,(A*k2*Power(L,4)*(-20160*Power(J3,2)*Power(1 + nu,2) + 336*A*J3*(1 + nu)*(9*k2*Power(L,2) + 2*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*k2*Power(L,2)*(-7*k2*Power(L,2) + 44*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,5,(Power(L,3)*(126*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 45360*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 241920*Power(J3,3)*Power(1 + nu,3) + Power(A,3)*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) - 6*Math.sqrt(J3/A)*(1 + nu)))*p)/(105.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
            }else{
                MyMatrix.putVal(0,0,0);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,0);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);

                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(-12*Power(A,2)*k2*Power(L,3)*(1 + nu)*(3*A*k2*Power(L,2) + 56*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,2,-(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(12*Power(A,2)*k2*Power(L,3)*(1 + nu)*(3*A*k2*Power(L,2) + 56*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,5,-(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));

                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,-(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,2,(-2*Power(A,3)*Power(k2,2)*Power(L,7)*(1 + nu)*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,5,(-2*Power(A,3)*Power(k2,2)*Power(L,7)*(1 + nu)*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));

                MyMatrix.putVal(3,0,0);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,0);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);

                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(12*Power(A,2)*k2*Power(L,3)*(1 + nu)*(3*A*k2*Power(L,2) + 56*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,2,(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(-12*Power(A,2)*k2*Power(L,3)*(1 + nu)*(3*A*k2*Power(L,2) + 56*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,5,(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));

                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,-(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,2,(-2*Power(A,3)*Power(k2,2)*Power(L,7)*(1 + nu)*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,(Power(A,2)*k2*Power(L,4)*(1 + nu)*(11*A*k2*Power(L,2) + 168*J3*(1 + nu))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,5,(-2*Power(A,3)*Power(k2,2)*Power(L,7)*(1 + nu)*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
            }
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getDM_4() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double p= ((ElasticMaterial) this.ElementMaterial).getDensity();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        double k2=As2/A;
        
        if(As2==0.){
            if(this.RotaryInertia){
                MyMatrix.putVal(0,0,(L*p)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(L*p)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);

                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,((13 + 21*Math.sqrt(J3/A))*L*p)/35.);
                MyMatrix.putVal(1,2,((22 + 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(L*(18 - 7*Math.sqrt(J3/A)*L)*p)/140.);
                MyMatrix.putVal(1,5,((-13 + 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);

                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,((22 + 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);
                MyMatrix.putVal(2,2,((1 + 7*Math.sqrt(J3/A))*Power(L,3)*p)/105.);
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,((13 - 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);
                MyMatrix.putVal(2,5,((-2 + Math.sqrt(J3/A))*Power(L,3)*p)/280.);

                MyMatrix.putVal(3,0,(L*p)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(L*p)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);

                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(L*(18 - 7*Math.sqrt(J3/A)*L)*p)/140.);
                MyMatrix.putVal(4,2,((13 - 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,((13 + 21*Math.sqrt(J3/A))*L*p)/35.);
                MyMatrix.putVal(4,5,-((22 + 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);

                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,((-13 + 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);
                MyMatrix.putVal(5,2,((-2 + Math.sqrt(J3/A))*Power(L,3)*p)/280.);
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-((22 + 21*Math.sqrt(J3/A))*Power(L,2)*p)/420.);
                MyMatrix.putVal(5,5,((1 + 7*Math.sqrt(J3/A))*Power(L,3)*p)/105.);
            }else{
                MyMatrix.putVal(0,0,(L*p)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(L*p)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);

                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(13*L*p)/35.);
                MyMatrix.putVal(1,2,(11*Power(L,2)*p)/210.);
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(9*L*p)/70.);
                MyMatrix.putVal(1,5,(-13*Power(L,2)*p)/420.);

                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(11*Power(L,2)*p)/210.);
                MyMatrix.putVal(2,2,(Power(L,3)*p)/105.);
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(13*Power(L,2)*p)/420.);
                MyMatrix.putVal(2,5,-(Power(L,3)*p)/140.);

                MyMatrix.putVal(3,0,(L*p)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(L*p)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);

                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(9*L*p)/70.);
                MyMatrix.putVal(4,2,(13*Power(L,2)*p)/420.);
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(13*L*p)/35.);
                MyMatrix.putVal(4,5,(-11*Power(L,2)*p)/210.);

                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,(-13*Power(L,2)*p)/420.);
                MyMatrix.putVal(5,2,-(Power(L,3)*p)/140.);
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,(-11*Power(L,2)*p)/210.);
                MyMatrix.putVal(5,5,(Power(L,3)*p)/105.);
            }
        }else {
            if(this.RotaryInertia){
                MyMatrix.putVal(0,0,(L*p)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(L*p)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);
                
                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(L*(13*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 161280*Power(J3,3)*Math.sqrt(J3/A)*Power(1 + nu,3) + 168*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(15*k2*Power(L,2) + 128*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) + 312*Math.sqrt(J3/A)*(1 + nu)))*p)/(35.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,2,(Power(L,2)*(22*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(-3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(5*k2*Power(L,2) + 34*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) + 528*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(L*(18*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 20160*Power(J3,3)*Power(1 + nu,2)*(3*k2*Power(L,3) + 16*Math.sqrt(J3/A)*(1 + nu)) + 336*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(-5*k2*Power(L,3) + 104*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*J3*Power(k2,2)*Power(L,4)*(-7*k2*Power(L,3) + 1296*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,5,-(Power(L,2)*(13*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(-5*k2*Power(L,2) + 26*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(-7*k2*Power(L,2) + 312*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                
                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(Power(L,2)*(22*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(-3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(5*k2*Power(L,2) + 34*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) + 528*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,2,(Power(L,3)*(Power(A,3)*(1 + 7*Math.sqrt(J3/A))*Power(k2,3)*Power(L,6) + 18*Power(A,2)*J3*(4 + 35*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4)*(1 + nu) - 1512*A*Power(J3,2)*(-1 + 10*Math.sqrt(J3/A))*k2*Power(L,2)*Power(1 + nu,2) + 12096*Power(J3,3)*(1 + 20*Math.sqrt(J3/A))*Power(1 + nu,3))*p)/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(Power(L,2)*(13*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(-5*k2*Power(L,2) + 26*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(-7*k2*Power(L,2) + 312*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,5,(Power(L,3)*(Power(A,3)*(-2 + Math.sqrt(J3/A))*Power(k2,3)*Power(L,6) + 16*Power(A,2)*J3*(-9 + 4*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4)*(1 + nu) + 2016*A*Power(J3,2)*(-2 + Math.sqrt(J3/A))*k2*Power(L,2)*Power(1 + nu,2) + 16128*Power(J3,3)*(-2 + Math.sqrt(J3/A))*Power(1 + nu,3))*p)/(280.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                
                MyMatrix.putVal(3,0,(L*p)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(L*p)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);
                
                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(L*(18*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 20160*Power(J3,3)*Power(1 + nu,2)*(3*k2*Power(L,3) + 16*Math.sqrt(J3/A)*(1 + nu)) + 336*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(-5*k2*Power(L,3) + 104*Math.sqrt(J3/A)*(1 + nu)) + Power(A,2)*J3*Power(k2,2)*Power(L,4)*(-7*k2*Power(L,3) + 1296*Math.sqrt(J3/A)*(1 + nu)))*p)/(140.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,2,(Power(L,2)*(13*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(-5*k2*Power(L,2) + 26*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(-7*k2*Power(L,2) + 312*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(L*(13*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 161280*Power(J3,3)*Math.sqrt(J3/A)*Power(1 + nu,3) + 168*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(15*k2*Power(L,2) + 128*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) + 312*Math.sqrt(J3/A)*(1 + nu)))*p)/(35.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,5,-(Power(L,2)*(22*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(-3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(5*k2*Power(L,2) + 34*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) + 528*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                
                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,-(Power(L,2)*(13*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(-5*k2*Power(L,2) + 26*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(-7*k2*Power(L,2) + 312*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,2,(Power(L,3)*(Power(A,3)*(-2 + Math.sqrt(J3/A))*Power(k2,3)*Power(L,6) + 16*Power(A,2)*J3*(-9 + 4*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4)*(1 + nu) + 2016*A*Power(J3,2)*(-2 + Math.sqrt(J3/A))*k2*Power(L,2)*Power(1 + nu,2) + 16128*Power(J3,3)*(-2 + Math.sqrt(J3/A))*Power(1 + nu,3))*p)/(280.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-(Power(L,2)*(22*Power(A,3)*Math.sqrt(J3/A)*Power(k2,3)*Power(L,6) + 60480*Power(J3,3)*Power(1 + nu,2)*(-3*k2*Power(L,2) + 4*Math.sqrt(J3/A)*(1 + nu)) + 1008*A*Power(J3,2)*k2*Power(L,2)*(1 + nu)*(5*k2*Power(L,2) + 34*Math.sqrt(J3/A)*(1 + nu)) + 3*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(7*k2*Power(L,2) + 528*Math.sqrt(J3/A)*(1 + nu)))*p)/(420.*Math.sqrt(J3/A)*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,5,(Power(L,3)*(Power(A,3)*(1 + 7*Math.sqrt(J3/A))*Power(k2,3)*Power(L,6) + 18*Power(A,2)*J3*(4 + 35*Math.sqrt(J3/A))*Power(k2,2)*Power(L,4)*(1 + nu) - 1512*A*Power(J3,2)*(-1 + 10*Math.sqrt(J3/A))*k2*Power(L,2)*Power(1 + nu,2) + 12096*Power(J3,3)*(1 + 20*Math.sqrt(J3/A))*Power(1 + nu,3))*p)/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
            }else{
                MyMatrix.putVal(0,0,(L*p)/3.);
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(0,2,0);
                MyMatrix.putVal(0,3,(L*p)/6.);
                MyMatrix.putVal(0,4,0);
                MyMatrix.putVal(0,5,0);

                MyMatrix.putVal(1,0,0);
                MyMatrix.putVal(1,1,(L*(13*Power(A,3)*Power(k2,3)*Power(L,6) + 936*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 21504*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 161280*Power(J3,3)*Power(1 + nu,3))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,2,(Power(L,2)*(11*Power(A,3)*Power(k2,3)*Power(L,6) + 792*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 17136*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 120960*Power(J3,3)*Power(1 + nu,3))*p)/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,3,0);
                MyMatrix.putVal(1,4,(3*L*(3*Power(A,3)*Power(k2,3)*Power(L,6) + 216*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 5824*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 53760*Power(J3,3)*Power(1 + nu,3))*p)/(70.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(1,5,-(Power(L,2)*(13*Power(A,3)*Power(k2,3)*Power(L,6) + 936*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 26208*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 241920*Power(J3,3)*Power(1 + nu,3))*p)/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));

                MyMatrix.putVal(2,0,0);
                MyMatrix.putVal(2,1,(Power(L,2)*(11*Power(A,3)*Power(k2,3)*Power(L,6) + 792*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 17136*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 120960*Power(J3,3)*Power(1 + nu,3))*p)/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,2,(Power(L,3)*(Power(A,3)*Power(k2,3)*Power(L,6) + 72*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 1512*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 12096*Power(J3,3)*Power(1 + nu,3))*p)/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,3,0);
                MyMatrix.putVal(2,4,(Power(L,2)*(13*Power(A,3)*Power(k2,3)*Power(L,6) + 936*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 26208*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 241920*Power(J3,3)*Power(1 + nu,3))*p)/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(2,5,-(Power(L,3)*(Power(A,3)*Power(k2,3)*Power(L,6) + 72*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 2016*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 16128*Power(J3,3)*Power(1 + nu,3))*p)/(140.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));

                MyMatrix.putVal(3,0,(L*p)/6.);
                MyMatrix.putVal(3,1,0);
                MyMatrix.putVal(3,2,0);
                MyMatrix.putVal(3,3,(L*p)/3.);
                MyMatrix.putVal(3,4,0);
                MyMatrix.putVal(3,5,0);

                MyMatrix.putVal(4,0,0);
                MyMatrix.putVal(4,1,(3*L*(3*Power(A,3)*Power(k2,3)*Power(L,6) + 216*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 5824*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 53760*Power(J3,3)*Power(1 + nu,3))*p)/(70.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,2,(Power(L,2)*(13*Power(A,3)*Power(k2,3)*Power(L,6) + 936*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 26208*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 241920*Power(J3,3)*Power(1 + nu,3))*p)/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,3,0);
                MyMatrix.putVal(4,4,(L*(13*Power(A,3)*Power(k2,3)*Power(L,6) + 936*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 21504*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 161280*Power(J3,3)*Power(1 + nu,3))*p)/(35.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(4,5,-(Power(L,2)*(11*Power(A,3)*Power(k2,3)*Power(L,6) + 792*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 17136*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 120960*Power(J3,3)*Power(1 + nu,3))*p)/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));

                MyMatrix.putVal(5,0,0);
                MyMatrix.putVal(5,1,-(Power(L,2)*(13*Power(A,3)*Power(k2,3)*Power(L,6) + 936*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 26208*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 241920*Power(J3,3)*Power(1 + nu,3))*p)/(420.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,2,-(Power(L,3)*(Power(A,3)*Power(k2,3)*Power(L,6) + 72*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 2016*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 16128*Power(J3,3)*Power(1 + nu,3))*p)/(140.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,3,0);
                MyMatrix.putVal(5,4,-(Power(L,2)*(11*Power(A,3)*Power(k2,3)*Power(L,6) + 792*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 17136*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 120960*Power(J3,3)*Power(1 + nu,3))*p)/(210.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
                MyMatrix.putVal(5,5,(Power(L,3)*(Power(A,3)*Power(k2,3)*Power(L,6) + 72*Power(A,2)*J3*Power(k2,2)*Power(L,4)*(1 + nu) + 1512*A*Power(J3,2)*k2*Power(L,2)*Power(1 + nu,2) + 12096*Power(J3,3)*Power(1 + nu,3))*p)/(105.*Power(A*k2*Power(L,2) + 24*J3*(1 + nu),3)));
            }
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    private double Power(double arg, double pwr){return Math.pow(arg, pwr);}
}
