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

import jmat.AbstractMatrix;
import java.util.Iterator;

/**
 *
 * @author pchr
 */
public class TBeam2d extends Element{
    private static int numberOfTBeam2d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    
    // some variables needed for elastoplastic analysis
    private double[]  Stress=new double[3];
    private double[]  Eplastic=new double[3];
    private double[]  Etotal=new double[3];
    private double[]  Xsi=new double[3];
    private double  YieldFunction=0.;
    private double  aPar=0.;
    private double  qPar_N=0.;
    private double  qPar_M=0.;
    
    private double[]  EplasticTrial=new double[3];
    private double  aParTrial=0.;
    private double  qPar_NTrial=0.;
    private double  qPar_MTrial=0.;
    
    
    boolean plast;
    
    // constructors
    public TBeam2d(int id, Node Node1, Node Node2, Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfTBeam2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        int[] dofs = new int[3];
        dofs[0]=1; dofs[1]=2; dofs[2]=6;
        Node1.setNdofs_ofNode(dofs);
        Node2.setNdofs_ofNode(dofs);
        
        this.EplasticTrial[0]=0.;this.EplasticTrial[1]=0.;this.EplasticTrial[2]=0.;
        dof_per_node=3;
        ndofs = 6;
        dimension=1;
    }

    
    // methods
    private double getL(){
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
    public AbstractMatrix getK() {
        /*
        MyMatrix = new Matrix(ndofs,ndofs);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        
        double[] km = new double[5];
        km[0]= E*A/L;
        km[1]= G*As2/L;
	km[2]= G*As2/2;
	km[3]= E*J3/L+G*As2*L/4;
	km[4]=-E*J3/L+G*As2*L/4;
        
        MyMatrix.putVal(0, 0,km[0]);
        MyMatrix.putVal(1, 0,0.);
        MyMatrix.putVal(2, 0,0.);
        MyMatrix.putVal(3, 0,-km[0]);
        MyMatrix.putVal(4, 0,0.);
        MyMatrix.putVal(5, 0,0.);
        
        MyMatrix.putVal(0, 1,0.);
        MyMatrix.putVal(1, 1,km[1]);
        MyMatrix.putVal(2, 1,km[2]);
        MyMatrix.putVal(3, 1,0.);
        MyMatrix.putVal(4, 1,-km[1]);
        MyMatrix.putVal(5, 1,km[2]);
        
        MyMatrix.putVal(0, 2,0.);
        MyMatrix.putVal(1, 2,km[2]);
        MyMatrix.putVal(2, 2,km[3]);
        MyMatrix.putVal(3, 2,0.);
        MyMatrix.putVal(4, 2,-km[2]);
        MyMatrix.putVal(5, 2,km[4]);
        
        MyMatrix.putVal(0, 3,-km[0]);
        MyMatrix.putVal(1, 3,0.);
        MyMatrix.putVal(2, 3,0.);
        MyMatrix.putVal(3, 3,km[0]);
        MyMatrix.putVal(4, 3,0.);
        MyMatrix.putVal(5, 3,0.);
        
        MyMatrix.putVal(0, 4,0.);
        MyMatrix.putVal(1, 4,-km[1]);
        MyMatrix.putVal(2, 4,-km[2]);
        MyMatrix.putVal(3, 4,0.);
        MyMatrix.putVal(4, 4,km[1]);
        MyMatrix.putVal(5, 4,-km[2]);
                
        MyMatrix.putVal(0, 5,0.);
        MyMatrix.putVal(1, 5,km[2]);
        MyMatrix.putVal(2, 5,km[4]);
        MyMatrix.putVal(3, 5,0.);
        MyMatrix.putVal(4, 5,-km[2]);
        MyMatrix.putVal(5, 5,km[3]);
        
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
        */
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double A  =this.theCrossSection.getA();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        
        double F2;
        if(As2==0.){
            F2=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F2=(12.*E*J3)/(G*As2*L*L); // Thas it for Timoshenko beam assumption
        }
        
        
        MyMatrix.putVal(0, 0,E*A/L);
        MyMatrix.putVal(1, 0,0.);
        MyMatrix.putVal(2, 0,0.);
        MyMatrix.putVal(3, 0,-E*A/L);
        MyMatrix.putVal(4, 0,0.);
        MyMatrix.putVal(5, 0,0.);
        
        MyMatrix.putVal(0, 1,0.);
        MyMatrix.putVal(1, 1,12.*E*J3/(L*L*L*(1+F2)));
        MyMatrix.putVal(2, 1,6.*E*J3/(L*L*(1+F2)));
        MyMatrix.putVal(3, 1,0.);
        MyMatrix.putVal(4, 1,-12.*E*J3/(L*L*L*(1+F2)));
        MyMatrix.putVal(5, 1,6.*E*J3/(L*L*(1+F2)));
        
        MyMatrix.putVal(0, 2,0.);
        MyMatrix.putVal(1, 2,6.*E*J3/(L*L*(1+F2)));
        MyMatrix.putVal(2, 2,(4.+F2)*E*J3/(L*(1+F2)));
        MyMatrix.putVal(3, 2,0.);
        MyMatrix.putVal(4, 2,-6.*E*J3/(L*L*(1+F2)));
        MyMatrix.putVal(5, 2,(2.-F2)*E*J3/(L*(1+F2)));
        
        MyMatrix.putVal(0, 3,-E*A/L);
        MyMatrix.putVal(1, 3,0.);
        MyMatrix.putVal(2, 3,0.);
        MyMatrix.putVal(3, 3,E*A/L);
        MyMatrix.putVal(4, 3,0.);
        MyMatrix.putVal(5, 3,0.);
        
        MyMatrix.putVal(0, 4,0.);
        MyMatrix.putVal(1, 4,-12.*E*J3/(L*L*L*(1+F2)));
        MyMatrix.putVal(2, 4,-6.*E*J3/(L*L*(1+F2)));
        MyMatrix.putVal(3, 4,0.);
        MyMatrix.putVal(4, 4,12.*E*J3/(L*L*L*(1+F2)));
        MyMatrix.putVal(5, 4,-6.*E*J3/(L*L*(1+F2)));
                
        MyMatrix.putVal(0, 5,0.);
        MyMatrix.putVal(1, 5,6.*E*J3/(L*L*(1+F2)));
        MyMatrix.putVal(2, 5,(2.-F2)*E*J3/(L*(1+F2)));
        MyMatrix.putVal(3, 5,0.);
        MyMatrix.putVal(4, 5,-6.*E*J3/(L*L*(1+F2)));
        MyMatrix.putVal(5, 5,(4.+F2)*E*J3/(L*(1+F2)));

        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }

    private AbstractMatrix getM(int consistent) {
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
        double c=m*L/((1+F2)*(1+F2));
        double[] mcoefs = new double[10];
        mcoefs[0]=13./35.+7*F2/10.+F2*F2/3.;  mcoefs[1]=11*L/210.+(11*F2/120.+F2*F2/24.)*L;
        mcoefs[2]=9/70.+3*F2/10.+F2*F2/6.; mcoefs[3]=13*L/420.+(3*F2/40.+F2*F2/24.)*L;
        mcoefs[4]=L*L/105.+(F2/60.+F2*F2/120.)*L*L; mcoefs[5]=L*L/140+(F2/60.+F2*F2/120.)*L*L;
        mcoefs[6]=6./5.; mcoefs[7]=(1./10.-F2/2.)*L;
        mcoefs[8]=(2./15.+F2/6.+F2*F2/3.)*L*L; mcoefs[9]=(-1./30-F2/6.+F2*F2/6.)*L*L;
        switch(consistent){
            case 0:
                MyMatrix.putVal(0, 0,m/3);
                MyMatrix.putVal(1, 0,0.);
                MyMatrix.putVal(2, 0,0.);
                MyMatrix.putVal(3, 0,m/6);
                MyMatrix.putVal(4, 0,0.);
                MyMatrix.putVal(5, 0,0.);
                
                MyMatrix.putVal(0,1,0);
                MyMatrix.putVal(1,1,c*(mcoefs[0]+mcoefs[6]*r*r/(L*L)));
                MyMatrix.putVal(2,1,c*(mcoefs[1]+mcoefs[7]*r*r/(L*L)));
                MyMatrix.putVal(3,1,0.);
                MyMatrix.putVal(4,1,c*(mcoefs[2]-mcoefs[6]*r*r/(L*L)));
                MyMatrix.putVal(5,1,c*(-mcoefs[3]+mcoefs[7]*r*r/(L*L)));
                
                MyMatrix.putVal(0,2,0.);
                MyMatrix.putVal(1,2,c*(mcoefs[1]+mcoefs[7]*r*r/(L*L)));
                MyMatrix.putVal(2,2,c*(mcoefs[4]+mcoefs[8]*r*r/(L*L)));
                MyMatrix.putVal(3,2,0.);
                MyMatrix.putVal(4,2,c*(mcoefs[3]-mcoefs[7]*r*r/(L*L)));
                MyMatrix.putVal(5,2,c*(-mcoefs[5]+mcoefs[9]*r*r/(L*L)));
                
                MyMatrix.putVal(0,3,m/6);
                MyMatrix.putVal(1,3,0.);
                MyMatrix.putVal(2,3,0.);
                MyMatrix.putVal(3,3,m/3);
                MyMatrix.putVal(4,3,0.);
                MyMatrix.putVal(5,3,0.);
                
                MyMatrix.putVal(0,4,0.);
                MyMatrix.putVal(1,4,c*(mcoefs[2]-mcoefs[6]*r*r/(L*L)));
                MyMatrix.putVal(2,4,c*(mcoefs[3]-mcoefs[7]*r*r/(L*L)));
                MyMatrix.putVal(3,4,0.);
                MyMatrix.putVal(4,4,c*(mcoefs[0]+mcoefs[6]*r*r/(L*L)));
                MyMatrix.putVal(5,4,c*(-mcoefs[1]-mcoefs[7]*r*r/(L*L)));
                
                MyMatrix.putVal(0,5,0.);
                MyMatrix.putVal(1,5,c*(-mcoefs[3]+mcoefs[7]*r*r/(L*L)));
                MyMatrix.putVal(2,5,c*(-mcoefs[5]+mcoefs[9]*r*r/(L*L)));
                MyMatrix.putVal(3,5,0.);
                MyMatrix.putVal(4,5,c*(-mcoefs[1]-mcoefs[7]*r*r/(L*L)));
                MyMatrix.putVal(5,5,c*(mcoefs[4]+mcoefs[8]*r*r/(L*L)));
                break;
            default:
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
                MyMatrix.putVal(2,2,0.);
                MyMatrix.putVal(3,2,0.);
                MyMatrix.putVal(4,2,0.);
                MyMatrix.putVal(5,2,m*1.0e-12);
                
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
                MyMatrix.putVal(5,5,m*1.0e-12);
                break;
                
        }
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getM() {
        int cons;
        if(this.consistentM){
            cons=0;
        } else {
            cons=1;
        }
        return this.getM(cons);
    }

    @Override
    int getElementNdofs() {
        return ndofs;
    }

    @Override
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
    
    public AbstractMatrix getM_v() {
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
        // from global to local disps
        disps=this.getT().times(disps);
        
        //compute total strain
        this.Etotal[0]=(-disps.get(0, 0)+disps.get(3, 0))/this.getL();
        this.Etotal[2]=(-disps.get(2, 0)+disps.get(5, 0))/this.getL();
        this.Etotal[1]=(-disps.get(1, 0)+disps.get(4, 0))/this.getL()
                      -(disps.get(5, 0)+disps.get(2, 0))/2.;
        // find type (class) of material (elastoplastic or elastic_ simple material)
        if(this.ElementMaterial.getType()==0){
            this.Stress[0]=((ElasticMaterial) this.ElementMaterial).getElasticity()*this.Etotal[0]*this.theCrossSection.getA();
            this.Stress[1]=((ElasticMaterial) this.ElementMaterial).getG()*this.Etotal[1]*this.theCrossSection.getAs2();
            this.Stress[2]=((ElasticMaterial) this.ElementMaterial).getElasticity()*this.Etotal[2]*this.theCrossSection.getJ3();
        }else{
            double kin=((ElastoPlastic) this.ElementMaterial).getKkinematic();
            double iso=((ElastoPlastic) this.ElementMaterial).getKisotropic();
            double E=((ElastoPlastic) this.ElementMaterial).getElasticity();
            double G=((ElastoPlastic) this.ElementMaterial).getG();
            double A=this.theCrossSection.getA();
            double As=this.theCrossSection.getAs2();
            double J=this.theCrossSection.getJ3();
            double Py=((ElastoPlastic) this.ElementMaterial).getYieldStress()[0];
            double My=((ElastoPlastic) this.ElementMaterial).getYieldStress()[1];
            double c=((BeamElastoPlastic) this.ElementMaterial).getC();
            double Dg;
            double DgIncr;
            double normRes;
            AbstractMatrix res  = new AbstractMatrix(5,1);
            AbstractMatrix vec  = new AbstractMatrix(5,1);
            AbstractMatrix Dint = new AbstractMatrix(5,1);
            AbstractMatrix Amat = new AbstractMatrix(5,5);
            AbstractMatrix Dmat = new AbstractMatrix(5,5);
            
            // trial state
            plast=false;
            
            EplasticTrial[0]=Eplastic[0]; EplasticTrial[1]=Eplastic[1]; EplasticTrial[2]=Eplastic[2];
            this.aParTrial=this.aPar;
            this.qPar_NTrial=this.qPar_N;
            this.qPar_MTrial=this.qPar_M;
            
            Dg=0.;
            
            this.Stress[0]=E*A*(this.Etotal[0]-this.EplasticTrial[0]);
            this.Stress[1]=G*As*(this.Etotal[1]-this.EplasticTrial[1]);
            this.Stress[2]=E*J*(this.Etotal[2]-this.EplasticTrial[2]);

            this.Xsi[0]=this.Stress[0]-kin*qPar_NTrial;
            this.Xsi[1]=this.Stress[1];
            this.Xsi[2]=this.Stress[2]-kin*qPar_MTrial;
            
            this.YieldFunction=(Xsi[0]/Py)*(Xsi[0]/Py)+(Xsi[2]/My)*(Xsi[2]/My)
                              +c*(Xsi[0]/Py)*(Xsi[2]/My)*(Xsi[0]/Py)*(Xsi[2]/My)
                              -(1.-iso*this.aParTrial/Py)*(1.-iso*this.aParTrial/Py);
            
            res.putVal(0, 0, Eplastic[0]-EplasticTrial[0]+Dg*( 2.*Xsi[0]/(Py*Py) + 2.*c*Xsi[0]/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)  ));
            res.putVal(1, 0, Eplastic[2]-EplasticTrial[2]+Dg*( 2.*Xsi[2]/(My*My) + 2.*c*Xsi[2]/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
            res.putVal(2, 0, qPar_N-qPar_NTrial+Dg*( 2.*Xsi[0]/(Py*Py) + 2.*c*Xsi[0]/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)  ));
            res.putVal(3, 0, qPar_M-qPar_MTrial+Dg*( 2.*Xsi[2]/(My*My) + 2.*c*Xsi[2]/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
            res.putVal(4, 0, aPar-aParTrial-Dg*2./Py*( 1.-iso*aParTrial/Py ));
                    
            normRes=res.norm1();
            int jj=0;
            while((YieldFunction>0. || normRes>1e-8)&&jj<100){
                ++jj;
                this.plast=true;
                Dmat.init();
                Amat.init();
                
                Amat.putVal(0, 0, 1./(E*A)+Dg*(2./(Py*Py)+2.*c/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)) );
                Amat.putVal(0, 1, Dg*(4.*c/(Py*Py)*Xsi[0]*Xsi[2]/(My*My)) );
                Amat.putVal(0, 2, Dg*(2./(Py*Py)+2.*c/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)));
                Amat.putVal(0, 3, Dg*(4.*c*Xsi[0]/(Py*Py)*Xsi[2]/(My*My)  ));
                
                Amat.putVal(1, 0, Dg*(4.*c/(Py*Py)*Xsi[0]*Xsi[2]/(My*My)) );
                Amat.putVal(1, 1, 1./(E*J)+Dg*( 2./(My*My) + 2.*c/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                Amat.putVal(1, 2, Dg*(4.*c*Xsi[0]/(Py*Py)*Xsi[2]/(My*My)  ));
                Amat.putVal(1, 3, Dg*( 2./(My*My) + 2.*c/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                
                Amat.putVal(2, 0, Dg*(2./(Py*Py)+2.*c/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)));
                Amat.putVal(2, 1, Dg*( 4.*c*Xsi[0]/(Py*Py)*Xsi[2]/(My*My)  ));
                Amat.putVal(2, 2, 1/kin+Dg*(2./(Py*Py)+2.*c/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)));
                Amat.putVal(2, 3, Dg*( 4.*c*Xsi[0]/(Py*Py)*Xsi[2]/(My*My)  ));
                
                Amat.putVal(3, 0, Dg*( 4.*c*Xsi[0]/(Py*Py)*Xsi[2]/(My*My)  ));
                Amat.putVal(3, 1, Dg*( 2./(My*My) + 2.*c/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                Amat.putVal(3, 2, Dg*( 4.*c*Xsi[0]/(Py*Py)*Xsi[2]/(My*My)  ));
                Amat.putVal(3, 3, 1./kin+Dg*( 2./(My*My) + 2.*c/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
               
                Amat.putVal(4, 4, 1/iso-Dg*2./(Py*Py));
               
                Amat=Amat.inverse();
                
                vec.putVal(0, 0, ( 2.*Xsi[0]/(Py*Py) + 2.*c*Xsi[0]/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)  ));
                vec.putVal(1, 0, ( 2.*Xsi[2]/(My*My) + 2.*c*Xsi[2]/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                vec.putVal(2, 0, ( 2.*Xsi[0]/(Py*Py) + 2.*c*Xsi[0]/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)  ));
                vec.putVal(3, 0, ( 2.*Xsi[2]/(My*My) + 2.*c*Xsi[2]/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                vec.putVal(4, 0, -2./Py*( 1.-iso*aParTrial/Py ));
                
                double Numerator=YieldFunction-(((vec.transpose()).times(Amat)).times(res)).get(0, 0);
                double Denominator=(((vec.transpose()).times(Amat)).times(vec)).get(0, 0);
                DgIncr=Numerator/Denominator;
                
                Dmat.putVal(0, 0, 1./(E*A));
                Dmat.putVal(1, 1, 1./(E*J));
                Dmat.putVal(2, 2, 1./kin);
                Dmat.putVal(3, 3, 1./kin);
                Dmat.putVal(4, 4, 1./iso);
                
                Dint=(Dmat.times(Amat)).times(res.plus(vec.times(DgIncr)));
                
                EplasticTrial[0]+=Dint.get(0, 0); 
                EplasticTrial[2]+=Dint.get(1, 0);

                qPar_NTrial+=Dint.get(2, 0);
                qPar_MTrial+=Dint.get(3, 0);
                aParTrial+=Dint.get(4, 0);

                Dg+=DgIncr;
                
                this.Stress[0]=E*A*(this.Etotal[0]-this.EplasticTrial[0]);
                this.Stress[1]=G*As*(this.Etotal[1]-this.EplasticTrial[1]);
                this.Stress[2]=E*J*(this.Etotal[2]-this.EplasticTrial[2]);

                this.Xsi[0]=this.Stress[0]-kin*qPar_NTrial;
                this.Xsi[1]=this.Stress[1];
                this.Xsi[2]=this.Stress[2]-kin*qPar_MTrial;

                this.YieldFunction=(Xsi[0]/Py)*(Xsi[0]/Py)+(Xsi[2]/My)*(Xsi[2]/My)
                                  +c*(Xsi[0]/Py)*(Xsi[2]/My)*(Xsi[0]/Py)*(Xsi[2]/My)
                                  -(1.-iso*this.aParTrial/Py)*(1.-iso*this.aParTrial/Py);

                res.putVal(0, 0, Eplastic[0]-EplasticTrial[0]+Dg*( 2.*Xsi[0]/(Py*Py) + 2.*c*Xsi[0]/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)  ));
                res.putVal(1, 0, Eplastic[2]-EplasticTrial[2]+Dg*( 2.*Xsi[2]/(My*My) + 2.*c*Xsi[2]/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                res.putVal(2, 0, qPar_N-qPar_NTrial+Dg*( 2.*Xsi[0]/(Py*Py) + 2.*c*Xsi[0]/(Py*Py)*Xsi[2]*Xsi[2]/(My*My)  ));
                res.putVal(3, 0, qPar_M-qPar_MTrial+Dg*( 2.*Xsi[2]/(My*My) + 2.*c*Xsi[2]/(My*My)*Xsi[0]*Xsi[0]/(Py*Py)  ));
                res.putVal(4, 0, aPar-aParTrial-Dg*2./Py*( 1.-iso*aParTrial/Py ));

                normRes=res.norm1();
            }
            
        }
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint.init();
        Fint.addVal(0, 0, -Stress[0]);
        Fint.addVal(1, 0, -Stress[1]);
        Fint.addVal(2, 0, -Stress[2]-Stress[1]*getL()/2.);
        
        Fint.addVal(3, 0, Stress[0]);
        Fint.addVal(4, 0, Stress[1]);
        Fint.addVal(5, 0, Stress[2]-Stress[1]*getL()/2.);
        Fint=this.getT().transpose().times(Fint);
        return Fint;
    }

    @Override
    void clear() {
        this.Eplastic[0]=0.; this.Eplastic[1]=0.; this.Eplastic[2]=0.;
        this.Etotal[0]=0.; this.Etotal[1]=0.; this.Etotal[2]=0.;
        this.Stress[0]=0.; this.Stress[1]=0.; this.Stress[2]=0.;
        this.Xsi[0]=0.; this.Xsi[1]=0.; this.Xsi[0]=0.;
        this.YieldFunction=0.;
        this.aPar=0.;
        this.qPar_N=0.;
        this.qPar_M=0.;
        this.EplasticTrial[0]=0.; this.EplasticTrial[1]=0.; this.EplasticTrial[2]=0.;
        this.plast=false;
    }
    
    @Override
    void commit() {
        Eplastic[0]=EplasticTrial[0]; Eplastic[1]=EplasticTrial[1]; Eplastic[2]=EplasticTrial[2];
        this.aPar=this.aParTrial;
        this.qPar_N=this.qPar_NTrial;
        this.qPar_M=this.qPar_MTrial;
    }

    @Override
    public boolean ElementPlastified() {
        return plast;
    }

    @Override
    public double getuKu() {
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
        disps=this.getT().times(disps);
        return disps.transpose().times(getK().times(disps)).get(0, 0);
    }
    
    @Override
    public double getvMv() {
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
        velcs=this.getT().times(velcs);
        return velcs.transpose().times(getM().times(velcs)).get(0, 0);
    }

//    @Override
//    public Node getNodeHierarchy(int nodeID) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }

    @Override
    public double[] getBVNorm(double coef) {
        
        double[] norm= new double[2];
        // to be fixed
        norm[0]=0.; norm[1]=0.;
//        double abs_u1=Math.sqrt(this.getNodeHierarchy(1).getDisp(1)*this.getNodeHierarchy(1).getDisp(1)+this.getNodeHierarchy(1).getDisp(2)*this.getNodeHierarchy(1).getDisp(2))/coef;
//        double abs_u2=Math.sqrt(this.getNodeHierarchy(2).getDisp(1)*this.getNodeHierarchy(2).getDisp(1)+this.getNodeHierarchy(2).getDisp(2)*this.getNodeHierarchy(2).getDisp(2))/coef;
//        norm[0]=0.5*(abs_u1+abs_u2)*this.getL();
        return norm;
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
}
