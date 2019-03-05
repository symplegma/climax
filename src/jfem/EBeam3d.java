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
import java.util.Iterator;
import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class EBeam3d extends Element {
    private static int numberOfEBeam3d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    private double[] offset;
    public boolean RotaryInertia=false;
    
    // constructors
    public EBeam3d(int id, Node Node1, Node Node2, Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfEBeam3d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        int[] dofs = new int[6];
        dofs[0]=1; dofs[1]=2; dofs[2]=3; dofs[3]=4; dofs[4]=5; dofs[5]=6;
        Node1.setNdofs_ofNode(dofs);
        Node2.setNdofs_ofNode(dofs);
        offset = new double[6];
        for(int i=0;i<offset.length;i++)offset[i]=0.0;
        ndofs = 12;
        dimension=1;
    }

// methods
    public int getNumberOfEBeam3d() {
        return numberOfEBeam3d;
    }

    @Override
    public boolean ElementPlastified() {
        return false;
    }
    
    private double getL(){
            double x1_=0.;
            double x2_=0.;
            double y1_=0.;
            double y2_=0.;
            double z1_=0.;
            double z2_=0.;
            int j=0;
            for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
                Node theNode = it.next();
                switch(j){
                    case 0:
                        x1_=theNode.getCoords()[0];
                        y1_=theNode.getCoords()[1];
                        z1_=theNode.getCoords()[2];
                        break;
                    case 1:
                        x2_=theNode.getCoords()[0];
                        y2_=theNode.getCoords()[1];
                        z2_=theNode.getCoords()[2];
                        break;
                }
                ++j;
            }
            
            double L=Math.sqrt((x2_+offset[3]-x1_-offset[0])*(x2_+offset[3]-x1_-offset[0])
		          +(y2_+offset[4]-y1_-offset[1])*(y2_+offset[4]-y1_-offset[1])
				  +(z2_+offset[5]-z1_-offset[2])*(z2_+offset[5]-z1_-offset[2]));
            return L;
    }

    @Override
    public AbstractMatrix getK() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double A  =this.theCrossSection.getA();
	double J1 =this.theCrossSection.getJ1();
	double J2 =this.theCrossSection.getJ2();
	double J3 =this.theCrossSection.getJ3();
	double As2=this.theCrossSection.getAs2();
	double As3=this.theCrossSection.getAs3();
        
        double E  =((ElasticMaterial) this.ElementMaterial).getElasticity();
        double nu =((ElasticMaterial) this.ElementMaterial).getPoisson();
        double G  =E/(2*(1+nu));
        
        double L=this.getL();
        
        double F2;
        if(As2<=1.e-12){
            F2=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F2=(12.*E*J2)/(G*As3*L*L); // Thas it for Timoshenko beam assumption
        }
        
        double F3;
        if(As3<=1.e-12){
            F3=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F3=(12.*E*J3)/(G*As2*L*L); // Thas it for Timoshenko beam assumption
        }
        
        double KN=E*A/L;	// Full length for axial!
	double KT=G*J1/L;	// Full length for torsion!

	double K12y=12*E*J2/((1+F2)*L*L*L);
	double K12z=12*E*J3/((1+F3)*L*L*L);
	
	double K6y=6*E*J2/((1+F2)*L*L);
	double K6z=6*E*J3/((1+F3)*L*L);

	double K2y=(2-F2)*E*J2/((1+F2)*L);
	double K2z=(2-F3)*E*J3/((1+F3)*L);
	double K4y=(4+F2)*E*J2/((1+F2)*L);
	double K4z=(4+F3)*E*J3/((1+F3)*L);
        
        MyMatrix.putVal(0,0, KN);
	MyMatrix.putVal(0,6,-KN);

	MyMatrix.putVal(1,1, K12z);
	MyMatrix.putVal(1,5,K6z);
	MyMatrix.putVal(1,7,-K12z);
	MyMatrix.putVal(1,11, K6z);

	MyMatrix.putVal(2,2, K12y);
	MyMatrix.putVal(2,4,-K6y);
	MyMatrix.putVal(2,8,-K12y);
	MyMatrix.putVal(2,10,-K6y);

	MyMatrix.putVal(3,3,KT);
	MyMatrix.putVal(3,9,-KT);

	MyMatrix.putVal(4,2,-K6y);
	MyMatrix.putVal(4,4,K4y);
	MyMatrix.putVal(4,8, K6y);
	MyMatrix.putVal(4,10,K2y);

	MyMatrix.putVal(5,1, K6z);
	MyMatrix.putVal(5,5, K4z);
	MyMatrix.putVal(5,7,-K6z);
	MyMatrix.putVal(5,11, K2z);

	MyMatrix.putVal(6,0,-KN);
	MyMatrix.putVal(6,6, KN);

	MyMatrix.putVal(7,1,-K12z);
	MyMatrix.putVal(7,5,-K6z);
	MyMatrix.putVal(7,7, K12z);
	MyMatrix.putVal(7,11,-K6z);

	MyMatrix.putVal(8,2,-K12y);
	MyMatrix.putVal(8,4, K6y);
	MyMatrix.putVal(8,8, K12y);
	MyMatrix.putVal(8,10, K6y);

	MyMatrix.putVal(9,3,-KT);
	MyMatrix.putVal(9,9, KT);

	MyMatrix.putVal(10,2,-K6y);
	MyMatrix.putVal(10,4, K2y);
	MyMatrix.putVal(10,8, K6y);
	MyMatrix.putVal(10,10,K4y);

	MyMatrix.putVal(11,1, K6z);
	MyMatrix.putVal(11,5, K2z);
	MyMatrix.putVal(11,7,-K6z);
	MyMatrix.putVal(11,11,K4z);

        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getT() {
            TrMatrix = new AbstractMatrix(ndofs,ndofs);
            double L=this.getL();
            double x1=0.;
            double x2=0.;
            double y1=0.;
            double y2=0.;
            double z1=0.;
            double z2=0.;
            int j=0;
            for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
                Node theNode = it.next();
                switch(j){
                    case 0:
                        x1=theNode.getCoords()[0];
                        y1=theNode.getCoords()[1];
                        z1=theNode.getCoords()[2];
                        break;
                    case 1:
                        x2=theNode.getCoords()[0];
                        y2=theNode.getCoords()[1];
                        z2=theNode.getCoords()[2];
                        break;
                }
                ++j;
            }
            
            double cx = x2 - x1 + offset[3] - offset[0];
            double cy = y2 - y1 + offset[4] - offset[1];
            double cz = z2 - z1 + offset[5] - offset[2];
            double L1=Math.sqrt( cx * cx + cy * cy );
            
            double dx = cx;
            double dy = cy;
            double dz = cz;
            cx = cx/L; 
            cy = cy/L;
            cz = cz/L;
            double alpha=this.theCrossSection.getAngle();
            double csa = Math.cos(alpha);
            double sna = Math.sin(alpha);
            
            double xol1,yol1;
            if(L1>1e-10){
		xol1=dx/L1;
		yol1=dy/L1;
            }else{
		xol1=1.0;
		yol1=0.0;
            }

            double T00=  cx;
            double T01=  cy;
            double T02=  cz;
            double T10= -yol1 * csa - xol1 * dz * sna / L;
            double T11= +xol1 * csa - yol1 * dz * sna / L;
            double T12=  L1 * sna / L;
            double T20=  yol1 * sna - xol1 * dz * csa / L;
            double T21= -xol1 * sna - yol1 * dz * csa / L;
            double T22= +L1 * csa / L;
            
            for(int i=0;i<4;i++){
		TrMatrix.putVal(0+3*i,0+3*i,T00);
		TrMatrix.putVal(0+3*i,1+3*i,T01);
		TrMatrix.putVal(0+3*i,2+3*i,T02);
		TrMatrix.putVal(1+3*i,0+3*i,T10);
		TrMatrix.putVal(1+3*i,1+3*i,T11);
		TrMatrix.putVal(1+3*i,2+3*i,T12);
		TrMatrix.putVal(2+3*i,0+3*i,T20);
		TrMatrix.putVal(2+3*i,1+3*i,T21);
		TrMatrix.putVal(2+3*i,2+3*i,T22);
            }
            
            AbstractMatrix Le = new AbstractMatrix(12,12,0.0);
            for(int i=0;i<12;i++) Le.putVal(i,i,1.);

            Le.putVal(0, 4,offset[2]);	Le.putVal(0,5,-offset[1]);
            Le.putVal(1, 3,-offset[2]);	Le.putVal(1,5, offset[0]);
            Le.putVal(2, 3, offset[1]);	Le.putVal(2,4,-offset[0]);
            Le.putVal(6,10, offset[5]);	Le.putVal(6,11,-offset[4]);
            Le.putVal(7, 9,-offset[5]);	Le.putVal(7,11, offset[3]);
            Le.putVal(8, 9, offset[4]);	Le.putVal(8,10,-offset[3]);
            
            return TrMatrix.times(Le);
    }
    
    public void setOffset(double x1,double y1, double z1,double x2,double y2, double z2){
	offset[0]=x1;
	offset[1]=y1;
	offset[2]=z1;
	offset[3]=x2;
	offset[4]=y2;
	offset[5]=z2;
	int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            switch(j){
                case 0:
                    x1=theNode.getCoords()[0];
                    y1=theNode.getCoords()[1];
                    z1=theNode.getCoords()[2];
                    break;
                case 1:
                    x2=theNode.getCoords()[0];
                    y2=theNode.getCoords()[1];
                    z2=theNode.getCoords()[2];
                    break;
            }
            ++j;
        }
    }
    
    @Override
    public AbstractMatrix getM() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs,0.0);
        double L=this.getL();
        double p=((ElasticMaterial) this.ElementMaterial).getDensity();
        double E=((ElasticMaterial) this.ElementMaterial).getElasticity();
        double G=((ElasticMaterial) this.ElementMaterial).getG();
        double A=this.theCrossSection.getA();
        double J1=this.theCrossSection.getJ1();
        
        double J3=this.theCrossSection.getJ3();
        double As2=this.theCrossSection.getAs2();
        
        double J2=this.theCrossSection.getJ2();
        double As3=this.theCrossSection.getAs3();
        
        double m=p*A*L;
        
        double F2;
        if(As2==0.){
            F2=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F2=(12.*E*J3)/(G*As2*L*L); // Thas it for Timoshenko beam assumption
        }
        double r2=Math.sqrt(J3/A);
        double c2=m/((1+F2)*(1+F2));
        double[] mcoefs2 = new double[10];
        mcoefs2[0]=13./35.+7*F2/10.+F2*F2/3.;  mcoefs2[1]=11*L/210.+(11*F2/120.+F2*F2/24.)*L;
        mcoefs2[2]=9/70.+3*F2/10.+F2*F2/6.; mcoefs2[3]=13*L/420.+(3*F2/40.+F2*F2/24.)*L;
        mcoefs2[4]=L*L/105.+(F2/60.+F2*F2/120.)*L*L; mcoefs2[5]=L*L/140+(F2/60.+F2*F2/120.)*L*L;
        mcoefs2[6]=6./5.; mcoefs2[7]=(1./10.-F2/2.)*L;
        mcoefs2[8]=(2./15.+F2/6.+F2*F2/3.)*L*L; mcoefs2[9]=(-1./30-F2/6.+F2*F2/6.)*L*L;
        double rlq2=r2*r2/(L*L); if(!RotaryInertia)rlq2=0.0;
        
        double F3;
        if(As3==0.){
            F3=0.; // Thas it for Bernoulli-Euler beam assumption
        }else {
            F3=(12.*E*J2)/(G*As3*L*L); // Thas it for Timoshenko beam assumption
        }
        double r3=Math.sqrt(J2/A);
        double c3=m/((1+F3)*(1+F3));
        double[] mcoefs3 = new double[10];
        mcoefs3[0]=13./35.+7*F3/10.+F3*F3/3.;  mcoefs3[1]=11*L/210.+(11*F3/120.+F3*F3/24.)*L;
        mcoefs3[2]=9/70.+3*F3/10.+F3*F3/6.; mcoefs3[3]=13*L/420.+(3*F3/40.+F3*F3/24.)*L;
        mcoefs3[4]=L*L/105.+(F3/60.+F3*F3/120.)*L*L; mcoefs3[5]=L*L/140+(F3/60.+F3*F3/120.)*L*L;
        mcoefs3[6]=6./5.; mcoefs3[7]=(1./10.-F3/2.)*L;
        mcoefs3[8]=(2./15.+F3/6.+F3*F3/3.)*L*L; mcoefs3[9]=(-1./30-F3/6.+F3*F3/6.)*L*L;
        double rlq3=r3*r3/(L*L); if(!RotaryInertia)rlq3=0.0;
        
        if(this.consistentM){
//            MyMatrix.putVal(0,0,(1./3.)*m);
//            MyMatrix.putVal(6,0,(1./6.)*m); MyMatrix.putVal(0,6,MyMatrix.get(6, 0));
//            
//            MyMatrix.putVal(1,1,(13./35+6.*J3/(5*A*L*L))*m);
//            MyMatrix.putVal(5,1,(11.*L/210.+J3/(10*A*L))*m);  MyMatrix.putVal(1,5,MyMatrix.get(5, 1));
//            MyMatrix.putVal(7,1,(9./70.-6.*J3/(5*A*L*L))*m);  MyMatrix.putVal(1,7,MyMatrix.get(7, 1));
//            MyMatrix.putVal(11,1,(-13.*L/420.+J3/(10*A*L))*m);MyMatrix.putVal(1,11,MyMatrix.get(11, 1));
//            
//            MyMatrix.putVal(2,2,(13./35+6.*J2/(5*A*L*L))*m);
//            MyMatrix.putVal(4,2,(-11.*L/210.-J2/(10*A*L))*m); MyMatrix.putVal(2,4,MyMatrix.get(4, 2));
//            MyMatrix.putVal(8,2,(9./70.-6.*J2/(5*A*L*L))*m);  MyMatrix.putVal(2,8,MyMatrix.get(8, 2));
//            MyMatrix.putVal(10,2,(13.*L/420.-J2/(10*A*L))*m); MyMatrix.putVal(2,10,MyMatrix.get(10, 2));
//            
//            MyMatrix.putVal(3,3,(J1/(3.*A))*m); 
//            MyMatrix.putVal(9,3,(J1/(6.*A))*m); MyMatrix.putVal(3,9,MyMatrix.get(9, 3));
//            
//            MyMatrix.putVal(4,4,(L*L/105.+2.*J2/(15*A))*m);
//            MyMatrix.putVal(8,4,(-13.*L/420.+J2/(10*A*L))*m); MyMatrix.putVal(4,8,MyMatrix.get(8, 4));
//            MyMatrix.putVal(10,4,(-L*L/140.-J2/(30*A))*m);    MyMatrix.putVal(4,10,MyMatrix.get(10, 4));
//            
//            MyMatrix.putVal(5,5,(L*L/105.+2.*J3/(15*A))*m);
//            MyMatrix.putVal(7,5,(13.*L/420.-J3/(10*A*L))*m); MyMatrix.putVal(5,7,MyMatrix.get(7, 5));
//            MyMatrix.putVal(11,5,(-L*L/140.-J3/(30*A))*m);   MyMatrix.putVal(5,11,MyMatrix.get(11, 5));
//            
//            MyMatrix.putVal(6,6,(1./3.)*m);
//            
//            MyMatrix.putVal(7,7,(13./35+6.*J3/(5*A*L*L))*m);
//            MyMatrix.putVal(11,7,(-11.*L/210.-J3/(10*A*L))*m); MyMatrix.putVal(7,11,MyMatrix.get(11, 7));
//            
//            MyMatrix.putVal(8,8,(13./35+6.*J2/(5*A*L*L))*m);
//            MyMatrix.putVal(10,8,(11.*L/210.+J2/(10*A*L))*m); MyMatrix.putVal(8,10,MyMatrix.get(10, 8));
//            
//            MyMatrix.putVal(9,9,(J1/(3.*A))*m);
//            
//            MyMatrix.putVal(10,10,(L*L/105.+2.*J2/(15*A))*m);
//            
//            MyMatrix.putVal(11,11,(L*L/105.+2.*J2/(15*A))*m);    
            
            
            MyMatrix.putVal(0, 0,m/3);
            MyMatrix.putVal(6, 0,m/6); MyMatrix.putVal(0,6,MyMatrix.get(6, 0));
            
            MyMatrix.putVal(1,1,c2*(mcoefs2[0]+mcoefs2[6]*rlq2));
            MyMatrix.putVal(5,1,c2*(mcoefs2[1]+mcoefs2[7]*rlq2)); MyMatrix.putVal(1,5,MyMatrix.get(5, 1));
            MyMatrix.putVal(7,1,c2*(mcoefs2[2]-mcoefs2[6]*rlq2)); MyMatrix.putVal(1,7,MyMatrix.get(7, 1));
            MyMatrix.putVal(11,1,c2*(-mcoefs2[3]+mcoefs2[7]*rlq2)); MyMatrix.putVal(1,11,MyMatrix.get(11, 1));
            
            MyMatrix.putVal(2,2,c3*(mcoefs3[0]+mcoefs3[6]*rlq3));
            MyMatrix.putVal(4,2,c3*(-mcoefs3[1]-mcoefs3[7]*rlq3)); MyMatrix.putVal(2,4,MyMatrix.get(4, 2));
            MyMatrix.putVal(8,2,c3*(mcoefs3[2]-mcoefs3[6]*rlq3));  MyMatrix.putVal(2,8,MyMatrix.get(8, 2));
            MyMatrix.putVal(10,2,c3*(mcoefs3[3]-mcoefs3[7]*rlq3)); MyMatrix.putVal(2,10,MyMatrix.get(10, 2));
            
            MyMatrix.putVal(3, 3,m*J1/(3*A));
            MyMatrix.putVal(9, 3,m*J1/(6*A)); MyMatrix.putVal(3,9,MyMatrix.get(9, 3));

            MyMatrix.putVal(4,4,c3*(mcoefs3[4]+mcoefs3[8]*rlq3));
            MyMatrix.putVal(8,4,c3*(-mcoefs3[3]+mcoefs3[7]*rlq3)); MyMatrix.putVal(4,8,MyMatrix.get(8, 4));
            MyMatrix.putVal(10,4,c3*(-mcoefs3[5]+mcoefs3[9]*rlq3)); MyMatrix.putVal(4,10,MyMatrix.get(10, 4));

            MyMatrix.putVal(5,5,c2*(mcoefs2[4]+mcoefs2[8]*rlq2));
            MyMatrix.putVal(7,5,c2*(mcoefs2[3]-mcoefs2[7]*rlq2)); MyMatrix.putVal(5,7,MyMatrix.get(7, 5));
            MyMatrix.putVal(11,5,c2*(-mcoefs2[5]+mcoefs2[9]*rlq2));MyMatrix.putVal(5,11,MyMatrix.get(11, 5));

            MyMatrix.putVal(6,6,m/3);

            MyMatrix.putVal(7,7,c2*(mcoefs2[0]+mcoefs2[6]*rlq2));
            MyMatrix.putVal(11,7,c2*(-mcoefs2[1]-mcoefs2[7]*rlq2)); MyMatrix.putVal(7,11,MyMatrix.get(11, 7));
            
            MyMatrix.putVal(8,8,c3*(mcoefs3[0]+mcoefs3[6]*rlq3));
            MyMatrix.putVal(10,8,c3*(mcoefs3[1]+mcoefs3[7]*rlq3)); MyMatrix.putVal(8,10,MyMatrix.get(10, 8));
            
            MyMatrix.putVal(9, 9,m*J1/(3*A));
            
            MyMatrix.putVal(10,10,c3*(mcoefs3[4]+mcoefs3[8]*rlq3));

            MyMatrix.putVal(11,11,c2*(mcoefs2[4]+mcoefs2[8]*rlq2));
            
        }else{
            MyMatrix.putVal(0,0,0.5*m);
            MyMatrix.putVal(1,1,0.5*m);
            MyMatrix.putVal(2,2,0.5*m);
            MyMatrix.putVal(6,6,0.5*m);
            MyMatrix.putVal(7,7,0.5*m);
            MyMatrix.putVal(8,8,0.5*m);

            MyMatrix.putVal(3,3,0.5*m*1.e-12);
            MyMatrix.putVal(4,4,0.5*m*1.e-12);
            MyMatrix.putVal(5,5,0.5*m*1.e-12);
            MyMatrix.putVal(9,9,0.5*m*1.e-12);
            MyMatrix.putVal(10,10,0.5*m*1.e-12);
            MyMatrix.putVal(11,11,0.5*m*1.e-12);
        }
        
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }

    @Override
    public int getElementNdofs() {
        return ndofs;
    }

    @Override
    int[] getFtable() {
        int[] theEFTable = new int[ndofs];
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theEFTable[j]=theNode.getFtable()[0];
            theEFTable[j+1]=theNode.getFtable()[1];
            theEFTable[j+2]=theNode.getFtable()[2];
            theEFTable[j+3]=theNode.getFtable()[3];
            theEFTable[j+4]=theNode.getFtable()[4];
            theEFTable[j+5]=theNode.getFtable()[5];
            j=j+6;
        }
        return theEFTable;
    }

    @Override
    public AbstractMatrix getF() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[2]);
            disps.addVal(i+3, 0, theNode.getDispsTrial()[3]);
            disps.addVal(i+4, 0, theNode.getDispsTrial()[4]);
            disps.addVal(i+5, 0, theNode.getDispsTrial()[5]);
            i+=6;
        }
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint=this.getK().times(disps);
        return Fint;
    }

    @Override
    void clear(){}

    @Override
    AbstractMatrix getM_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i  , 0, theNode.getVelcsTrial()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
            velcs.addVal(i+2, 0, theNode.getVelcsTrial()[2]);
            velcs.addVal(i+3, 0, theNode.getVelcsTrial()[3]);
            velcs.addVal(i+4, 0, theNode.getVelcsTrial()[4]);
            velcs.addVal(i+5, 0, theNode.getVelcsTrial()[5]);
            i+=6;
        }
        return getM().times(velcs);
    }

    @Override
    AbstractMatrix getM_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i  , 0, theNode.getAcclsTrial()[0]);
            accls.addVal(i+1, 0, theNode.getAcclsTrial()[1]);
            accls.addVal(i+2, 0, theNode.getAcclsTrial()[2]);
            accls.addVal(i+3, 0, theNode.getAcclsTrial()[3]);
            accls.addVal(i+4, 0, theNode.getAcclsTrial()[4]);
            accls.addVal(i+5, 0, theNode.getAcclsTrial()[5]);
            i+=6;
        }
        return getM().times(accls);
    }

    @Override
    AbstractMatrix getM_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[2]);
            disps.addVal(i+3, 0, theNode.getDispsTrial()[3]);
            disps.addVal(i+4, 0, theNode.getDispsTrial()[4]);
            disps.addVal(i+5, 0, theNode.getDispsTrial()[5]);
            i+=6;
        }
        return getM().times(disps);
    }

    @Override
    AbstractMatrix getK_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            disps.addVal(i+2, 0, theNode.getDispsTrial()[2]);
            disps.addVal(i+3, 0, theNode.getDispsTrial()[3]);
            disps.addVal(i+4, 0, theNode.getDispsTrial()[4]);
            disps.addVal(i+5, 0, theNode.getDispsTrial()[5]);
            i+=6;
        }
        return getK().times(disps);
    }
    
    @Override
    AbstractMatrix getK_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i  , 0, theNode.getVelcsTrial()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
            velcs.addVal(i+2, 0, theNode.getVelcsTrial()[2]);
            velcs.addVal(i+3, 0, theNode.getVelcsTrial()[3]);
            velcs.addVal(i+4, 0, theNode.getVelcsTrial()[4]);
            velcs.addVal(i+5, 0, theNode.getVelcsTrial()[5]);
            i+=6;
        }
        return getK().times(velcs);
    }

    @Override
    AbstractMatrix getK_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i  , 0, theNode.getAcclsTrial()[0]);
            accls.addVal(i+1, 0, theNode.getAcclsTrial()[1]);
            accls.addVal(i+2, 0, theNode.getAcclsTrial()[2]);
            accls.addVal(i+3, 0, theNode.getAcclsTrial()[3]);
            accls.addVal(i+4, 0, theNode.getAcclsTrial()[4]);
            accls.addVal(i+5, 0, theNode.getAcclsTrial()[5]);
            i+=6;
        }
        return getK().times(accls);
    }

    @Override
    void commit(){}

    @Override
    public double getuKu() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDisp(1));
            disps.addVal(i+1, 0, theNode.getDisp(2));
            disps.addVal(i+2, 0, theNode.getDisp(3));
            disps.addVal(i+3, 0, theNode.getDisp(4));
            disps.addVal(i+4, 0, theNode.getDisp(5));
            disps.addVal(i+5, 0, theNode.getDisp(6));
            i+=6;
        }
        disps=this.getT().times(disps);
        return disps.transpose().times(getK().times(disps)).get(0, 0);
    }

    @Override
    public double getvMv() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelc(1));
            velcs.addVal(i+1, 0, theNode.getVelc(2));
            velcs.addVal(i+2, 0, theNode.getVelc(3));
            velcs.addVal(i+3, 0, theNode.getVelc(4));
            velcs.addVal(i+4, 0, theNode.getVelc(5));
            velcs.addVal(i+5, 0, theNode.getVelc(6));
            i+=6;
        }
        velcs=this.getT().times(velcs);
        return velcs.transpose().times(getM().times(velcs)).get(0, 0);
    }

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getVolume(){
        return this.getL()*this.theCrossSection.getA();
    }
    
    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public void setRotaryInertiaContrinution(boolean bool){this.RotaryInertia=bool;}
}
