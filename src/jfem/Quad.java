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
abstract public class Quad extends Element{
    protected static AbstractMatrix MyMatrix = new AbstractMatrix(8,8);
    private static AbstractMatrix Matrix_2_4 = new AbstractMatrix(2,4);
    private static AbstractMatrix Matrix_4_2 = new AbstractMatrix(4,2);
    private static AbstractMatrix Matrix_3_4 = new AbstractMatrix(3,4);
    private static AbstractMatrix Matrix_4_8 = new AbstractMatrix(4,8);
    private static AbstractMatrix Matrix_2_2 = new AbstractMatrix(2,2);
    private static AbstractMatrix Matrix_3_8 = new AbstractMatrix(3,8);
    protected static AbstractMatrix Matrix_3_3 = new AbstractMatrix(3,3);
    protected final static int GaussPoints=2; // per direction of integration
    
    // default constructor
    public Quad(){
        ndofs = 8;
        dimension=2;
    }
    
    public AbstractMatrix getT() {
        MyMatrix=new AbstractMatrix(ndofs, ndofs);
        MyMatrix.init(0.0);
        for(int i=0;i<ndofs;i++){MyMatrix.set(i, i, 1.0);}
        return MyMatrix;
    }
    
    protected double getArea(){
        double Area=0.;
        for(int i=1; i<=GaussPoints; i++){
            for(int j=1; j<=GaussPoints; j++){
                Area+=this.getDetJ(this.getGaussCoord(i), this.getGaussCoord(j))
                        *this.getGaussWeight(i)*this.getGaussWeight(j);
            }
        }
        return Area;
    }
    
    public double getVolume(){
        return this.getArea()*this.theCrossSection.getThickness();
    }
    
    protected AbstractMatrix getB1(double xsi, double eta){
        double detJ=this.getArea();
        for(int i=0; i<3; i++){
            for(int j=0; j<4; j++){
                Matrix_3_4.putVal(i, j, 0.0);
            }
        }
        Matrix_2_2=this.getJ(xsi, eta);
        double J11,J12,J21,J22;
        J11=Matrix_2_2.get(0, 0);
        J12=Matrix_2_2.get(0, 1);
        J21=Matrix_2_2.get(1, 0);
        J22=Matrix_2_2.get(1, 1);
        
        Matrix_3_4.putVal(0, 0, J22);
        Matrix_3_4.putVal(0, 1, -J12);
        
        Matrix_3_4.putVal(1, 2, -J21);
        Matrix_3_4.putVal(1, 3, J11);
        
        Matrix_3_4.putVal(2, 0, -J21);
        Matrix_3_4.putVal(2, 1, J11);
        Matrix_3_4.putVal(2, 2, J22);
        Matrix_3_4.putVal(2, 3, -J12);
        return Matrix_3_4.times(1./detJ);
    }
    
    protected AbstractMatrix getB2(double xsi, double eta){
        for(int i=0; i<4; i++){
            for(int j=0; j<8; j++){
                Matrix_4_8.putVal(i, j, 0.0);
            }
        }
        
        Matrix_4_8.putVal(0, 0,this.ShapeFunction_xsi(1, xsi, eta) );
        Matrix_4_8.putVal(0, 2,this.ShapeFunction_xsi(2, xsi, eta) );
        Matrix_4_8.putVal(0, 4,this.ShapeFunction_xsi(3, xsi, eta) );
        Matrix_4_8.putVal(0, 6,this.ShapeFunction_xsi(4, xsi, eta) );
        
        Matrix_4_8.putVal(1, 0,this.ShapeFunction_eta(1, xsi, eta) );
        Matrix_4_8.putVal(1, 2,this.ShapeFunction_eta(2, xsi, eta) );
        Matrix_4_8.putVal(1, 4,this.ShapeFunction_eta(3, xsi, eta) );
        Matrix_4_8.putVal(1, 6,this.ShapeFunction_eta(4, xsi, eta) );
        
        Matrix_4_8.putVal(2, 1,this.ShapeFunction_xsi(1, xsi, eta) );
        Matrix_4_8.putVal(2, 3,this.ShapeFunction_xsi(2, xsi, eta) );
        Matrix_4_8.putVal(2, 5,this.ShapeFunction_xsi(3, xsi, eta) );
        Matrix_4_8.putVal(2, 7,this.ShapeFunction_xsi(4, xsi, eta) );
        
        Matrix_4_8.putVal(3, 1,this.ShapeFunction_eta(1, xsi, eta) );
        Matrix_4_8.putVal(3, 3,this.ShapeFunction_eta(2, xsi, eta) );
        Matrix_4_8.putVal(3, 5,this.ShapeFunction_eta(3, xsi, eta) );
        Matrix_4_8.putVal(3, 7,this.ShapeFunction_eta(4, xsi, eta) );
        
        return Matrix_4_8;
    }
    
    protected AbstractMatrix getB(double xsi, double eta){
        Matrix_3_8=getB1(xsi,eta).times(getB2(xsi,eta));
        return Matrix_3_8;
    }
    
    protected double getGaussCoord(int wG){
        double coord;
        switch(wG){
            case 1: coord=+Math.sqrt(3.)/3.; break;
            case 2: coord=-Math.sqrt(3.)/3.; break;
            default:coord=0.; break;
        }
        return coord;
    }
    
    protected double getGaussWeight(int wG){
        double weight;
        switch(wG){
//            case 1: weight=+1.0/4.; break;
//            case 2: weight=+1.0/4.; break;
            case 1: weight=+1.0; break;
            case 2: weight=+1.0; break;
            default:weight=0.; break;
        }
        return weight;
    }
    
    @Override
    public double ShapeFunction(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=0.25*(1-xsi)*(1-eta); break;
            case 2: N=0.25*(1+xsi)*(1-eta); break;
            case 3: N=0.25*(1+xsi)*(1+eta); break;
            case 4: N=0.25*(1-xsi)*(1+eta); break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    protected double ShapeFunction_xsi(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=-0.25*(1-eta); break;
            case 2: N= 0.25*(1-eta); break;
            case 3: N= 0.25*(1+eta); break;
            case 4: N=-0.25*(1+eta); break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    protected double ShapeFunction_eta(int wsf, double xsi, double eta){
        double N;
        switch(wsf){
            case 1: N=-0.25*(1-xsi); break;
            case 2: N=-0.25*(1+xsi); break;
            case 3: N= 0.25*(1+xsi); break;
            case 4: N= 0.25*(1-xsi); break;
            default:N=0. ;  break;
        }
        return N;
    }
    
    protected double getDetJ(double xsi, double eta){
        double detJ;
        for(int i=0; i<4; i++){
            Matrix_2_4.putVal(0, i, ShapeFunction_xsi(i+1,xsi,eta));
            Matrix_2_4.putVal(1, i, ShapeFunction_eta(i+1,xsi,eta));
        }
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            Matrix_4_2.putVal(i, 0, theNode.getCoords()[0]);
            Matrix_4_2.putVal(i, 1, theNode.getCoords()[1]);
            ++i;
        }
        detJ=(Matrix_2_4.times(Matrix_4_2)).get(0, 0)*(Matrix_2_4.times(Matrix_4_2)).get(1, 1)
            -(Matrix_2_4.times(Matrix_4_2)).get(0, 1)*(Matrix_2_4.times(Matrix_4_2)).get(1, 0);
        
        return detJ;
    }
    
    protected AbstractMatrix getJ(double xsi, double eta){
        for(int i=0; i<4; i++){
            Matrix_2_4.putVal(0, i, ShapeFunction_xsi(i+1,xsi,eta));
            Matrix_2_4.putVal(1, i, ShapeFunction_eta(i+1,xsi,eta));
        }
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            Matrix_4_2.putVal(i, 0, theNode.getCoords()[0]);
            Matrix_4_2.putVal(i, 1, theNode.getCoords()[1]);
            ++i;
        }        
        return Matrix_2_4.times(Matrix_4_2);
    }
    
    protected AbstractMatrix getM(int consistent) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double p;
        if(this.ElementMaterial.getType()==4){
            p=((AcousticMaterial) this.ElementMaterial).getMatConstant();
        }else{
            p=((ElasticMaterial) this.ElementMaterial).getDensity();
        }
        double h=this.theCrossSection.getThickness();
        double A=getArea();
        double m=p*A*h;
        MyMatrix.init();
        switch(consistent){
            case 0:
                // that is for consistent
                double x1,x2,w1,w2;
                for(int i=0; i<ndofs; i++){
                    for(int j=0; j<ndofs; j++){
                        MyMatrix.putVal(i, j, 0.0);
                    }
                }
                for (int i = 0; i < ndofs; i+=2) {
                    int sfi = 0;
                    switch(i){
                        case 0:sfi=1;break;
                        case 1:sfi=1;break;
                        case 2:sfi=2;break;
                        case 3:sfi=2;break;
                        case 4:sfi=3;break;
                        case 5:sfi=3;break;
                        case 6:sfi=4;break;
                        case 7:sfi=4;break;
                    }
                    for(int j=i; j< ndofs ; j+=2){
                        int sfj = 0;
                        for(int nx=1; nx<=GaussPoints; nx++){
                            w1=this.getGaussWeight(nx);
                            for(int ny=1; ny<=GaussPoints; ny++){
                                w2=this.getGaussWeight(ny);
                                x1=this.getGaussCoord(nx);
                                x2=this.getGaussCoord(ny);
                                double detJ=getDetJ(x1, x2);
                                switch(j){
                                    case 0:sfj=1;break;
                                    case 1:sfj=1;break;
                                    case 2:sfj=2;break;
                                    case 3:sfj=2;break;
                                    case 4:sfj=3;break;
                                    case 5:sfj=3;break;
                                    case 6:sfj=4;break;
                                    case 7:sfj=4;break;
                                }
                                MyMatrix.set(i, j, MyMatrix.get(i, j)+h*p*detJ*w1*w2*ShapeFunction(sfi, x1, x2)*ShapeFunction(sfj, x1, x2));
                            }
                        }
                        MyMatrix.set(j, i, MyMatrix.get(i, j));
                    }
                }
                
                for (int i = 1; i < ndofs; i+=2) {
                    int sfi = 0;
                    switch(i){
                        case 0:sfi=1;break;
                        case 1:sfi=1;break;
                        case 2:sfi=2;break;
                        case 3:sfi=2;break;
                        case 4:sfi=3;break;
                        case 5:sfi=3;break;
                        case 6:sfi=4;break;
                        case 7:sfi=4;break;
                    }
                    for(int j=i; j< ndofs ; j+=2){
                        int sfj = 0;
                        for(int nx=1; nx<=GaussPoints; nx++){
                            w1=this.getGaussWeight(nx);
                            for(int ny=1; ny<=GaussPoints; ny++){
                                w2=this.getGaussWeight(ny);
                                x1=this.getGaussCoord(nx);
                                x2=this.getGaussCoord(ny);
                                double detJ=getDetJ(x1, x2);
                                switch(j){
                                    case 0:sfj=1;break;
                                    case 1:sfj=1;break;
                                    case 2:sfj=2;break;
                                    case 3:sfj=2;break;
                                    case 4:sfj=3;break;
                                    case 5:sfj=3;break;
                                    case 6:sfj=4;break;
                                    case 7:sfj=4;break;
                                }
                                MyMatrix.set(i, j, MyMatrix.get(i, j)+h*p*detJ*w1*w2*ShapeFunction(sfi, x1, x2)*ShapeFunction(sfj, x1, x2));
                            }
                        }
                        MyMatrix.set(j, i, MyMatrix.get(i, j));
                    }
                }
                break;
            default:
                m=m/4.;
                for (int i = 0; i < ndofs; i++) {
                    MyMatrix.putVal(i, i,m);
                }
                break;
        }
//        MyMatrix=this.getT().transpose().times(MyMatrix);
//        return MyMatrix.times(this.getT());
        return MyMatrix;
    }

    @Override
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
    public int getElementNdofs() {
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
            j=j+2;
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
            i+=2;
        }
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        

        Fint=this.getK().times(disps);
        return Fint;
    }

    @Override
    void clear() {
    }

    @Override
    public AbstractMatrix getM_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsTrial()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
            i+=2;
        }
        // from global to local disps
        velcs=this.getT().times(velcs);
        return getM().times(velcs);
    }
    
    public AbstractMatrix getM_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i, 0, theNode.getAcclsTrial()[0]);
            accls.addVal(i+1, 0, theNode.getAcclsTrial()[1]);
            i+=2;
        }
        // from global to local disps
        accls=this.getT().times(accls);
        return getM().times(accls);
    }
    
    public AbstractMatrix getM_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            i+=2;
        }
        // from global to local disps
        disps=this.getT().times(disps);
        return getM().times(disps);
    }
    
    public AbstractMatrix getK_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
            i+=2;
        }
        // from global to local disps
        disps=this.getT().times(disps);
        return getK().times(disps);
    }
    
    @Override
    public AbstractMatrix getK_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsTrial()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsTrial()[1]);
            i+=2;
        }
        // from global to local disps
        velcs=this.getT().times(velcs);
        return getK().times(velcs);
    }
    
    public AbstractMatrix getK_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i, 0, theNode.getAcclsTrial()[0]);
            accls.addVal(i+1, 0, theNode.getAcclsTrial()[1]);
            i+=2;
        }
        // from global to local disps
        accls=this.getT().times(accls);
        return getK().times(accls);
    }
    
    @Override
    void commit() {
    }

    public boolean ElementPlastified() {
        return false;
    }

    public double getuKu() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsConvg()[0]);
            disps.addVal(i+1, 0, theNode.getDispsConvg()[1]);
            i+=2;
        }
        // from global to local disps
        disps=this.getT().times(disps);
        return disps.transpose().times(getK().times(disps)).get(0, 0);
    }
    
    public double getvMv() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsConvg()[0]);
            velcs.addVal(i+1, 0, theNode.getVelcsConvg()[1]);
            i+=2;
        }
        // from global to local disps
        velcs=this.getT().times(velcs);
        return velcs.transpose().times(getM().times(velcs)).get(0, 0);
    }
    
//    @Override
//    public Matrix getK() {
//        double x1,x2,w1,w2;
//        double h;
//        for(int i=0; i<8; i++){
//            for(int j=0; j<8; j++){
//                MyMatrix.putVal(i, j, 0.0);
//            }
//        }
//        for(int i=1; i<=GaussPoints; i++){
//            for(int j=1; j<=GaussPoints; j++){
//                x1=this.getGaussCoord(i);
//                x2=this.getGaussCoord(j);
//                w1=this.getGaussWeight(i);
//                w2=this.getGaussWeight(j);
//                h=this.theCrossSection.getThickness();
//                double detJ=getDetJ(x1, x2);
//                MyMatrix=MyMatrix.plus( 
//                 ((getB(x1, x2)).transpose()).times(getEmat()).times(getB(x1, x2)).times(h).times(detJ*16.).times(w1).times(w2)
//                 );
//            }
//        }
//        return MyMatrix;
//    }
    
    @Override
    public AbstractMatrix getK() {
        double x1,x2,w1,w2;
        double h;
        for(int i=0; i<8; i++){
            for(int j=0; j<8; j++){
                MyMatrix.putVal(i, j, 0.0);
            }
        }
        double e11,e12,e13,e22,e23,e33;
        e11 = this.getEmat().get(0,0);
        e12 = this.getEmat().get(0,1);
        e13 = this.getEmat().get(0,2);
        e22 = this.getEmat().get(1,1);
        e23 = this.getEmat().get(1,2); 
        e33 = this.getEmat().get(2,2);

        for(int i=1; i<=GaussPoints; i++){
            for(int j=1; j<=GaussPoints; j++){
                x1=this.getGaussCoord(i);
                x2=this.getGaussCoord(j);
                w1=this.getGaussWeight(i);
                w2=this.getGaussWeight(j);
                h=this.theCrossSection.getThickness();
                double detJ=getDetJ(x1, x2);
                double J11=this.getJ(x1, x2).get(0,0); double J12=this.getJ(x1, x2).get(0,1);
                double J21=this.getJ(x1, x2).get(1,0); double J22=this.getJ(x1, x2).get(1,1);
                detJ=J11*J22-J12*J21;
                double N1x=ShapeFunction_xsi(1,x1,x2); double N1y=ShapeFunction_eta(1,x1,x2);
                double N2x=ShapeFunction_xsi(2,x1,x2); double N2y=ShapeFunction_eta(2,x1,x2);
                double N3x=ShapeFunction_xsi(3,x1,x2); double N3y=ShapeFunction_eta(3,x1,x2);
                double N4x=ShapeFunction_xsi(4,x1,x2); double N4y=ShapeFunction_eta(4,x1,x2);
                double coef=h*detJ*w1*w2;
                
                MyMatrix.set(0, 0, MyMatrix.get(0, 0)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e13*((N1x*J22)/detJ-(N1y*J12)/detJ)+e33*((N1y*J11)/detJ-(N1x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e11*((N1x*J22)/detJ-(N1y*J12)/detJ)+e13*((N1y*J11)/detJ-(N1x*J21)/detJ))));
                MyMatrix.set(0, 1, MyMatrix.get(0, 1)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e33*((N1x*J22)/detJ-(N1y*J12)/detJ)+e23*((N1y*J11)/detJ-(N1x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N1x*J22)/detJ-(N1y*J12)/detJ)+e12*((N1y*J11)/detJ-(N1x*J21)/detJ))));
                MyMatrix.set(0, 2, MyMatrix.get(0, 2)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e13*((N2x*J22)/detJ-(N2y*J12)/detJ)+e33*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e11*((N2x*J22)/detJ-(N2y*J12)/detJ)+e13*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(0, 3, MyMatrix.get(0, 3)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e33*((N2x*J22)/detJ-(N2y*J12)/detJ)+e23*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N2x*J22)/detJ-(N2y*J12)/detJ)+e12*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(0, 4, MyMatrix.get(0, 4)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e33*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e11*((N3x*J22)/detJ-(N3y*J12)/detJ)+e13*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(0, 5, MyMatrix.get(0, 5)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e33*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e12*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(0, 6, MyMatrix.get(0, 6)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e11*((N4x*J22)/detJ-(N4y*J12)/detJ)+e13*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(0, 7, MyMatrix.get(0, 7)+coef*(((N1y*J11)/detJ-(N1x*J21)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e12*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                
                MyMatrix.set(1, 1, MyMatrix.get(1, 1)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e33*((N1x*J22)/detJ-(N1y*J12)/detJ)+e23*((N1y*J11)/detJ-(N1x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e23*((N1x*J22)/detJ-(N1y*J12)/detJ)+e22*((N1y*J11)/detJ-(N1x*J21)/detJ))));
                MyMatrix.set(1, 2, MyMatrix.get(1, 2)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N2x*J22)/detJ-(N2y*J12)/detJ)+e33*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e12*((N2x*J22)/detJ-(N2y*J12)/detJ)+e23*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(1, 3, MyMatrix.get(1, 3)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e33*((N2x*J22)/detJ-(N2y*J12)/detJ)+e23*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e23*((N2x*J22)/detJ-(N2y*J12)/detJ)+e22*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(1, 4, MyMatrix.get(1, 4)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e33*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e12*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(1, 5, MyMatrix.get(1, 5)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e33*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e23*((N3x*J22)/detJ-(N3y*J12)/detJ)+e22*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(1, 6, MyMatrix.get(1, 6)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e12*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(1, 7, MyMatrix.get(1, 7)+coef*(((N1x*J22)/detJ-(N1y*J12)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N1y*J11)/detJ-(N1x*J21)/detJ)*(e23*((N4x*J22)/detJ-(N4y*J12)/detJ)+e22*((N4y*J11)/detJ-(N4x*J21)/detJ))));
            
                MyMatrix.set(2, 2, MyMatrix.get(2, 2)+coef*(((N2y*J11)/detJ-(N2x*J21)/detJ)*(e13*((N2x*J22)/detJ-(N2y*J12)/detJ)+e33*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N2x*J22)/detJ-(N2y*J12)/detJ)*(e11*((N2x*J22)/detJ-(N2y*J12)/detJ)+e13*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(2, 3, MyMatrix.get(2, 3)+coef*(((N2y*J11)/detJ-(N2x*J21)/detJ)*(e33*((N2x*J22)/detJ-(N2y*J12)/detJ)+e23*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N2x*J22)/detJ-(N2y*J12)/detJ)*(e13*((N2x*J22)/detJ-(N2y*J12)/detJ)+e12*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(2, 4, MyMatrix.get(2, 4)+coef*(((N2y*J11)/detJ-(N2x*J21)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e33*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N2x*J22)/detJ-(N2y*J12)/detJ)*(e11*((N3x*J22)/detJ-(N3y*J12)/detJ)+e13*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(2, 5, MyMatrix.get(2, 5)+coef*(((N2y*J11)/detJ-(N2x*J21)/detJ)*(e33*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N2x*J22)/detJ-(N2y*J12)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e12*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(2, 6, MyMatrix.get(2, 6)+coef*(((N2y*J11)/detJ-(N2x*J21)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N2x*J22)/detJ-(N2y*J12)/detJ)*(e11*((N4x*J22)/detJ-(N4y*J12)/detJ)+e13*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(2, 7, MyMatrix.get(2, 7)+coef*(((N2y*J11)/detJ-(N2x*J21)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N2x*J22)/detJ-(N2y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e12*((N4y*J11)/detJ-(N4x*J21)/detJ))));
            
                MyMatrix.set(3, 3, MyMatrix.get(3, 3)+coef*(((N2x*J22)/detJ-(N2y*J12)/detJ)*(e33*((N2x*J22)/detJ-(N2y*J12)/detJ)+e23*((N2y*J11)/detJ-(N2x*J21)/detJ))+((N2y*J11)/detJ-(N2x*J21)/detJ)*(e23*((N2x*J22)/detJ-(N2y*J12)/detJ)+e22*((N2y*J11)/detJ-(N2x*J21)/detJ))));
                MyMatrix.set(3, 4, MyMatrix.get(3, 4)+coef*(((N2x*J22)/detJ-(N2y*J12)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e33*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N2y*J11)/detJ-(N2x*J21)/detJ)*(e12*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(3, 5, MyMatrix.get(3, 5)+coef*(((N2x*J22)/detJ-(N2y*J12)/detJ)*(e33*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N2y*J11)/detJ-(N2x*J21)/detJ)*(e23*((N3x*J22)/detJ-(N3y*J12)/detJ)+e22*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(3, 6, MyMatrix.get(3, 6)+coef*(((N2x*J22)/detJ-(N2y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N2y*J11)/detJ-(N2x*J21)/detJ)*(e12*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(3, 7, MyMatrix.get(3, 7)+coef*(((N2x*J22)/detJ-(N2y*J12)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N2y*J11)/detJ-(N2x*J21)/detJ)*(e23*((N4x*J22)/detJ-(N4y*J12)/detJ)+e22*((N4y*J11)/detJ-(N4x*J21)/detJ))));
            
                MyMatrix.set(4, 4, MyMatrix.get(4, 4)+coef*(((N3y*J11)/detJ-(N3x*J21)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e33*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N3x*J22)/detJ-(N3y*J12)/detJ)*(e11*((N3x*J22)/detJ-(N3y*J12)/detJ)+e13*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(4, 5, MyMatrix.get(4, 5)+coef*(((N3y*J11)/detJ-(N3x*J21)/detJ)*(e33*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N3x*J22)/detJ-(N3y*J12)/detJ)*(e13*((N3x*J22)/detJ-(N3y*J12)/detJ)+e12*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(4, 6, MyMatrix.get(4, 6)+coef*(((N3y*J11)/detJ-(N3x*J21)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N3x*J22)/detJ-(N3y*J12)/detJ)*(e11*((N4x*J22)/detJ-(N4y*J12)/detJ)+e13*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(4, 7, MyMatrix.get(4, 7)+coef*(((N3y*J11)/detJ-(N3x*J21)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N3x*J22)/detJ-(N3y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e12*((N4y*J11)/detJ-(N4x*J21)/detJ))));
            
                MyMatrix.set(5, 5, MyMatrix.get(5, 5)+coef*(((N3x*J22)/detJ-(N3y*J12)/detJ)*(e33*((N3x*J22)/detJ-(N3y*J12)/detJ)+e23*((N3y*J11)/detJ-(N3x*J21)/detJ))+((N3y*J11)/detJ-(N3x*J21)/detJ)*(e23*((N3x*J22)/detJ-(N3y*J12)/detJ)+e22*((N3y*J11)/detJ-(N3x*J21)/detJ))));
                MyMatrix.set(5, 6, MyMatrix.get(5, 6)+coef*(((N3x*J22)/detJ-(N3y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N3y*J11)/detJ-(N3x*J21)/detJ)*(e12*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(5, 7, MyMatrix.get(5, 7)+coef*(((N3x*J22)/detJ-(N3y*J12)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N3y*J11)/detJ-(N3x*J21)/detJ)*(e23*((N4x*J22)/detJ-(N4y*J12)/detJ)+e22*((N4y*J11)/detJ-(N4x*J21)/detJ))));
            
                MyMatrix.set(6, 6, MyMatrix.get(6, 6)+coef*(((N4y*J11)/detJ-(N4x*J21)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e33*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N4x*J22)/detJ-(N4y*J12)/detJ)*(e11*((N4x*J22)/detJ-(N4y*J12)/detJ)+e13*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                MyMatrix.set(6, 7, MyMatrix.get(6, 7)+coef*(((N4y*J11)/detJ-(N4x*J21)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N4x*J22)/detJ-(N4y*J12)/detJ)*(e13*((N4x*J22)/detJ-(N4y*J12)/detJ)+e12*((N4y*J11)/detJ-(N4x*J21)/detJ))));
            
                MyMatrix.set(7, 7, MyMatrix.get(7, 7)+coef*(((N4x*J22)/detJ-(N4y*J12)/detJ)*(e33*((N4x*J22)/detJ-(N4y*J12)/detJ)+e23*((N4y*J11)/detJ-(N4x*J21)/detJ))+((N4y*J11)/detJ-(N4x*J21)/detJ)*(e23*((N4x*J22)/detJ-(N4y*J12)/detJ)+e22*((N4y*J11)/detJ-(N4x*J21)/detJ))));
                
            }
        }
        MyMatrix.set(1,0,MyMatrix.get(0, 1));
        MyMatrix.set(2,0,MyMatrix.get(0, 2)); MyMatrix.set(2,1,MyMatrix.get(1, 2));
        MyMatrix.set(3,0,MyMatrix.get(0, 3)); MyMatrix.set(3,1,MyMatrix.get(1, 3)); MyMatrix.set(3,2,MyMatrix.get(2, 3));
        MyMatrix.set(4,0,MyMatrix.get(0, 4)); MyMatrix.set(4,1,MyMatrix.get(1, 4)); MyMatrix.set(4,2,MyMatrix.get(2, 4)); MyMatrix.set(4,3,MyMatrix.get(3, 4));
        MyMatrix.set(5,0,MyMatrix.get(0, 5)); MyMatrix.set(5,1,MyMatrix.get(1, 5)); MyMatrix.set(5,2,MyMatrix.get(2, 5)); MyMatrix.set(5,3,MyMatrix.get(3, 5)); MyMatrix.set(5,4,MyMatrix.get(4, 5));
        MyMatrix.set(6,0,MyMatrix.get(0, 6)); MyMatrix.set(6,1,MyMatrix.get(1, 6)); MyMatrix.set(6,2,MyMatrix.get(2, 6)); MyMatrix.set(6,3,MyMatrix.get(3, 6)); MyMatrix.set(6,4,MyMatrix.get(4, 6)); MyMatrix.set(6,5,MyMatrix.get(5, 6));
        MyMatrix.set(7,0,MyMatrix.get(0, 7)); MyMatrix.set(7,1,MyMatrix.get(1, 7)); MyMatrix.set(7,2,MyMatrix.get(2, 7)); MyMatrix.set(7,3,MyMatrix.get(3, 7)); MyMatrix.set(7,4,MyMatrix.get(4, 7)); MyMatrix.set(7,5,MyMatrix.get(5, 7)); MyMatrix.set(7,6,MyMatrix.get(6, 7));
        return MyMatrix;
    }
    
    abstract AbstractMatrix getEmat();

    @Override
    public double[] getBVNorm(double coef) {
        // to be fixed ....
        double[] norm= new double[2];
        norm[0]=0.; norm[1]=0.;
        double abs_u1=Math.sqrt(this.getNodeHierarchy(1).getDisp(1)*this.getNodeHierarchy(1).getDisp(1)+this.getNodeHierarchy(1).getDisp(2)*this.getNodeHierarchy(1).getDisp(2))/coef;
        double abs_u2=Math.sqrt(this.getNodeHierarchy(2).getDisp(1)*this.getNodeHierarchy(2).getDisp(1)+this.getNodeHierarchy(2).getDisp(2)*this.getNodeHierarchy(2).getDisp(2))/coef;
        double abs_u3=Math.sqrt(this.getNodeHierarchy(3).getDisp(1)*this.getNodeHierarchy(3).getDisp(1)+this.getNodeHierarchy(3).getDisp(2)*this.getNodeHierarchy(3).getDisp(2))/coef;
        double abs_u4=Math.sqrt(this.getNodeHierarchy(4).getDisp(1)*this.getNodeHierarchy(4).getDisp(1)+this.getNodeHierarchy(4).getDisp(2)*this.getNodeHierarchy(4).getDisp(2))/coef;

        double u1_x,u1_y,u2_x,u2_y,u3_x,u3_y,u4_x,u4_y;
        double xsi=-1.; double ysi =-1.;
        u1_x=ShapeFunction_xsi(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_xsi(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_xsi(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_xsi(4,xsi,ysi)*abs_u4/coef;
        u1_y=ShapeFunction_eta(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_eta(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_eta(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_eta(4,xsi,ysi)*abs_u4/coef;
        
        xsi=1.; ysi =-1.;
        u2_x=ShapeFunction_xsi(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_xsi(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_xsi(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_xsi(4,xsi,ysi)*abs_u4/coef;
        u2_y=ShapeFunction_eta(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_eta(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_eta(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_eta(4,xsi,ysi)*abs_u4/coef;
        
        xsi=1.; ysi =1.;
        u3_x=ShapeFunction_xsi(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_xsi(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_xsi(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_xsi(4,xsi,ysi)*abs_u4/coef;
        u3_y=ShapeFunction_eta(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_eta(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_eta(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_eta(4,xsi,ysi)*abs_u4/coef;
        
        xsi=-1.; ysi =1.;
        u4_x=ShapeFunction_xsi(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_xsi(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_xsi(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_xsi(4,xsi,ysi)*abs_u4/coef;
        u4_y=ShapeFunction_eta(1,xsi,ysi)*abs_u1/coef
                +ShapeFunction_eta(2,xsi,ysi)*abs_u2/coef
                +ShapeFunction_eta(3,xsi,ysi)*abs_u3/coef
                +ShapeFunction_eta(4,xsi,ysi)*abs_u4/coef;
        double Du1=Math.sqrt(u1_x*u1_x+u1_y*u1_y);
        double Du2=Math.sqrt(u2_x*u2_x+u2_y*u2_y);
        double Du3=Math.sqrt(u3_x*u3_x+u3_y*u3_y);
        double Du4=Math.sqrt(u4_x*u4_x+u4_y*u4_y);
        for(int i=1; i<=GaussPoints; i++){
            for(int j=1; j<=GaussPoints; j++){
                double x1=this.getGaussCoord(i);
                double x2=this.getGaussCoord(j);
                double w1=this.getGaussWeight(i);
                double w2=this.getGaussWeight(j);
                
                double N1=ShapeFunction(1,x1,x2);
                double N2=ShapeFunction(2,x1,x2);
                double N3=ShapeFunction(3,x1,x2);
                double N4=ShapeFunction(4,x1,x2);
                
                double detJ=getDetJ(x1, x2);
                
                norm[0]+=(N1*abs_u1+ N2*abs_u2+ N3*abs_u3+ N4*abs_u4)*detJ*w1*w2;
                norm[1]+=(N1*abs_u1+ N2*abs_u2+ N3*abs_u3+ N4*abs_u4+N1*Du1+ N2*Du2+ N3*Du3+ N4*Du4)*detJ*w1*w2;
            }
        }
        return norm;
    }
    
    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
     public double[] getStress(double[] at, int LC,int step){
         double xsi=at[0];
         double eta=at[1];
         AbstractMatrix d = new AbstractMatrix(this.ndofs,1,0.0);
         int count=0;
         for(int i=1;i<=this.numNodes;i++){
             for(int j=0;j<this.dof_per_node;j++){
                 d.set(count, 0, getNodeHierarchy(i).getLoadCaseDisps(LC,step)[j]);
                 count++;
             }
         }
         AbstractMatrix M = getEmat().times(getB(xsi, eta).times(d));
         double[] stress = new double[3];
         stress[0]=M.get(0, 0);
         stress[1]=M.get(1, 0);
         stress[2]=M.get(2, 0);
         return stress;
     }
}
