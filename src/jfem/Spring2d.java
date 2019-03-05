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
public class Spring2d extends Element {
    private static int numberOfSpring2d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    private boolean sens=false;
    // some variables needed for elastoplastic analysis
    MaterialElasticPlasticPoint theMaterialPoint;
    private double c,s;
    
    // constructors
    public Spring2d(){
        
    }
    
    public Spring2d(int id, Node Node1, Node Node2, Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfSpring2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        int[] dofs = new int[3];
        dofs[0]=1; dofs[1]=2; dofs[2]=6;
        Node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, Node1.getID());
        Node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, Node2.getID());

        double[] coords = new double[2];
        coords[0]=Node1.getCoords()[0]+Node2.getCoords()[0];
        coords[1]=Node1.getCoords()[1]+Node2.getCoords()[1];
        theMaterialPoint = new MaterialElasticPlasticPoint();
        theMaterialPoint.setMaterial(Mater);
        theMaterialPoint.initMatPoint(1);
        theMaterialPoint.setCoords(coords,1.0);
        c=0.0;
        s=1.0;
        double L=getL();
        if(L>0.00001){
            c=(Node2.getCoords()[0]-Node1.getCoords()[0])/L;
            s=(Node2.getCoords()[1]-Node1.getCoords()[1])/L;
        }
        ndofs = 6;
        
    }
    
    public Spring2d(int id, Node Node1, Node Node2, Material Mater, CrossSection Sect, double theta){
        this.id=id;
        ++numberOfSpring2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        int[] dofs = new int[3];
        dofs[0]=1; dofs[1]=2; dofs[2]=6;
        Node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, Node1.getID());
        Node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, Node2.getID());

        double[] coords = new double[2];
        coords[0]=Node1.getCoords()[0]+Node2.getCoords()[0];
        coords[1]=Node1.getCoords()[1]+Node2.getCoords()[1];
        theMaterialPoint = new MaterialElasticPlasticPoint();
        theMaterialPoint.setMaterial(Mater);
        theMaterialPoint.initMatPoint(1);
        theMaterialPoint.setCoords(coords,1.0);
        c=Math.cos(theta);
        s=Math.sin(theta);
        ndofs = 6;
        
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
    
    // methods
    public void setOrientation(double cval, double sval){
        this.c=cval;
        this.s=sval;
    }
    public void setOrientation(double theta){
        this.c=Math.cos(theta);
        this.s=Math.sin(theta);
    }
    
    public int getNumberOfSpring2d() {
        return numberOfSpring2d;
    }
    
    public int getElementNdofs() {
        return ndofs;
    }
    
    public AbstractMatrix getK() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double E;
        
        if(this.ElementMaterial.getType()==0){
            E=((ElasticMaterial) ElementMaterial).getElasticity();
        } else {
            ElastoPlastic aMaterial =(ElastoPlastic) this.ElementMaterial;
            if(this.ElementPlastified()){
                E=( aMaterial.getElasticity()*(aMaterial.getKisotropic() +aMaterial.getKkinematic()) )/
                  ( aMaterial.getElasticity()+aMaterial.getKisotropic()+aMaterial.getKkinematic() );
            } else {
                E=aMaterial.getElasticity();
            }
        }
        double Ax=this.theCrossSection.getA();
        double Ay=this.theCrossSection.getAs2();
        double Az=this.theCrossSection.getAs3();
        double kx=E*Ax;
        double ky=E*Ay;
        double kz=E*Az;
        
        MyMatrix.putVal(0, 0,kx);
        MyMatrix.putVal(1, 0,0.);
        MyMatrix.putVal(2, 0,0.);
        MyMatrix.putVal(3, 0,-kx);
        MyMatrix.putVal(4, 0,0.);
        MyMatrix.putVal(5, 0,0.);
        
        MyMatrix.putVal(0, 1,0.);
        MyMatrix.putVal(1, 1,ky);
        MyMatrix.putVal(2, 1,0.);
        MyMatrix.putVal(3, 1,0.);
        MyMatrix.putVal(4, 1,-ky);
        MyMatrix.putVal(5, 1,0.);
        
        MyMatrix.putVal(0, 2,0.);
        MyMatrix.putVal(1, 2,0.);
        MyMatrix.putVal(2, 2,kz);
        MyMatrix.putVal(3, 2,0.);
        MyMatrix.putVal(4, 2,0.);
        MyMatrix.putVal(5, 2,-kz);
        
        MyMatrix.putVal(0, 3,-kx);
        MyMatrix.putVal(1, 3,0.);
        MyMatrix.putVal(2, 3,0.);
        MyMatrix.putVal(3, 3,kx);
        MyMatrix.putVal(4, 3,0.);
        MyMatrix.putVal(5, 3,0.);
        
        MyMatrix.putVal(0, 4,0.);
        MyMatrix.putVal(1, 4,-ky);
        MyMatrix.putVal(2, 4,0.);
        MyMatrix.putVal(3, 4,0.);
        MyMatrix.putVal(4, 4,ky);
        MyMatrix.putVal(5, 4,0.);
        
        MyMatrix.putVal(0, 5,0.);
        MyMatrix.putVal(1, 5,0.);
        MyMatrix.putVal(2, 5,-kz);
        MyMatrix.putVal(3, 5,0.);
        MyMatrix.putVal(4, 5,0.);
        MyMatrix.putVal(5, 5,kz);

        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    public AbstractMatrix getT() {
        TrMatrix = new AbstractMatrix(ndofs,ndofs);
        
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
    
    public AbstractMatrix getF() {
//        double[] vec = new double[1];
//        double[] vec2 = new double[1];
//        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
//        int i=0;
//        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
//            Node theNode = it.next();
//            disps.addVal(i, 0, theNode.getDispsTrial()[0]);
//            disps.addVal(i+1, 0, theNode.getDispsTrial()[1]);
//            disps.addVal(i+2, 0, theNode.getDispsTrial()[5]);
//            i+=3;
//        }
//        // from global to local disps
//        disps=this.getT().times(disps);
//
//        //ex
//        vec[0]=(-disps.get(0, 0)+disps.get(2, 0));
//        vec2[0]=this.theMaterialPoint.giveStrainTakeStress(vec)[0];
//        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
//        Fint.addVal(0, 0, -vec2[0]*this.theCrossSection.getA());
//        Fint.addVal(2, 0, vec2[0]*this.theCrossSection.getA());
//        Fint=this.getT().transpose().times(Fint);
//        return Fint;
        return this.getK_u();
    }

    @Override
    public boolean ElementPlastified() {
        return this.theMaterialPoint.AskIfIsPlastified();
    }

    @Override
    public AbstractMatrix getM() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public int[] getFtable(){
        int[] theEFTable = new int[ndofs];
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theEFTable[j]=theNode.getFtable()[0];
            theEFTable[j+1]=theNode.getFtable()[1];
            theEFTable[j+2]=theNode.getFtable()[2];
            j=j+3;
        }
        return theEFTable;
    }

    @Override
    void clear() {
        theMaterialPoint.clear();
    }

    @Override
    AbstractMatrix getM_v() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getM_a() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getM_u() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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

    @Override
    void commit() {
        theMaterialPoint.commit();
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getK_v() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getK_a() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getVolume(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
