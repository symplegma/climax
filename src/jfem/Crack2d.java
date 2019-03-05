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
public class Crack2d extends Element {
    private static int numberOfCrack2d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    private boolean sens=false;
    private double ky2kx=1.;
    // some variables needed for elastoplastic analysis
    MaterialElasticPlasticPoint theMaterialPoint;
    
    // constructors
    public Crack2d(){
        
    }
    
    public Crack2d(int id, Node Node1, Node Node2, Node Node3, Node Node4, Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfCrack2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        this.putNode(Node3);
        this.putNode(Node4);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        Node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, Node1.getID());
        Node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, Node2.getID());
        Node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, Node3.getID());
        Node4.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(4, Node4.getID());

        double[] coords = new double[2];
        coords[0]=(Node1.getCoords()[0]+Node2.getCoords()[0])/2.;
        coords[1]=(Node1.getCoords()[1]+Node2.getCoords()[1])/2.;
        theMaterialPoint = new MaterialElasticPlasticPoint();
        theMaterialPoint.setMaterial(Mater);
        theMaterialPoint.initMatPoint(2);
        theMaterialPoint.setCoords(coords,1.0);
        ndofs = 8;
    }
    
    // methods
    public void set_ky2kx(double val){this.ky2kx=val;}
    
    public int getNumberOfCrack2d() {
        return numberOfCrack2d;
    }
    
    public int getElementNdofs() {
        return ndofs;
    }
    
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
    
    public AbstractMatrix getK() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double L=this.getL();
        double E;
        
        double z=1.;
        if(this.damage)z=this.getNodeHierarchy(1).getDamageVariableCur()*z;
        if(this.ElementMaterial.getType()==0){
            E=((ElasticMaterial) ElementMaterial).getElasticity()*z;
        } else {
            ElastoPlastic aMaterial =(ElastoPlastic) this.ElementMaterial;
            if(this.ElementPlastified()){
                E=( aMaterial.getElasticity()*z*(aMaterial.getKisotropic() +aMaterial.getKkinematic()) )/
                  ( aMaterial.getElasticity()*z+aMaterial.getKisotropic()+aMaterial.getKkinematic() );
            } else {
                E=aMaterial.getElasticity()*z;
            }
        }
        double b=this.theCrossSection.getThickness();
        double kx=E*b*L;
        double ky=kx*ky2kx;
        
        MyMatrix.putVal(0, 0,kx/3.);
        MyMatrix.putVal(1, 0,0.);
        MyMatrix.putVal(2, 0,kx/6.);
        MyMatrix.putVal(3, 0,0.);
        MyMatrix.putVal(4, 0,-kx/3.);
        MyMatrix.putVal(5, 0,0.);
        MyMatrix.putVal(6, 0,-kx/6.);
        MyMatrix.putVal(7, 0,0.);
        
        MyMatrix.putVal(0, 1,0.);
        MyMatrix.putVal(1, 1,ky/3.);
        MyMatrix.putVal(2, 1,0.);
        MyMatrix.putVal(3, 1,ky/6.);
        MyMatrix.putVal(4, 1,0.);
        MyMatrix.putVal(5, 1,-ky/6.);
        MyMatrix.putVal(6, 1,0.);
        MyMatrix.putVal(7, 1,-ky/3.);
        
        MyMatrix.putVal(0, 2,kx/6.);
        MyMatrix.putVal(1, 2,0.);
        MyMatrix.putVal(2, 2,kx/3.);
        MyMatrix.putVal(3, 2,0.);
        MyMatrix.putVal(4, 2,-kx/3.);
        MyMatrix.putVal(5, 2,0.);
        MyMatrix.putVal(6, 2,-kx/6.);
        MyMatrix.putVal(7, 2,0.);
        
        MyMatrix.putVal(0, 3,0.);
        MyMatrix.putVal(1, 3,ky/6.);
        MyMatrix.putVal(2, 3,0.);
        MyMatrix.putVal(3, 3,ky/3.);
        MyMatrix.putVal(4, 3,0.);
        MyMatrix.putVal(5, 3,-ky/3.);
        MyMatrix.putVal(6, 3,0.);
        MyMatrix.putVal(7, 3,-ky/6.);
        
        MyMatrix.putVal(0, 4,-kx/6.);
        MyMatrix.putVal(1, 4,0.);
        MyMatrix.putVal(2, 4,-kx/3.);
        MyMatrix.putVal(3, 4,0.);
        MyMatrix.putVal(4, 4,kx/3.);
        MyMatrix.putVal(5, 4,0.);
        MyMatrix.putVal(6, 4,kx/6.);
        MyMatrix.putVal(7, 4,0.);
        
        MyMatrix.putVal(0, 5,0.);
        MyMatrix.putVal(1, 5,-ky/6.);
        MyMatrix.putVal(2, 5,0.);
        MyMatrix.putVal(3, 5,-ky/3.);
        MyMatrix.putVal(4, 5,0.);
        MyMatrix.putVal(5, 5,ky/3.);
        MyMatrix.putVal(6, 5,0.);
        MyMatrix.putVal(7, 5,ky/6.);
        
        MyMatrix.putVal(0, 6,-kx/3.);
        MyMatrix.putVal(1, 6,0.);
        MyMatrix.putVal(2, 6,--kx/6.);
        MyMatrix.putVal(3, 6,0.);
        MyMatrix.putVal(4, 6,kx/6.);
        MyMatrix.putVal(5, 6,0.);
        MyMatrix.putVal(6, 6,kx/3.);
        MyMatrix.putVal(7, 6,0.);
        
        MyMatrix.putVal(0, 7,0.);
        MyMatrix.putVal(1, 7,-ky/3.);
        MyMatrix.putVal(2, 7,0.);
        MyMatrix.putVal(3, 7,-ky/6.);
        MyMatrix.putVal(4, 7,0.);
        MyMatrix.putVal(5, 7,ky/6.);
        MyMatrix.putVal(6, 7,0.);
        MyMatrix.putVal(7, 7,ky/3.);
        
        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
//        return MyMatrix;
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
        TrMatrix.putVal(6, 0,0.);
        TrMatrix.putVal(7, 0,0.);
        
        TrMatrix.putVal(0, 1,s);
        TrMatrix.putVal(1, 1,c);
        TrMatrix.putVal(2, 1,0.);
        TrMatrix.putVal(3, 1,0.);
        TrMatrix.putVal(4, 1,0.);
        TrMatrix.putVal(5, 1,0.);
        TrMatrix.putVal(6, 1,0.);
        TrMatrix.putVal(7, 1,0.);
        
        TrMatrix.putVal(0, 2,0.);
        TrMatrix.putVal(1, 2,0.);
        TrMatrix.putVal(2, 2,c);
        TrMatrix.putVal(3, 2,-s);
        TrMatrix.putVal(4, 2,0.);
        TrMatrix.putVal(5, 2,0.);
        TrMatrix.putVal(6, 2,0.);
        TrMatrix.putVal(7, 2,0.);
        
        TrMatrix.putVal(0, 3,0.);
        TrMatrix.putVal(1, 3,0.);
        TrMatrix.putVal(2, 3,s);
        TrMatrix.putVal(3, 3,c);
        TrMatrix.putVal(4, 3,0.);
        TrMatrix.putVal(5, 3,0.);
        TrMatrix.putVal(6, 3,0.);
        TrMatrix.putVal(7, 3,0.);
        
        
        TrMatrix.putVal(0, 4,0.);
        TrMatrix.putVal(1, 4,0.);
        TrMatrix.putVal(2, 4,0.);
        TrMatrix.putVal(3, 4,0.);
        TrMatrix.putVal(4, 4,c);
        TrMatrix.putVal(5, 4,-s);
        TrMatrix.putVal(6, 4,0.);
        TrMatrix.putVal(7, 4,0.);
        
        TrMatrix.putVal(0, 5,0.);
        TrMatrix.putVal(1, 5,0.);
        TrMatrix.putVal(2, 5,0.);
        TrMatrix.putVal(3, 5,0.);
        TrMatrix.putVal(4, 5,s);
        TrMatrix.putVal(5, 5,c);
        TrMatrix.putVal(6, 5,0.);
        TrMatrix.putVal(7, 5,0.);
        
        TrMatrix.putVal(0, 6,0.);
        TrMatrix.putVal(1, 6,0.);
        TrMatrix.putVal(2, 6,0.);
        TrMatrix.putVal(3, 6,0.);
        TrMatrix.putVal(4, 6,0.);
        TrMatrix.putVal(5, 6,0.);
        TrMatrix.putVal(6, 6,c);
        TrMatrix.putVal(7, 6,-s);
        
        TrMatrix.putVal(0, 7,0.);
        TrMatrix.putVal(1, 7,0.);
        TrMatrix.putVal(2, 7,0.);
        TrMatrix.putVal(3, 7,0.);
        TrMatrix.putVal(4, 7,0.);
        TrMatrix.putVal(5, 7,0.);
        TrMatrix.putVal(6, 7,s);
        TrMatrix.putVal(7, 7,c);

        return TrMatrix;
    }
    
    private AbstractMatrix getM(int consistent) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        MyMatrix.init();
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
            i+=2;
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
            i+=2;
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
            i+=2;
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
            i+=2;
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
            i+=2;
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
            j=j+2;
        }
        return theEFTable;
    }
//
    public AbstractMatrix getF() {
        double[] vec = new double[2];
        double[] vec2 = new double[2];
        double z=1.;
        if(this.damage)z=this.getNodeHierarchy(1).getDamageVariableCur();
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

        //ex
        vec[0]=(disps.get(6, 0)-disps.get(0, 0)+disps.get(4, 0)-disps.get(2, 0))/2.0;
        vec[1]=(disps.get(7, 0)-disps.get(1, 0)+disps.get(5, 0)-disps.get(3, 0))/2.0;
        vec2=this.theMaterialPoint.giveStrainTakeStress(vec);
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        double a=this.theCrossSection.getThickness()*this.getL()/2.;
//        vec[0]=disps.get(6, 0)-disps.get(0, 0);
//        vec[1]=disps.get(7, 0)-disps.get(1, 0);
//        vec2=this.theMaterialPoint.giveStrainTakeStress(vec);
        Fint.addVal(0, 0, -a*vec2[0]);
        Fint.addVal(1, 0, -a*vec2[1]);
        Fint.addVal(6, 0, a*vec2[0]);
        Fint.addVal(7, 0, a*vec2[1]);
//        vec[0]=disps.get(4, 0)-disps.get(2, 0);
//        vec[1]=disps.get(5, 0)-disps.get(3, 0);
//        vec2=this.theMaterialPoint.giveStrainTakeStress(vec);
        Fint.addVal(2, 0, -a*vec2[0]);
        Fint.addVal(3, 0, -a*vec2[1]);
        Fint.addVal(4, 0, a*vec2[0]);
        Fint.addVal(5, 0, a*vec2[1]);
        Fint=this.getT().transpose().times(Fint);
        return Fint.times(z);
//        return this.getK_u();
    }
    
    @Override
    void clear() {
        theMaterialPoint.clear();
    }

    @Override
    void commit() {
        theMaterialPoint.commit();
    }
    
    @Override
    public void commitDamage() {
        double zprv=this.getNodeHierarchy(1).getDamageVariablePrv();
        double zcur=this.getNodeHierarchy(1).getDamageVariableCur();
        if(zprv>=1.0e-12 && zcur<=1.0e-12)this.getNodeHierarchy(1).setDamageVariablePrv(0.0);
        if(zprv<1.0e-12)this.getNodeHierarchy(1).setDamageVariableCur(0.0);
    }
    
    @Override
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
    
    @Override
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

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean ElementPlastified() {
        return this.theMaterialPoint.AskIfIsPlastified();
    }
    
    @Override
    public double DamageDissipation(){
        double diss=((ElasticMaterial) this.ElementMaterial).getFracTough()*this.getL()*(
                -this.getNodeHierarchy(1).getDamageVariableCur()
                +this.getNodeHierarchy(1).getDamageVariablePrv()
                );
        return diss;
    }
    
    @Override
    public boolean checkForDamage(){
        boolean answ=false;
        if(this.damage){
        double nrgvalprv,nrgvalcur;
        double zprv=this.getNodeHierarchy(1).getDamageVariablePrv();
        nrgvalprv = this.getuKu()/2. + this.DamageDissipation();
        nrgvalcur = nrgvalprv;
        if(zprv>=1.0e-12){
            ElementNodes.get(1).setDamageVariableCur(0.0);
            answ=true;
            nrgvalcur = this.getuKu()/2. + this.DamageDissipation();
            if(nrgvalcur>nrgvalprv){
                this.getNodeHierarchy(1).setDamageVariableCur(zprv);
                answ=false;
            }else{
//                System.out.println("element "+id+" damaged");
            }
        }}
        return answ;
    }
    
    public double getVolume(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
