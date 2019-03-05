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
public class Truss2d extends Element {
    private static int numberOfTruss2d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    private boolean sens=false;
    
    // some variables needed for elastoplastic analysis
    MaterialElasticPlasticPoint theMaterialPoint;
    /*private double  Stress=0.;
    private double  Eplastic=0.;
    private double  Etotal=0.;
    private double  Xsi=0.;
    private double  YieldFunction=0.;
    private double  aPar=0.;
    private double  qPar=0.;
    
    private double  EplasticTrial=0.;
    private double  aParTrial=0.;
    private double  qParTrial=0.;*/

    // constructors
    public Truss2d(){
        dimension=1;
    }
    
    public Truss2d(int id, Node Node1, Node Node2, Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfTruss2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        Node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, Node1.getID());
        Node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, Node2.getID());

        double[] coords = new double[2];
        coords[0]=Node1.getCoords()[0]+Node2.getCoords()[0];
        coords[1]=Node1.getCoords()[1]+Node2.getCoords()[1];
        theMaterialPoint = new MaterialElasticPlasticPoint();
        theMaterialPoint.setMaterial(Mater);
        theMaterialPoint.initMatPoint(1);
        theMaterialPoint.setCoords(coords,1.0);
        this.dof_per_node=2;
        ndofs = 4;
        dimension=1;
    }
    
    // methods
    public int getNumberOfTruss2d() {
        return numberOfTruss2d;
    }
    
    public int getElementNdofs() {
        return ndofs;
    }
    
    public AbstractMatrix getK() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double L=this.getL();
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
        double A=this.theCrossSection.getA();
        double k=E*A/L;
        if(this.sens==true){
            k=E/L;
        }
        
        MyMatrix.putVal(0, 0,k);
        MyMatrix.putVal(1, 0,0.);
        MyMatrix.putVal(2, 0,-k);
        MyMatrix.putVal(3, 0,0.);
        
        MyMatrix.putVal(0, 1,0.);
        MyMatrix.putVal(1, 1,0.);
        MyMatrix.putVal(2, 1,0.);
        MyMatrix.putVal(3, 1,0.);
        
        MyMatrix.putVal(0, 2,-k);
        MyMatrix.putVal(1, 2,0.);
        MyMatrix.putVal(2, 2,k);
        MyMatrix.putVal(3, 2,0.);
        
        MyMatrix.putVal(0, 3,0.);
        MyMatrix.putVal(1, 3,0.);
        MyMatrix.putVal(2, 3,0.);
        MyMatrix.putVal(3, 3,0.);

        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }
    
    private AbstractMatrix getM(int consistent) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double L=this.getL();
        double p=((ElasticMaterial) this.ElementMaterial).getDensity();
        double A=this.theCrossSection.getA();
        double m=p*A*L;
        switch(consistent){
            case 0:
                m=m/6;
                MyMatrix.putVal(0, 0,2.*m);
                MyMatrix.putVal(1, 0,0.);
                MyMatrix.putVal(2, 0,m);
                MyMatrix.putVal(3, 0,0.);
                
                MyMatrix.putVal(0, 1,0.);
                MyMatrix.putVal(1, 1,2.*m);
                MyMatrix.putVal(2, 1,0.);
                MyMatrix.putVal(3, 1,m);
                
                MyMatrix.putVal(0, 2,m);
                MyMatrix.putVal(1, 2,0.);
                MyMatrix.putVal(2, 2,2.*m);
                MyMatrix.putVal(3, 2,0.);
                
                MyMatrix.putVal(0, 3,0.);
                MyMatrix.putVal(1, 3,m);
                MyMatrix.putVal(2, 3,0.);
                MyMatrix.putVal(3, 3,2.*m);
                break;
            default:
                m=m/2;
                MyMatrix.putVal(0, 0,m);
                MyMatrix.putVal(1, 0,0.);
                MyMatrix.putVal(2, 0,0.);
                MyMatrix.putVal(3, 0,0.);
                
                MyMatrix.putVal(0, 1,0.);
                MyMatrix.putVal(1, 1,m);
                MyMatrix.putVal(2, 1,0.);
                MyMatrix.putVal(3, 1,0.);
                
                MyMatrix.putVal(0, 2,0.);
                MyMatrix.putVal(1, 2,0.);
                MyMatrix.putVal(2, 2,m);
                MyMatrix.putVal(3, 2,0.);
                
                MyMatrix.putVal(0, 3,0.);
                MyMatrix.putVal(1, 3,0.);
                MyMatrix.putVal(2, 3,0.);
                MyMatrix.putVal(3, 3,m);
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
        
        TrMatrix.putVal(0, 1,s);
        TrMatrix.putVal(1, 1,c);
        TrMatrix.putVal(2, 1,0.);
        TrMatrix.putVal(3, 1,0.);
        
        TrMatrix.putVal(0, 2,0.);
        TrMatrix.putVal(1, 2,0.);
        TrMatrix.putVal(2, 2,c);
        TrMatrix.putVal(3, 2,-s);
        
        TrMatrix.putVal(0, 3,0.);
        TrMatrix.putVal(1, 3,0.);
        TrMatrix.putVal(2, 3,s);
        TrMatrix.putVal(3, 3,c);
        

        return TrMatrix;
    }

    public AbstractMatrix getF() {
        double[] vec = new double[1];
        double[] vec2 = new double[1];
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
        vec[0]=(-disps.get(0, 0)+disps.get(2, 0))/this.getL();
        vec2[0]=this.theMaterialPoint.giveStrainTakeStress(vec)[0];
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint.addVal(0, 0, -vec2[0]*this.theCrossSection.getA());
        Fint.addVal(2, 0, vec2[0]*this.theCrossSection.getA());
        Fint=this.getT().transpose().times(Fint);
        return Fint;
        //return this.getK_u();
    }

    @Override
    void clear() {
        theMaterialPoint.clear();
    }

    @Override
    void commit() {
        theMaterialPoint.commit();
    }
    
    public void setSens(boolean flag){
        this.sens=flag;
    }

    @Override
    public boolean ElementPlastified() {
        return this.theMaterialPoint.AskIfIsPlastified();
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
        double[] norm= new double[2];
        norm[0]=0.; norm[1]=0.;
        double abs_u1=Math.sqrt(this.getNodeHierarchy(1).getDisp(1)*this.getNodeHierarchy(1).getDisp(1)+this.getNodeHierarchy(1).getDisp(2)*this.getNodeHierarchy(1).getDisp(2))/coef;
        double abs_u2=Math.sqrt(this.getNodeHierarchy(2).getDisp(1)*this.getNodeHierarchy(2).getDisp(1)+this.getNodeHierarchy(2).getDisp(2)*this.getNodeHierarchy(2).getDisp(2))/coef;
        norm[0]=0.5*(abs_u1+abs_u2)*this.getL();
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
    public AbstractMatrix getF(int LC, int step){
        double[] vec = new double[1];
        double[] vec2 = new double[1];
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getLoadCaseDisps(LC,step)[0]);
            disps.addVal(i+1, 0, theNode.getLoadCaseDisps(LC,step)[1]);
            i+=2;
        }
        // from global to local disps
        disps=this.getT().times(disps);

        //ex
        vec[0]=(-disps.get(0, 0)+disps.get(2, 0))/this.getL();
        vec2[0]=this.theMaterialPoint.giveStrainTakeStress(vec)[0];
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1);
        Fint.addVal(0, 0, -vec2[0]*this.theCrossSection.getA());
        Fint.addVal(2, 0, vec2[0]*this.theCrossSection.getA());
        //Fint=this.getT().transpose().times(Fint);
        return Fint;
    }
    
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
                Fequivalent.set(2, 0, val*L/2.0);
                break;
            case 2:
                Fequivalent.set(1, 0, val*L/2.0);
                Fequivalent.set(3, 0, val*L/2.0);
                break;
            default:
                System.err.println("Not defined equivalent load for element "+this.getClass().getName()+", in respect to dof "+wd);
        }
        return Fequivalent;
    }
    
    @Override
     public double[] getStress(double[] at, int LC,int step){
         AbstractMatrix d = new AbstractMatrix(this.ndofs,1,0.0);
         int count=0;
         for(int i=1;i<=this.numNodes;i++){
             for(int j=0;j<this.dof_per_node;j++){
                 d.set(count, 0, getNodeHierarchy(i).getLoadCaseDisps(LC,step)[j]);
                 count++;
             }
         }
         AbstractMatrix M = this.getT().times(d);
         double[] stress = new double[1];
         stress[0]=(M.get(2, 0)-M.get(0, 0))*((ElasticMaterial)this.ElementMaterial).getElasticity()/this.getL();
         return stress;
     }
    
    @Override
    public double getStressVonMisses(double[] at, int LC,int step){
        double[] Stress = this.getStress(at, LC, step);
        return Stress[0];
    }

}