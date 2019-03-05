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
import java.util.TreeMap;
import java.util.Map;

/**
 *
 * @author pchr
 */
public class TNTruss2d extends Element {
    private static int numberOfTNTruss2d = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    private boolean sens=false;
    private static AbstractMatrix Matrix_1_6 = new AbstractMatrix(1,6);
    private static AbstractMatrix Matrix_2_6 = new AbstractMatrix(2,6);
    // some variables needed for elastoplastic analysis
    Map<Integer,MaterialElasticPlasticPoint> theMaterialPoints = new TreeMap<Integer,MaterialElasticPlasticPoint>();
    private static MaterialElasticPlasticPoint aMaterialPoint;
    // this element gives a zero Jacobian when the midside node is on 1/3 or 2/3 of the total length

    // constructors
    public TNTruss2d(){

    }

    public TNTruss2d(int id, Node Node1, Node Node2, Node Node3, Material Mater, CrossSection Sect, int numMP){
        this.id=id;
        ++numberOfTNTruss2d;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        this.putNode(Node3);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        Node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, Node1.getID());
        Node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, Node2.getID());
        Node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, Node3.getID());
        this.theGaussData = new GaussData(numMP);
        this.numMP=numMP;

        double[] coords = new double[1];
        for(int i=1;i<=numMP;i++){
            aMaterialPoint = new MaterialElasticPlasticPoint();
            coords[0]=theGaussData.getGaussCoordinate(i-1);
            aMaterialPoint.setCoords(coords,theGaussData.getGaussWeight(i-1));
            aMaterialPoint.setMaterial(Mater);
            aMaterialPoint.initMatPoint(1);
            this.theMaterialPoints.put(aMaterialPoint.getID(), aMaterialPoint);
        }
        ndofs = 6;
        dimension=1;
    }

    public TNTruss2d(Node Node1, Node Node2, Node Node3, Material Mater, CrossSection Sect){
        ++numberOfTNTruss2d;
        this.numMP=3;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        this.putNode(Node2);
        this.putNode(Node3);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        Node1.setNdofs_ofNode(dofs);
        Node2.setNdofs_ofNode(dofs);
        Node3.setNdofs_ofNode(dofs);
        this.theGaussData = new GaussData(numMP);

        double[] coords = new double[1];
        for(int i=1;i<=numMP;i++){
            aMaterialPoint = new MaterialElasticPlasticPoint();
            coords[0]=theGaussData.getGaussCoordinate(i-1);
            aMaterialPoint.setCoords(coords,theGaussData.getGaussWeight(i-1));
            aMaterialPoint.setMaterial(Mater);
            aMaterialPoint.initMatPoint(1);
            this.theMaterialPoints.put(aMaterialPoint.getID(), aMaterialPoint);
        }
    }

    // methods
    public int getNumberOfTNTruss2d() {
        return numberOfTNTruss2d;
    }

    public int getElementNdofs() {
        return ndofs;
    }

    public double ShapeFunction(int wsf, double xsi){
        double N;
        switch(wsf){
            case 1: N=0.5*xsi*(xsi-1.); break;
            case 2: N=1.-xsi*xsi; break;
            case 3: N=0.5*xsi*(xsi+1.); break;
            default:N=0. ;  break;
        }
        return N;
    }

    public double ShapeFunction_xsi(int wsf, double xsi){
        double N;
        switch(wsf){
            case 1: N=xsi-0.5; break;
            case 2: N=-2.*xsi; break;
            case 3: N=xsi+0.5; break;
            default:N=0. ;  break;
        }
        return N;
    }

    public double getX_xsi(double xsi){
        double val=0.;
        for(int i=1;i<=3;i++){
            val+=ShapeFunction_xsi(i,xsi)*this.ElementNodes.get(i).getCoords()[0];
        }
        return val;
    }

    public double getY_xsi(double xsi){
        double val=0.;
        for(int i=1;i<=3;i++){
            val+=ShapeFunction_xsi(i,xsi)*this.ElementNodes.get(i).getCoords()[1];
        }
        return val;
    }

    private double getEpsilon(double xsi){
        double val=0.;
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int cc=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(cc, 0, theNode.getDispsTrial()[0]);
            disps.addVal(cc+1, 0, theNode.getDispsTrial()[1]);
            cc+=2;
        }
        // from global to local disps
        disps=this.getT().times(disps);
        val=ShapeFunction_xsi(1,xsi)*disps.get(0, 0)
                +ShapeFunction_xsi(2,xsi)*disps.get(2, 0)
                +ShapeFunction_xsi(3,xsi)*disps.get(4, 0);
        val=val/getDetJ(xsi);
        return val;
    }

    private double getDetJ(double xsi){
        double val=Math.sqrt(Math.pow(getY_xsi(xsi), 2)+Math.pow(getX_xsi(xsi), 2));
        if(val<=0.)System.err.println("warning: element "+this.getID()+" with zero or negative Jacobian at ξ="+xsi);
        return val;
    }

    private AbstractMatrix getB(double xsi){
        Matrix_1_6.init(0.0);

        Matrix_1_6.putVal(0, 0, this.ShapeFunction_xsi(1, xsi));
        Matrix_1_6.putVal(0, 2, this.ShapeFunction_xsi(2, xsi));
        Matrix_1_6.putVal(0, 4, this.ShapeFunction_xsi(3, xsi));

        return Matrix_1_6;
    }

    private AbstractMatrix getH(double xsi){
        Matrix_2_6.init(0.0);

        Matrix_2_6.putVal(0, 0, this.ShapeFunction(1, xsi));
        Matrix_2_6.putVal(0, 2, this.ShapeFunction(2, xsi));
        Matrix_2_6.putVal(0, 4, this.ShapeFunction(3, xsi));

        Matrix_2_6.putVal(1, 1, this.ShapeFunction(1, xsi));
        Matrix_2_6.putVal(1, 3, this.ShapeFunction(2, xsi));
        Matrix_2_6.putVal(1, 5, this.ShapeFunction(3, xsi));

        return Matrix_2_6;
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
        double A=this.theCrossSection.getA();

        MyMatrix.init();
        double xsi,w;

        /*for(int i=1; i<=this.numMP;i++){
            if(this.ElementMaterial.getType()!=0){

            }
            xsi=theGaussData.getGaussCoordinate(i-1);
            MyMatrix=MyMatrix.plus( ( (getB(xsi)).transpose().times(getB(xsi)) ).times(theGaussData.getGaussWeight(i-1)).times(1./this.getDetJ(xsi)).times(E*A) );
        }*/

        for(Iterator<MaterialElasticPlasticPoint> it=this.theMaterialPoints.values().iterator(); it.hasNext();){
            aMaterialPoint = it.next();
            if(this.ElementMaterial.getType()!=0){
                ElastoPlastic aMaterial =(ElastoPlastic) this.ElementMaterial;
                if(aMaterialPoint.AskIfIsPlastified()){
                    E=( aMaterial.getElasticity()*(aMaterial.getKisotropic() +aMaterial.getKkinematic()) )/
                      ( aMaterial.getElasticity()+aMaterial.getKisotropic()+aMaterial.getKkinematic() );
                }else{
                    E=aMaterial.getElasticity();
                }
            }else{
                E=((ElasticMaterial) ElementMaterial).getElasticity();
            }
            xsi = aMaterialPoint.getGaussCoord(0);
            w = aMaterialPoint.getGaussWeight();
            MyMatrix=MyMatrix.plus( ( (getB(xsi)).transpose().times(getB(xsi)) ).times(w).times(1./this.getDetJ(xsi)).times(E*A) );
        }


        MyMatrix=this.getT().transpose().times(MyMatrix);
        return MyMatrix.times(this.getT());
    }

    private AbstractMatrix getM(int consistent) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double p=((ElasticMaterial) this.ElementMaterial).getDensity();
        double A=this.theCrossSection.getA();

        MyMatrix.init(0.0);
        double xsi,w;


        for(Iterator<MaterialElasticPlasticPoint> it=this.theMaterialPoints.values().iterator(); it.hasNext();){
            aMaterialPoint = it.next();

            xsi = aMaterialPoint.getGaussCoord(0);
            w = aMaterialPoint.getGaussWeight();

            MyMatrix=MyMatrix.plus( ( (getH(xsi)).transpose().times(getH(xsi)) ).times(w).times(1./this.getDetJ(xsi)).times(p*A) );
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

    public double getL(){
        double L=0.;
        for(int i=0; i<this.numMP; i++){
            L+=this.getDetJ(this.theGaussData.getGaussCoordinate(i))*this.theGaussData.getGaussWeight(i);
        }
        return L;
    }

    public AbstractMatrix getT() {
        TrMatrix = new AbstractMatrix(ndofs,ndofs);
        TrMatrix.init();

        double L=this.getL();

        double c;
        double s;
        double alpha;
        double xsi;

        xsi=-1.;
        alpha=this.getDetJ(xsi);
        c=getX_xsi(xsi)/alpha;
        s=getY_xsi(xsi)/alpha;
        TrMatrix.putVal(0, 0,c);
        TrMatrix.putVal(1, 0,-s);
        TrMatrix.putVal(0, 1,s);
        TrMatrix.putVal(1, 1,c);
        
        xsi=0.;
        alpha=this.getDetJ(xsi);
        c=getX_xsi(xsi)/alpha;
        s=getY_xsi(xsi)/alpha;
        TrMatrix.putVal(2, 2,c);
        TrMatrix.putVal(3, 2,-s);
        TrMatrix.putVal(2, 3,s);
        TrMatrix.putVal(3, 3,c);
        
        xsi=1.;
        alpha=this.getDetJ(xsi);
        c=getX_xsi(xsi)/alpha;
        s=getY_xsi(xsi)/alpha;
        TrMatrix.putVal(4, 4,c);
        TrMatrix.putVal(5, 4,-s);
        TrMatrix.putVal(4, 5,s);
        TrMatrix.putVal(5, 5,c);

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
        double xsi,w;

        double A =this.theCrossSection.getA();
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1,0.0);

        for(Iterator<MaterialElasticPlasticPoint> it=this.theMaterialPoints.values().iterator(); it.hasNext();){
            aMaterialPoint = it.next();
            xsi = aMaterialPoint.getGaussCoord(0);
            vec[0]=getEpsilon(xsi);
            vec2[0]=aMaterialPoint.giveStrainTakeStress(vec)[0];
            
            w = aMaterialPoint.getGaussWeight();
            Fint=Fint.plus( ( (getB(xsi)).transpose() ).times(w*A*vec2[0]));
        }
        Fint=this.getT().transpose().times(Fint);
        return Fint;
    }

    @Override
    void clear() {
        for(Iterator<MaterialElasticPlasticPoint> it=this.theMaterialPoints.values().iterator(); it.hasNext();){
            it.next().clear();
        }
    }

    @Override
    void commit() {
        for(Iterator<MaterialElasticPlasticPoint> it=this.theMaterialPoints.values().iterator(); it.hasNext();){
            it.next().commit();
        }
    }

    public void setSens(boolean flag){
        this.sens=flag;
    }

    @Override
    public boolean ElementPlastified() {
        boolean check=false;
        for(Iterator<MaterialElasticPlasticPoint> it=this.theMaterialPoints.values().iterator(); it.hasNext();){
            if(it.next().AskIfIsPlastified())check=true;
        }
        return check;
    }

    @Override
    public double getuKu() {
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
            i+=2;
        }
        // from global to local disps
        velcs=this.getT().times(velcs);
        return velcs.transpose().times(getK().times(velcs)).get(0, 0);
    }

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

}