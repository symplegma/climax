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
public class OneDOF extends Element {
    private static int numberOfOneDOF = 0;
    private static AbstractMatrix MyMatrix;
    private static AbstractMatrix TrMatrix;
    private boolean sens=false;
    private int whichDOF;
    private double disp_pr=0.,force_pr=0.;
    private double nrg=-1.;
    
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

    boolean plast;
    // constructors
    public OneDOF(){
        
    }
    
    public OneDOF(int id, Node Node1, Material Mater, CrossSection Sect, int whichDOF){
        this.id=id;
        ++numberOfOneDOF;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(Node1);
        int[] dofs = new int[1];
        dofs[0]=whichDOF;
        Node1.setNdofs_ofNode(dofs); 
        this.theNodesHierarchy.put(1, Node1.getID());
        this.whichDOF=whichDOF;

        theMaterialPoint = new MaterialElasticPlasticPoint();
        theMaterialPoint.setMaterial(Mater);
        theMaterialPoint.initMatPoint(1);
        theMaterialPoint.setCoords(Node1.getCoords(),1.0);
        this.dof_per_node=1;
        ndofs = 1;
    }
    
    // methods
    
    public int getNumberOfOneDOF() {
        return numberOfOneDOF;
    }
    
    public int getElementNdofs() {
        return ndofs;
    }
    
    public AbstractMatrix getK() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);

        double E;
        double z=1.;
        if(this.damage)z=this.getNodeHierarchy(1).getDamageVariableCur();
        if(this.ElementMaterial.getType()==0){
            E=((ElasticMaterial) this.ElementMaterial).getElasticity()*z;
        } else {
            if(plast){
                E=( ((ElastoPlastic) this.ElementMaterial).getElasticity()*z*(((ElastoPlastic) this.ElementMaterial).getKisotropic() +((ElastoPlastic) this.ElementMaterial).getKkinematic()) )/
                  ( ((ElastoPlastic) this.ElementMaterial).getElasticity()*z+((ElastoPlastic) this.ElementMaterial).getKisotropic()+((ElastoPlastic) this.ElementMaterial).getKkinematic() );
            } else {
                E=((ElastoPlastic) this.ElementMaterial).getElasticity()*z;
            }
        }
        double A=this.theCrossSection.getA();
        double k=E*A;
        if(this.sens==true){
            k=E;
        }
        
        MyMatrix.putVal(0, 0,k);
        return MyMatrix;
    }
    
    public AbstractMatrix getK_kin() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);

        double E=0.0;
        if(this.ElementMaterial.getType()!=0){
            E=((ElastoPlastic) this.ElementMaterial).getKkinematic();
        }
        double A=this.theCrossSection.getA();
        double k=E*A;
        
        MyMatrix.putVal(0, 0,k);
        return MyMatrix;
    }
    
    private AbstractMatrix getM(int consistent) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double p=((ElasticMaterial) this.ElementMaterial).getDensity();
        double A=this.theCrossSection.getA();
        double z=1.;
        if(this.damage)z=this.getNodeHierarchy(1).getDamageVariableCur();
        double m=p*A*z;
        switch(consistent){
            case 0:
                MyMatrix.putVal(0, 0,m);
                break;
            default:
                MyMatrix.putVal(0, 0,m);
                break;
                
        }
        return MyMatrix;
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
            velcs.addVal(i, 0, theNode.getVelcsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return getM().times(velcs);
    }
    
    public AbstractMatrix getM_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i, 0, theNode.getAcclsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return getM().times(accls);
    }
    
    public AbstractMatrix getM_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return getM().times(disps);
    }
    
    public AbstractMatrix getK_u() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return getK().times(disps);
    }
    
    public AbstractMatrix getK_v() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return getK().times(velcs);
    }
    
    public AbstractMatrix getK_a() {
        AbstractMatrix accls = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            accls.addVal(i, 0, theNode.getAcclsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return getK().times(accls);
    }
    
    public int[] getFtable(){
        int[] theEFTable = new int[ndofs];
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theEFTable[j]=theNode.getFtable()[whichDOF-1];
            j=j+1;
        }
        return theEFTable;
    }
    
    @Override
    public AbstractMatrix getF(int LC, int step){
        double[] vec = new double[1];
        double[] vec2 = new double[1];
        double z=1.;
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getLoadCaseDisps(LC,step)[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
//        disps=this.getT().times(disps);

        //ex
        vec[0]=(disps.get(0, 0));
        vec2[0]=this.theMaterialPoint.giveStrainTakeStress(vec)[0];
        
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1, 0.0);
        Fint.addVal(0, 0, vec2[0]*this.theCrossSection.getA()*z);
        //Fint=this.getT().transpose().times(Fint);
        return Fint;
    }

    public AbstractMatrix getF() {
        //double[][] disps = new double[2][2];
        double[] vec = new double[1];
        double[] vec2 = new double[1];
        double z=1.;
        if(this.damage)z=this.getNodeHierarchy(1).getDamageVariableCur();
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
//        disps=this.getT().times(disps);
        //ex
        vec[0]=(disps.get(0, 0));
        vec2[0]=this.theMaterialPoint.giveStrainTakeStress(vec)[0];
        // find type (class) of material (elastoplastic or elastic_ simple material)
        /*if(this.ElementMaterial.getType()==0){
            this.Stress=((ElasticMaterial) this.ElementMaterial).getElasticity()*this.Etotal;
        }else{
            this.Stress=((ElasticMaterial) this.ElementMaterial).getElasticity()*this.Etotal;
        }*/
        AbstractMatrix Fint = new AbstractMatrix(ndofs,1, 0.0);
        Fint.addVal(0, 0, vec2[0]*this.theCrossSection.getA()*z);
        return Fint;
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
    
    public void setSens(boolean flag){
        this.sens=flag;
    }

    @Override
    public boolean ElementPlastified() {
        plast=false;
        if(theMaterialPoint.AskIfIsPlastified())plast=true;
        return plast;
    }

    @Override
//    public double getuKu() {
//        double[] vec = new double[1];
//        double val=0.;
//        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
//        AbstractMatrix plasts = new AbstractMatrix(ndofs,1);
//        int i=0;
//        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
//            Node theNode = it.next();
//            disps.addVal(i, 0, theNode.getDispsConvg()[this.whichDOF-1]);
//            plasts.addVal(i, 0, this.theMaterialPoint.getPlasticStrain()[0]);
//            i+=1;
//        }
//        // mipos edo den perno to sosto plastic strain !!!!!
////        vec[0]=(disps.get(0, 0));
////        double plastic = this.theMaterialPoint.giveStraingetPlasticStrain(vec)[0];
////        double plastic = this.theMaterialPoint.getPlasticStrain()[0];
////        disps.addVal(0, 0, -plastic);
//        val=(disps.minus(plasts)).transpose().times(getK().times((disps.minus(plasts)))).get(0, 0);
////        val+=disps.transpose().times(getK_kin().times(plastic)).get(0, 0);
//        
////        disps=this.getT().times(disps);
////        val = disps.transpose().times(getK().times(disps)).get(0, 0);
////        if(this.plast) val+=0.5*((ElastoPlastic) this.ElementMaterial).getYieldStress()[0]*((ElastoPlastic) this.ElementMaterial).getYieldStress()[0]/getK().get(0, 0);
//        return val;
//    }
    
    public double getuKu() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsConvg()[this.whichDOF-1]);
            i+=1;
        }
        if(this.nrg<0){
            nrg=this.getF().get(0, 0)*disps.get(0, 0);
            disp_pr=disps.get(0, 0);
            force_pr=this.getF().get(0, 0);
        }else{
            nrg+=(disps.get(0, 0)-disp_pr)*(getF().get(0, 0)+force_pr);
            disp_pr=disps.get(0, 0);
            force_pr=this.getF().get(0, 0);
        }
        return nrg*this.getNodeHierarchy(1).getDamageVariableCur();
    }
    
    @Override
    public double PlasticNRG(){
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        disps.addVal(0, 0, Math.abs(this.theMaterialPoint.getPlasticStrainInc()[0]));
        // from global to local disps
        double val=0.;
        if(this.ElementMaterial.getType()!=0){val=((ElastoPlastic) this.ElementMaterial).getYieldStress()[0]*disps.get(0, 0);}
        return val*this.theCrossSection.getA();
    }
    
    @Override
    public double getvMv() {
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            velcs.addVal(i, 0, theNode.getVelcsConvg()[this.whichDOF-1]);
            i+=1;
        }
        // from global to local disps
        return velcs.transpose().times(getM().times(velcs)).get(0, 0);
    }

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public double DamageDissipation(){
        double diss=((ElasticMaterial) this.ElementMaterial).getFracTough()*this.theCrossSection.getA()*(
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
            double Str=this.getF().get(0, 0);
            double StrCret=((ElasticMaterial)this.getMaterial()).getDamageStress();
            //System.out.println(Str + " "+StrCret);
            if(zprv>=1.0e-12 && Str>=StrCret){
                ElementNodes.get(1).setDamageVariableCur(0.0);
                answ=true;
                nrgvalcur = this.getuKu()/2. + this.DamageDissipation();
                if(nrgvalcur>=nrgvalprv){
                    this.getNodeHierarchy(1).setDamageVariableCur(zprv);
                    answ=false;
                }else{
    //                System.out.println("element "+id+" damaged");
                }
            }
        }
        return answ;
    }
    
    public AbstractMatrix getT() {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        MyMatrix.init(0.0);
        for(int i=0;i<ndofs;i++){
            MyMatrix.set(i, i, 1.0);
        }
        return MyMatrix;
    }
    
    public double getVolume(){
        return this.theCrossSection.getA();
    }
    
    @Override
    public double[] getNormal(int nid) {
        double[] normal = new double[3];
        normal[0]=0.0;
        normal[1]=0.0;
        normal[2]=0.0;
        return normal;
    }
    
    @Override
    public AbstractMatrix getM_upre() {
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            disps.addVal(i, 0, theNode.getDispsTrial()[0]-theNode.getDeltaDisps()[0]);
            i+=1;
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
            i+=1;
        }
        // from global to local disps
//        disps=this.getT().times(disps);
        return getK().times(disps);
    }
}
