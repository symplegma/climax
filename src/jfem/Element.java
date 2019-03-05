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
import java.util.Map;
import java.util.TreeMap;
import jmat.AbstractMatrix;
import java.util.Iterator;

/**
 *
 * @author pchr
 */
 abstract public class Element {
    protected int id;
    private static int numberOfElements = 0;
    //protected Vector ElementNodes;
    protected Map<Integer,Node> ElementNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,Integer> theNodesHierarchy = new TreeMap<Integer,Integer>();
    protected int numNodes=0;
    protected Material ElementMaterial;
    protected boolean consistentM=false;
    protected GaussData theGaussData;
    protected int numMP;
    protected CrossSection theCrossSection;
    protected boolean damage=false;
    protected int dof_per_node=0;
    protected int ndofs=0;
    protected int dimension=0;
    //protected boolean isPlastified=false;
    
    // constructors
    public Element(){}
        
    // methods
    public int getDimension(){return this.dimension;}
    public int get_dof_per_node(){return dof_per_node;}
    public void setDamage(boolean b){this.damage=b;}
    
    public boolean getDamage(){return this.damage;}
    
    public int getnumberOfElements() {
        return numberOfElements;
    }
    
    public int getID() {
        return id;
    }
    
    public int get_ndofs(){return this.ndofs;}
    
    abstract public boolean ElementPlastified();
    
    protected void putNode(Node aNode) {
        ElementNodes.put(++numNodes, aNode);
        aNode.putElement(this);
    }
    
    public void activateConsMass(){
        this.consistentM=true;
    }
    
    public void deactivateConsMass(){
        this.consistentM=false;
    }
    
    public boolean getConsMass(){
        return this.consistentM;
    }

    public int getNumNodes(){
        return this.ElementNodes.size();
    }

    public void print(){
        System.out.print("elem "+this.id+": ");
        for(Iterator<Node>it=this.ElementNodes.values().iterator();it.hasNext();){
            System.out.print(it.next().getID()+" ");
        }
//        System.out.println(((EBeam2d)this).getL());
        System.out.println(this.ElementMaterial.getID());
    }

    public Node getNodeHierarchy(int nodeHier) {
        return this.ElementNodes.get(nodeHier);
    }
    
    public int getHierOfNode(int theNode_id){
        int seq=0;
        int i=0;
        for(Iterator<Integer> it=this.theNodesHierarchy.values().iterator(); it.hasNext();){
            ++i;
            if(it.next()==theNode_id){seq=i;}
        }
        return seq;
    }

    public void setMaterial(Material Mater){
        this.ElementMaterial = Mater;
    }

    public Material getMaterial(){
        return this.ElementMaterial;
    }

    public void setCrossSection(CrossSection Sect){
        this.theCrossSection = Sect;
    }

    public CrossSection getCrossSection(){
        return this.theCrossSection;
    }
    
    public double DamageDissipation(){
        return 0.0;
    }
    
    public double PlasticNRG(){
        return 0.0;
    }
    
    public void commitDamage(){}
    
    public boolean checkForDamage(){return false;}
    
    public double getuCv(double am, double ak){
        double val=0.;
        
        AbstractMatrix disps = new AbstractMatrix(ndofs,1);
        AbstractMatrix velcs = new AbstractMatrix(ndofs,1);
        int i=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
//            if(dof_per_node<=0)System.err.println("warning element id="+id+" with "+dof_per_node+" dofs per node");
            for(int dn=0;dn<this.dof_per_node;dn++){
//                disps.addVal(i+dn, 0, theNode.getDispsConvg()[dn]);
                disps.addVal(i+dn, 0, theNode.getDeltaDisps()[dn]);
                velcs.addVal(i+dn, 0, theNode.getVelcsConvg()[dn]);
            }
            i+=dof_per_node;
        }
        // from global to local disps
        disps=this.getT().times(disps);
        velcs=this.getT().times(velcs);
//        val=velcs.transpose().times(( (getK().times(ak)).plus(getM().times(am)) ).times(velcs)).get(0, 0);
        val=disps.transpose().times(( (getK().times(ak)).plus(getM().times(am)) ).times(velcs)).get(0, 0);
        return val;
    }
    
    
    abstract public AbstractMatrix getK();
    abstract public AbstractMatrix getM();
    abstract int getElementNdofs();
    abstract int[] getFtable();
    abstract public AbstractMatrix getF();
    abstract void clear();
    abstract AbstractMatrix getM_v();
    abstract AbstractMatrix getM_a();
    abstract AbstractMatrix getM_u();
//    abstract AbstractMatrix getM_upre();
    abstract AbstractMatrix getK_u();
//    abstract AbstractMatrix getK_upre();
    abstract AbstractMatrix getK_v();
    abstract AbstractMatrix getK_a();
    abstract void commit();
    abstract public double getuKu();
    abstract public double getvMv();
    abstract public double[] getBVNorm(double coef);
    abstract public AbstractMatrix getT();
    abstract public double getVolume();
    abstract public double[] getNormal(int nid);
    
    public AbstractMatrix getFequivalent(LoadDist theLoad){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getDK(int param, int id){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getDM(int param, int id){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getDK(int param){
        return getDK(param, 0);
    }
    
    public AbstractMatrix getDM(int param){
        return getDM(param, 0);
    }
    
    public AbstractMatrix getF(int LC, int step){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getM_upre(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getK_upre(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double ShapeFunction(int wsf, double xsi, double eta){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double ShapeFunction(int wsf, double xsi){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double[] getStress(double[] at, int LC,int step){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getStressVonMisses(double[] at, int LC,int step){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getStressVonMisses(double[] at){
        return getStressVonMisses(at, 1,1);
    }
    
    public double[] getStress(double[] at){
        return getStress(at,1,1);
    }
    
    public double[] getCoordsAt(double[] at){
        double[] ret = null;
        switch(at.length){
            case 1:
                ret= new double[3]; ret[0]=0.0; ret[1]=0.0; ret[2]=0.0;
                for(int i=1;i<=this.getNumNodes();i++){
                    ret[0]+=this.getNodeHierarchy(i).getCoords()[0]*this.ShapeFunction(i,at[0]);
                    ret[1]+=this.getNodeHierarchy(i).getCoords()[1]*this.ShapeFunction(i,at[0]);
                }
                break;
            case 2:
                ret= new double[3]; ret[0]=0.0; ret[1]=0.0; ret[2]=0.0;
                for(int i=1;i<=this.getNumNodes();i++){
                    ret[0]+=this.getNodeHierarchy(i).getCoords()[0]*this.ShapeFunction(i,at[0], at[1]);
                    ret[1]+=this.getNodeHierarchy(i).getCoords()[1]*this.ShapeFunction(i,at[0], at[1]);
                }
                break;
            case 3:
                break;
        }
        return ret;
    }
}
