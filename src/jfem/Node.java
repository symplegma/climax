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
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;
/**
 *
 * @author pchr
 */
public class Node  {
    private int id;
    private double[] coords ;
    private static int numberOfNodes = 0;
    static final int pdofs=7;
    private int dim;
    private static int ndofs=0;
    private int ndofs_ofNode=0;
    private int[] Ftable = new int[pdofs];
    private double[] dispsTrial;
    private double[] dispsConvg;
    private double[] velcsTrial;
    private double[] velcsConvg;
    private double[] acclsTrial;
    private double[] acclsConvg;
    private double[] accumulatedDisp;
    private double[] DeltaDisp;
    private Map<Integer,double[][]> dispsLC = new TreeMap<Integer,double[][]>();
    private Map<Integer,double[][]> velcsLC = new TreeMap<Integer,double[][]>();
    private double damageVariablePrv =1.0;
    private double damageVariableCur =1.0;
    private Map<Integer,Element> ConnectedElements = new TreeMap<Integer,Element>();
    
    private static double[] tempd;
    
    public Node(){}
    
    public Node(int id,double v[]){
        this.dim = v.length;
        this.coords = new double[3];
        for (int i=0;i<this.dim; i++){
            this.coords[i] = v[i];
        }
        if(v.length<2){this.coords[1]=0. ; this.coords[2]=0. ;}
        if(v.length<3){this.coords[2]=0. ;}
        this.id =id;
        ++numberOfNodes;
        this.dispsTrial= new double[pdofs];
        this.dispsConvg= new double[pdofs];
        this.velcsTrial= new double[pdofs];
        this.velcsConvg= new double[pdofs];
        this.acclsTrial= new double[pdofs];
        this.acclsConvg= new double[pdofs];
        this.accumulatedDisp= new double[pdofs];
        this.DeltaDisp= new double[pdofs];
        for(int i=0; i<pdofs; i++){
            Ftable[i]=0;
            accumulatedDisp[i]=0.;
            DeltaDisp[i]=0.;
        }
        
    }
    
    public Node(int id,double x, double y, double z){
        this.dim = 3;
        this.coords = new double[3];
        this.coords[0] = x;
        this.coords[1] = y;
        this.coords[2] = z;
        this.id =id;
        ++numberOfNodes;
        this.dispsTrial= new double[pdofs];
        this.dispsConvg= new double[pdofs];
        this.velcsTrial= new double[pdofs];
        this.velcsConvg= new double[pdofs];
        this.acclsTrial= new double[pdofs];
        this.acclsConvg= new double[pdofs];
        this.accumulatedDisp= new double[pdofs];
        this.DeltaDisp= new double[pdofs];
        for(int i=0; i<pdofs; i++){
            Ftable[i]=0;
            accumulatedDisp[i]=0.;
            DeltaDisp[i]=0.;
        }
        
    }
    
    public Node(int id,double x, double y){
        this.dim = 3;
        this.coords = new double[3];
        this.coords[0] = x;
        this.coords[1] = y;
        this.coords[2] = 0.0;
        this.id =id;
        ++numberOfNodes;
        this.dispsTrial= new double[pdofs];
        this.dispsConvg= new double[pdofs];
        this.velcsTrial= new double[pdofs];
        this.velcsConvg= new double[pdofs];
        this.acclsTrial= new double[pdofs];
        this.acclsConvg= new double[pdofs];
        this.accumulatedDisp= new double[pdofs];
        this.DeltaDisp= new double[pdofs];
        for(int i=0; i<pdofs; i++){
            Ftable[i]=0;
            accumulatedDisp[i]=0.;
            DeltaDisp[i]=0.;
        }
    }
    
    public Node(int id,double x){
        this.dim = 3;
        this.coords = new double[3];
        this.coords[0] = x;
        this.coords[1] = 0.0;
        this.coords[2] = 0.0;
        this.id =id;
        ++numberOfNodes;
        this.dispsTrial= new double[pdofs];
        this.dispsConvg= new double[pdofs];
        this.velcsTrial= new double[pdofs];
        this.velcsConvg= new double[pdofs];
        this.acclsTrial= new double[pdofs];
        this.acclsConvg= new double[pdofs];
        this.accumulatedDisp= new double[pdofs];
        this.DeltaDisp= new double[pdofs];
        for(int i=0; i<pdofs; i++){
            Ftable[i]=0;
            accumulatedDisp[i]=0.;
            DeltaDisp[i]=0.;
        }
    }
    
    public void setDamageVariablePrv(double val){this.damageVariablePrv=val;}
    
    public double getDamageVariablePrv(){return this.damageVariablePrv;}
    
    public void setDamageVariableCur(double val){this.damageVariableCur=val;}
    
    public double getDamageVariableCur(){return this.damageVariableCur;}
    
    // METHODS
    
    public void putElement(Element aElement) {
        this.ConnectedElements.put(aElement.getID(), aElement);
    }
    
    public Map getElements(){
        return this.ConnectedElements;
    }
    
    public int getNumElement(){
        return this.ConnectedElements.size();
    }
    
    public int getID() {
        return id;
    }
    
    public double[] getCoords() {
        return coords;
    }
    
    public double X() {
        return coords[0];
    }
    
    public double Y() {
        return coords[1];
    }
    
    public double Z() {
        return coords[2];
    }
    
    public int getDim() {
        return dim;
    }
    
    public int getNdofs_ofNode() {
        return ndofs_ofNode;
    }
    
    public int getNdofs() {
        return ndofs;
    }
    
    public int[] getFtable() {
        int[] theFtable = new int[this.Ftable.length];
        for(int i=0; i < this.Ftable.length; i++){
            theFtable[i]=this.Ftable[i];
        }
        return theFtable;
    }
    
    public int getnumberOfNodes() {
        return numberOfNodes;
    }
    
    public void setNdofs_ofNode(int[] dofs) {
        for(int i=0; i<dofs.length; i++ ){
            if(this.Ftable[dofs[i]-1]==0){
                this.Ftable[dofs[i]-1]=++ndofs;
                ++ndofs_ofNode;
            }
        }
    }
    
    public void updateDispsStep(double[] x, double coef){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.dispsTrial[i]+=coef*x[Ftable[i]-1];
                this.DeltaDisp[i]+=coef*x[Ftable[i]-1];
            }
        }
    }
    
    public void setDispsStepFromDerivatives(double coef_v,double coef_a){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.DeltaDisp[i]=coef_v*this.velcsTrial[i]+coef_a*this.acclsTrial[i];
            }
        }
    }
    
    public void updateDisps(double[] x){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.dispsTrial[i]=x[Ftable[i]-1];
            }
        }
    }
    
    public void accumulateDisps(double[] x){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.accumulatedDisp[i]+=x[Ftable[i]-1];
            }
        }
    }
    
    public void updateVelcs(double[] x){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.velcsTrial[i]=x[Ftable[i]-1];
            }
        }
    }
    
    public void updateVelcsStep(double[] x, double coef){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.velcsTrial[i]+=coef*x[Ftable[i]-1];
                //this.velcsTrial[i]=this.velcsTrial[i]+coef*x[Ftable[i]-1];
            }
        }
    }
    
    public void updateAccls(double[] x){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.acclsTrial[i]=x[Ftable[i]-1];
            }
        }
    }
    
    public void updateAcclsStep(double[] x, double coef){
        for(int i=0; i<pdofs; i++){
            if(Ftable[i]!=0){
                this.acclsTrial[i]+=coef*x[Ftable[i]-1];
                //this.acclsTrial[i]=this.acclsConvg[i]+coef*x[Ftable[i]-1];
            }
        }
    }
    
    public double[] getDispsConvg(){
        return this.dispsConvg;
    }
    
    public double[] getDeltaDisps(){
        return this.DeltaDisp;
    }
    
    public double[] getDispsAccum(){
        return this.accumulatedDisp;
    }
    
    public double[] getDispsTrial(){
        return this.dispsTrial;
    }
    
    public double[] getVelcsConvg(){
        return this.velcsConvg;
    }
    
    public double[] getVelcsTrial(){
        return this.velcsTrial;
    }
    
    public double[] getAcclsConvg(){
        return this.acclsConvg;
    }
    
    public double[] getAcclsTrial(){
        return this.acclsTrial;
    }
    
    ///////////////////////////////////////////////////////////
    public AbstractMatrix getDispsConvg_mat(){
        AbstractMatrix mat=new AbstractMatrix(dispsConvg.length,1);
    	for(int i=0; i<dispsConvg.length; i++){
            mat.set(i, 0, dispsConvg[i]);
        }
        return mat;
    }
    
    public AbstractMatrix getDispsTrial_mat(){
        AbstractMatrix mat=new AbstractMatrix(dispsTrial.length,1);
    	for(int i=0; i<dispsTrial.length; i++){
            mat.set(i, 0, dispsTrial[i]);
        }
        return mat;
    }
    
    public AbstractMatrix getVelcsConvg_mat(){
        AbstractMatrix mat=new AbstractMatrix(velcsConvg.length,1);
    	for(int i=0; i<velcsConvg.length; i++){
            mat.set(i, 0, velcsConvg[i]);
        }
        return mat;
    }    

    public AbstractMatrix getVelcsTrial_mat(){
        AbstractMatrix mat=new AbstractMatrix(velcsTrial.length,1);
    	for(int i=0; i<velcsTrial.length; i++){
            mat.set(i, 0, velcsTrial[i]);
        }
        return mat;
    }    
    
    public AbstractMatrix getAcclsConvg_mat(){
        AbstractMatrix mat=new AbstractMatrix(acclsConvg.length,1);
    	for(int i=0; i<acclsConvg.length; i++){
            mat.set(i, 0, acclsConvg[i]);
        }
        return mat;
    } 
    
    public AbstractMatrix getAcclsTrial_mat(){
        AbstractMatrix mat=new AbstractMatrix(acclsTrial.length,1);
    	for(int i=0; i<acclsTrial.length; i++){
            mat.set(i, 0, acclsTrial[i]);
        }
        return mat;
    } 
    ///////////////////////////////////////////////////////////
    
    
    public void commitDisps(){
        this.dispsConvg=this.dispsTrial;
    }
    
    public void printDisps(Analysis theAnalysis){
        theAnalysis.AnalysisFile.print(this.id);
        theAnalysis.AnalysisFile.print(" ");
        for(int i=0; i<pdofs; i++){
            theAnalysis.AnalysisFile.format("%16.8f", dispsConvg[i]);
            theAnalysis.AnalysisFile.print(" ");
        }
        theAnalysis.AnalysisFile.println();
    }
    
    public void commitVelcs(){
        this.velcsConvg=this.velcsTrial;
    }
    
    public void printVelcs(Analysis theAnalysis){
        theAnalysis.AnalysisFile.print(this.id);
        theAnalysis.AnalysisFile.print(" ");
        for(int i=0; i<pdofs; i++){
            theAnalysis.AnalysisFile.format("%16.8f", velcsConvg[i]);
            theAnalysis.AnalysisFile.print(" ");
        }
        theAnalysis.AnalysisFile.println();
    }

    public void commitAccls(){
        this.acclsConvg=this.acclsTrial;
    }
    
    public void printAccls(Analysis theAnalysis){
        theAnalysis.AnalysisFile.print(this.id);
        theAnalysis.AnalysisFile.print(" ");
        for(int i=0; i<pdofs; i++){
            theAnalysis.AnalysisFile.format("%16.8f", acclsConvg[i]);
            theAnalysis.AnalysisFile.print(" ");
        }
        theAnalysis.AnalysisFile.println();
    }
    
    public void clear(){
        for(int i=0; i<dispsTrial.length; i++){
            this.dispsTrial[i]=0.;
            this.dispsConvg[i]=0.;
            this.velcsTrial[i]=0.;
            this.velcsConvg[i]=0.;
            this.acclsTrial[i]=0.;
            this.acclsConvg[i]=0.;
        }
        
    }
    
    public void clearDeltaDisps(){
        for(int i=0; i<this.DeltaDisp.length; i++){
            this.DeltaDisp[i]=0.;
        }
    }
    
    public void clearVelcsTrial(){
        for(int i=0; i<velcsTrial.length; i++){
            this.velcsTrial[i]=0.;
        }
    }
    
    public void clearAcclsTrial(){
        for(int i=0; i<acclsTrial.length; i++){
            this.acclsTrial[i]=0.;
        }
    }
    
    public void clearAccumulatedDips(){
        for(int i=0; i<acclsTrial.length; i++){
            this.accumulatedDisp[i]=0.;
        }
    }
    
    public double getDisp(int wdof){
        return dispsConvg[wdof-1];
//        return this.DeltaDisp[wdof-1];
    }
    
    public double getVelc(int wdof){
        return velcsConvg[wdof-1];
    }
    
    public double getAccl(int wdof){
        return acclsConvg[wdof-1];
    }
    
    public AbstractMatrix getDispsAccumulated_mat(){
        AbstractMatrix mat=new AbstractMatrix(accumulatedDisp.length,1);
    	for(int i=0; i<accumulatedDisp.length; i++){
            mat.set(i, 0, accumulatedDisp[i]);
        }
        return mat;
    }

    public void print(){
        DecimalFormat Places = new DecimalFormat("0.000000000000000000");
        String vv;
        System.out.print(id+" ");
        for(int i=0; i<this.coords.length; i++){
            vv= Places.format(this.coords[i]);
            System.out.print(vv+" ");
        }
        System.out.println();
    }
    
    public void putLoadCase(int LC, int numSteps){
        double[][] convdispLC = new double[pdofs][numSteps+1];
        double[][] convvelcLC = new double[pdofs][numSteps+1];
        this.dispsLC.put(LC, convdispLC);
        this.velcsLC.put(LC, convvelcLC);
    }
    
    public void commitLoadCaseDisp(int LC, int step){
        for(int i=0;i<pdofs;i++){
            this.dispsLC.get(LC)[i][step]=this.dispsTrial[i];
            this.velcsLC.get(LC)[i][step]=this.velcsTrial[i];
        }
    }
    
    public double[] getLoadCaseDisps(int LC, int step){
        tempd = new double[pdofs];
        for(int i=0;i<pdofs;i++){
            tempd[i]=this.dispsLC.get(LC)[i][step];
        }
        return tempd;
    }
    
    public double[] getLoadCaseVelcs(int LC, int step){
        tempd = new double[pdofs];
        for(int i=0;i<pdofs;i++){
            tempd[i]=this.velcsLC.get(LC)[i][step];
        }
        return tempd;
    }
    
    public double[] getLoadCaseDisps(){return getLoadCaseDisps(1, 1);}
    
    public double[] getLoadCaseDisps(int LC){return getLoadCaseDisps(LC, 1);}
    
    public double[] getLoadCaseVelcs(){return getLoadCaseVelcs(1, 1);}
    
    public double[] getLoadCaseVelcs(int LC){return getLoadCaseVelcs(LC, 1);}
    
    public double[] getLoadCaseDispEvolution(int wdof, int LC){
        int len=this.dispsLC.get(LC)[wdof].length;
        double[] respevolution = new double[len];
        for(int i=0;i<len;i++)respevolution[i]=dispsLC.get(LC)[wdof][i];
        return respevolution;
    }
    
    public double[] getLoadCaseDispEvolution(int wdof){
        return getLoadCaseDispEvolution(wdof, 1);
    }
    
    public double[] getLoadCaseVelcEvolution(int wdof, int LC){
        int len=this.velcsLC.get(LC)[wdof].length;
        double[] respevolution = new double[len];
        for(int i=0;i<len;i++)respevolution[i]=velcsLC.get(LC)[wdof][i];
        return respevolution;
    }
    
    public double[] getLoadCaseVelcEvolution(int wdof){
        return getLoadCaseVelcEvolution(wdof, 1);
    }
    
    public double[] getLoadCaseEvolution(int wdof, int LC){
        if(wdof<=2){
            return getLoadCaseDispEvolution(wdof, LC);
        }else{
            return getLoadCaseVelcEvolution(wdof-3, LC);
        }
    }
    
    public void initDOFS(int n){ndofs=n;}
    public void initDOFS(){initDOFS(0);}
}
