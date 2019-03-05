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
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import jmat.AbstractChDecomposition;
import jmat.AbstractMatrix;
import jmat.AbstractLUDecomposition;
import jmat.AbstractQRDecomposition;
import jmat.AbstractEigenvalueDecomposition;

/**
 *
 * @author pchr
 */
public class SOE {
    private AbstractMatrix K;
    private AbstractMatrix M;
    private AbstractMatrix Fint;
    private AbstractMatrix Fext;
    private AbstractMatrix x;
    private AbstractMatrix xu;
    private AbstractMatrix xv;
    private AbstractMatrix xa;
    private AbstractLUDecomposition LUofA ;
    private boolean isLUfactored=false;
    private AbstractQRDecomposition QRofA ;
    private boolean isQRfactored=false;
    private AbstractChDecomposition CDofA ;
    private boolean isCDfactored=false;
    private AbstractEigenvalueDecomposition EVofInvM_A ;
    double[] eig_vals;
    private int numDOFS=0;
    private Integer[] eig_index;
    
    public enum Solver {
    QR, LU, CD
    }
    static private Solver EquationSolver = Solver.LU;
    
    public static Solver getEquationSolve(){return SOE.EquationSolver;}
    
    public static void setEquationSolver(Solver theSolver){SOE.EquationSolver=theSolver;}

    //constructor
    public SOE(){
    }
    
    public SOE(int ndofs){
        this.numDOFS=ndofs;
        this.K = new AbstractMatrix(ndofs,ndofs);
        this.Fint = new AbstractMatrix(ndofs,1);
        this.Fext = new AbstractMatrix(ndofs,1);
        for (int i=0; i<ndofs; i++){
            Fint.set(i, 0, 0.0);
            Fext.set(i, 0, 0.0);
        }
        this.x = new AbstractMatrix(ndofs,1);
    }
    
    public SOE(int ndofs, int mdofs){
        this.numDOFS=ndofs;
        this.K = new AbstractMatrix(ndofs,ndofs);
        this.Fint = new AbstractMatrix(ndofs,1);
        this.Fext = new AbstractMatrix(ndofs,1);
        for (int i=0; i<ndofs; i++){
            Fint.set(i, 0, 0.0);
            Fext.set(i, 0, 0.0);
        }
        this.x = new AbstractMatrix(ndofs,1);
        this.xu = new AbstractMatrix(ndofs,1);
        this.xv = new AbstractMatrix(ndofs,1);
        this.xa = new AbstractMatrix(ndofs,1);
        
        this.M = new AbstractMatrix(ndofs,ndofs);
        this.eig_vals = new double[ndofs];
    }
    
    public void solve(){
        switch(EquationSolver){
            case LU:
                solveLU();
                break;
            case QR:
                solveQR();
                break;
            case CD:
                solveCD();
                break;
            default:
                System.err.println("error in SOE");
                System.exit(1);
                break;
        }
    }
    
    private void solveLU(){
        if(isLUfactored){
            //System.out.println("NO DECOMPOSITION");
             x=LUofA.solve(Fext.minus(Fint));
         } else {
            //System.out.println("WITH DECOMPOSITION");
             this.LUofA = new AbstractLUDecomposition(K);
             x=LUofA.solve(Fext.minus(Fint));
             isLUfactored=true;
         }
    }
    
    private void solveQR(){
        if(isQRfactored){
            //System.out.println("NO DECOMPOSITION");
             x=QRofA.solve(Fext.minus(Fint));
         } else {
            //System.out.println("WITH DECOMPOSITION");
             this.QRofA = new AbstractQRDecomposition(K);
             x=QRofA.solve(Fext.minus(Fint));
             isQRfactored=true;
         }
    }
    
    private void solveCD(){
       if(isCDfactored){
            //System.out.println("NO DECOMPOSITION");
             x=CDofA.solve(Fext.minus(Fint));
         } else {
            //System.out.println("WITH DECOMPOSITION");
             this.CDofA = new AbstractChDecomposition(K);
             x=CDofA.solve(Fext.minus(Fint));
             isCDfactored=true;
         }
    }
    
    
    public void addToK(AbstractMatrix mat, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             for(int j=0; j<s; j++){
                 K.addVal(theFtable[i]-1, theFtable[j]-1, coef*mat.get(i, j));
             }
         }
    }
    
    public void addToM(AbstractMatrix mat, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             for(int j=0; j<s; j++){
                 M.addVal(theFtable[i]-1, theFtable[j]-1, coef*mat.get(i, j));
             }
         }
    }

    public void addToFext(AbstractMatrix vec, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             Fext.addVal(theFtable[i]-1, 0, coef*vec.get(i, 0));
         }
    }
    
    public void rplToFext(AbstractMatrix vec, int[] theFtable){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             Fext.putVal(theFtable[i]-1, 0, vec.get(i, 0));
         }
    }
    
    public void addToFint(AbstractMatrix vec, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             Fint.addVal(theFtable[i]-1, 0, coef*vec.get(i, 0));
         }
    }
    
    public void rplToFint(AbstractMatrix vec, int[] theFtable){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             Fint.putVal(theFtable[i]-1, 0, vec.get(i, 0));
         }
    }    
    
    public double[] getX(){
        int order=this.x.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=x.get(i, 0);
         }
         return solution;
     }
    
    //==========================================
    public double[] getXv(){
        int order=this.xv.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=xv.get(i, 0);
         }
         return solution;
     }
    
    public void addToXv(AbstractMatrix vec, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             if(theFtable[i]!=0){xv.addVal(theFtable[i]-1, 0, coef*vec.get(i, 0));}
         }
    }
    
    public void rplToXv(AbstractMatrix vec, int[] theFtable){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             if(theFtable[i]!=0){xv.putVal(theFtable[i]-1, 0, vec.get(i, 0));}
         }
    }
    //==========================================
    public double[] getXa(){
        int order=this.xa.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=xa.get(i, 0);
         }
         return solution;
     }
    
    public void addToXa(AbstractMatrix vec, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             if(theFtable[i]!=0){xa.addVal(theFtable[i]-1, 0, coef*vec.get(i, 0));}
         }
    }
    
    public void rplToXa(AbstractMatrix vec, int[] theFtable){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             if(theFtable[i]!=0){xa.putVal(theFtable[i]-1, 0, vec.get(i, 0));}
         }
    }
    
    //==========================================
    public double[] getXu(){
        int order=this.xu.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=xu.get(i, 0);
         }
         return solution;
     }
    
    public void addToXu(AbstractMatrix vec, int[] theFtable, double coef){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             if(theFtable[i]!=0){xu.addVal(theFtable[i]-1, 0, coef*vec.get(i, 0));}
         }
    }
    
    public void rplToXu(AbstractMatrix vec, int[] theFtable){
         int s=theFtable.length;
         for(int i=0; i<s; i++){
             if(theFtable[i]!=0){xu.putVal(theFtable[i]-1, 0, vec.get(i, 0));}
         }
    }
    
    public double[] getEigenVector(int wEV){
        int order=this.x.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=K.get(i, eig_index[wEV]);
         }
         return solution;
     }
    
    public double getEigenValues(int wEV){
        return this.eig_vals[eig_index[wEV]];
    }
    
    public double[] getFext(){
        int order=this.Fext.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=Fext.get(i, 0);
         }
         return solution;
     }
    
    public double[] getFint(){
        int order=this.Fint.getRowDimension();
         double[] solution = new double[order];
         for(int i=0; i<order; i++){
             solution[i]=Fint.get(i, 0);
         }
         return solution;
     }

    
    public double[] getResidual(){
        int order=this.Fint.getRowDimension();
        double[] solution = new double[order];
        for(int i=0; i<order; i++){
             solution[i]=Fext.get(i, 0)-Fint.get(i, 0);
//             System.out.println(Fext.get(i, 0)+" "+Fint.get(i, 0));
         }
        return solution;
     }
    
    public void StaticRun(){
//        System.out.println("K");
//        K.print(12, 6);
//        System.out.println("Fext-Fint");
//        (Fext.minus(Fint)).print(12, 6);
        this.solve();
//         if(isLUfactored){
//             //x=LUofA.solve(Fext.minus(Fint));
//             x=LUofA.solve(Fext.minus(Fint));
//         } else {
//             this.LUofA = new AbstractLUDecomposition(K);
//             this.isLUfactored=true;
//             //x=LUofA.solve(Fext.minus(Fint));
//             x=LUofA.solve(Fext.minus(Fint));
//         }
//        System.out.println("Solution");
//        x.print(12, 6);
         
     }
    
    public void TransientRunInit(){
        switch(EquationSolver){
            case LU:
                this.LUofA = new AbstractLUDecomposition(M);
                x=LUofA.solve(Fext.minus(Fint));
                break;
            case QR:
                this.QRofA = new AbstractQRDecomposition(M);
                x=QRofA.solve(Fext.minus(Fint));
                break;
            case CD:
                this.CDofA = new AbstractChDecomposition(M);
                x=CDofA.solve(Fext.minus(Fint));
                break;
            default:
                System.err.println("error in SOE");
                System.exit(1);
                break;
        }
     }
    
    public void TransientRun(double ak, double am, double ac, double a0, double a1){
         //x=A.solve(b);
        switch(EquationSolver){
            case LU:
                if(isLUfactored){
                    x=LUofA.solve(Fext.minus(Fint));
                    //System.out.println(Fext.get(0, 0)+"  "+Fint.get(0, 0)+"  "+Fext.minus(Fint).get(0, 0)+"  "+(K.times(ak).plus(M.times(am))).get(0, 0)+ ": "+Fext.minus(Fint).get(0, 0)/(K.times(ak).plus(M.times(am))).get(0, 0));
                } else {
                    this.LUofA = new AbstractLUDecomposition( ((K.times(ak).plus(M.times(am))).plus(K.times(ac*a1))).plus(M.times(ac*a0)) );
                    this.isLUfactored=true;
                    x=LUofA.solve(Fext.minus(Fint));
                }
                break;
            case QR:
                if(isQRfactored){
                    x=QRofA.solve(Fext.minus(Fint));
                    //System.out.println(Fext.get(0, 0)+"  "+Fint.get(0, 0)+"  "+Fext.minus(Fint).get(0, 0)+"  "+(K.times(ak).plus(M.times(am))).get(0, 0)+ ": "+Fext.minus(Fint).get(0, 0)/(K.times(ak).plus(M.times(am))).get(0, 0));
                } else {
                    this.QRofA = new AbstractQRDecomposition( ((K.times(ak).plus(M.times(am))).plus(K.times(ac*a1))).plus(M.times(ac*a0)) );
                    this.isQRfactored=true;
                    x=QRofA.solve(Fext.minus(Fint));
                }
                break;
            case CD:
                if(isCDfactored){
                    x=CDofA.solve(Fext.minus(Fint));
                    //System.out.println(Fext.get(0, 0)+"  "+Fint.get(0, 0)+"  "+Fext.minus(Fint).get(0, 0)+"  "+(K.times(ak).plus(M.times(am))).get(0, 0)+ ": "+Fext.minus(Fint).get(0, 0)/(K.times(ak).plus(M.times(am))).get(0, 0));
                } else {
                    this.CDofA = new AbstractChDecomposition( ((K.times(ak).plus(M.times(am))).plus(K.times(ac*a1))).plus(M.times(ac*a0)) );
                    this.isCDfactored=true;
                    x=CDofA.solve(Fext.minus(Fint));
                }
                break;
            default:
                System.err.println("error in SOE");
                System.exit(1);
                break;
        }
        
    }
    
    public void TransientRun(double ak, double am){
        switch(EquationSolver){
            case LU:
                if(isLUfactored){
                    x=LUofA.solve(Fext.minus(Fint));
                    //System.out.println(Fext.get(0, 0)+"  "+Fint.get(0, 0)+"  "+Fext.minus(Fint).get(0, 0)+"  "+(K.times(ak).plus(M.times(am))).get(0, 0)+ ": "+Fext.minus(Fint).get(0, 0)/(K.times(ak).plus(M.times(am))).get(0, 0));
                } else {
                    this.LUofA = new AbstractLUDecomposition( K.times(ak).plus(M.times(am)) );
                    this.isLUfactored=true;
                    x=LUofA.solve(Fext.minus(Fint));
                    //System.out.println(Fext.get(0, 0)+"  "+Fint.get(0, 0)+"  "+Fext.minus(Fint).get(0, 0)+"  "+(K.times(ak).plus(M.times(am))).get(0, 0)+ ": "+Fext.minus(Fint).get(0, 0)/(K.times(ak).plus(M.times(am))).get(0, 0));
                    //System.out.println(K.times(ak).plus(M.times(am)).cond());
                }
                break;
            case QR:
                if(isQRfactored){
                    x=QRofA.solve(Fext.minus(Fint));
                } else {
                    this.QRofA = new AbstractQRDecomposition( K.times(ak).plus(M.times(am)) );
                    this.isQRfactored=true;
                    x=QRofA.solve(Fext.minus(Fint));
                }
                break;
            case CD:
                if(isCDfactored){
                    x=CDofA.solve(Fext.minus(Fint));
                } else {
                    this.CDofA = new AbstractChDecomposition( K.times(ak).plus(M.times(am)) );
                    this.isCDfactored=true;
                    x=CDofA.solve(Fext.minus(Fint));
                }
                break;
            default:
                System.err.println("error in SOE");
                System.exit(1);
                break;
        }
        
         //x=A.solve(b);
         
    }
    
    public void EigenRun(){
             M=M.inverse();
             K=M.times(K);
             this.EVofInvM_A = new AbstractEigenvalueDecomposition(K);
             eig_vals=EVofInvM_A.getRealEigenvalues();
//             for(int i=0;i<eig_vals.length;i++)System.out.println("eigenval "+(i+1)+" = "+eig_vals[i]);
//             Arrays.sort(eig_vals);
             // http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-get-a-sorted-list-of-indices-of-an-array
             Map<Double, Integer> map = new TreeMap<Double, Integer>();
             for (int i = 0; i < eig_vals.length; ++i) {
                 map.put(eig_vals[i], i);
             }
             Collection<Integer> indices = map.values();
             eig_index = indices.toArray(new Integer[0]);             
             if(eig_index.length!=eig_vals.length){
                 System.out.println("WARNING! nDOFs= "+this.numDOFS+" num of eig_vals= "+eig_vals.length+" num indices= "+eig_index.length);
             }
             K=EVofInvM_A.getV();
     }
    
    public int getEigIndexLength(){return eig_index.length;}
    
    public void clearFext(){
        for(int i=0; i<this.numDOFS; i++){
            Fext.set(i, 0, 0.);
        }
    }
    
    public void clearFint(){
        for(int i=0; i<this.numDOFS; i++){
            Fint.set(i, 0, 0.);
        }
    }
    
    public void clearX(){
        for(int i=0; i<this.numDOFS; i++){
            x.set(i, 0, 0.);
        }
    }
    
    public void clearXv(){
        for(int i=0; i<this.numDOFS; i++){
            xv.set(i, 0, 0.);
        }
    }
    
    public void clearXa(){
        for(int i=0; i<this.numDOFS; i++){
            xa.set(i, 0, 0.);
        }
    }
    
    public void clearK(){
        for(int i=0; i<this.numDOFS; i++){
            for(int j=0; j<this.numDOFS ; j++){
                K.set(i, j, 0.);
            }
        }
        isLUfactored=false;
    }
    
    public void clearM(){
        for(int i=0; i<this.numDOFS; i++){
            for(int j=0; j<this.numDOFS ; j++){
                M.set(i, j, 0.);
            }
        }
        isLUfactored=false;
    }
    
    public void clearXu(){
        for(int i=0; i<this.numDOFS; i++){
            xu.set(i, 0, 0.);
        }
    }
    
    public double getNormFext(){
        return Fext.norm1();
    }
    
    public void printK(){
        this.K.print(12, 12);
    }
    
    public void printM(){
        if(M==null){
            System.out.println("Mass matrix has not be defined");
        }else{
            this.M.print(12, 12);
        }
    }
    
    public AbstractMatrix getM(){return this.M;}
    
    public void printF(){
        this.Fext.print(12, 12);
    }    
}
