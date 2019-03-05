/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;
import jmat.AbstractMatrix;
import jmat.AbstractLUDecomposition;
import jmat.AbstractQRDecomposition;
/**
 *
 * @author pchr
 */
public class SOE {
    private AbstractMatrix A;
    private AbstractMatrix x;
    private AbstractMatrix b;
    private AbstractLUDecomposition LUofA ;
    private boolean isLUfactored=false;
    private AbstractQRDecomposition QRofA ;
    private boolean isQRfactored=false;
    private boolean Sparse=false;
    
    public enum Solver {
    QR, LU
    }
    static protected Solver EquationSolver = Solver.LU;
    
    public Solver getEquationSolve(){return SOE.EquationSolver;}
    
    public void setEquationSolver(Solver theSolver){SOE.EquationSolver=theSolver;}
    
    // constructor
    public SOE(){}
    
    public SOE(int order){
        this.A = new AbstractMatrix(order,order); this.A.init();
        this.x = new AbstractMatrix(order,1); this.x.init();
        this.b = new AbstractMatrix(order,1); this.b.init();
    }
    
    public SOE(int order, int steps){
        this.A = new AbstractMatrix(order,order); this.A.init();
        this.x = new AbstractMatrix(order,steps); this.x.init();
        this.b = new AbstractMatrix(order,1); this.b.init();
    }
    
    // methods
    public void init(int order){
        //System.out.println("The order of the system is :"+order+"x"+order);
        this.A = new AbstractMatrix(order,order); this.A.init();
        this.x = new AbstractMatrix(order,1); this.x.init();
        this.b = new AbstractMatrix(order,1); this.b.init();
    }
    
    public void init(int order, int steps){
        System.out.println("The order of the system is :"+order+"x"+order);
        this.A = new AbstractMatrix(order,order); this.A.init();
        this.x = new AbstractMatrix(order,steps); 
        //this.x.init(1.0);
        this.x.init();
        this.b = new AbstractMatrix(order,1); this.b.init();
    }
    
    public void solve(){
        switch(EquationSolver){
            case LU:
                solveLU();
                break;
            case QR:
                solveQR();
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
             x=LUofA.solve(b);
         } else {
            //System.out.println("WITH DECOMPOSITION");
             this.LUofA = new AbstractLUDecomposition(A);
             x=LUofA.solve(b);
             isLUfactored=true;
         }
    }
    
    private void solveQR(){
        if(isQRfactored){
            //System.out.println("NO DECOMPOSITION");
             x=QRofA.solve(b);
         } else {
            //System.out.println("WITH DECOMPOSITION");
             this.QRofA = new AbstractQRDecomposition(A);
             x=QRofA.solve(b);
             isQRfactored=true;
         }
    }

    public void setB(int row, double value){
        b.putVal(row, 0, value);
    }
    
    public void setA(int row, int clmn, double value){
        if(isLUfactored){isLUfactored=false;}
        A.putVal(row, clmn, value);
    }
    
    public void addA(int row, int clmn, double value){
        if(isLUfactored){isLUfactored=false;}
        A.addVal(row, clmn, value);
    }

    public void zeroA(){
        if(isLUfactored){isLUfactored=false;}
        A.init();
    }
    
    public void addB(int row, double value){
        b.addVal(row, 0, value);
    }
    
    public AbstractMatrix getB(){
        return this.b;
    }
    
    public AbstractMatrix getA(){
        return this.A;
    }
    
    public AbstractMatrix getX(){
        return this.x;
    }
    
    public void zeroB(){
        b.init();
    }
    
    public void setSparse(){this.Sparse=true;}
    
    public boolean getSparsity(){return this.Sparse;}
}
