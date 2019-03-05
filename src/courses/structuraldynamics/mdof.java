/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package courses.structuraldynamics;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import mathman.DoubleFunction;
import mathman.LocalPowerDoubleFunction;
import mathman.ZeroDoubleFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

/**
 *
 * @author pchr
 */
public class mdof implements FirstOrderDifferentialEquations{
    protected static int num_mdofs=0;
    double[][] K,M,C;
    double[] xis;
    int NDOFS,id;
    boolean diagonalMass=false; 
    
    DoubleFunction[] df;
    private double[] u0;
    private double[] v0;
    double[][] u,v,a;
    
    EigenDecomposition eigen=null;
    Integer[] eig_index;
    double tol=1.0e-2;
    int maxIter=100;
    
    public mdof(double[][] k, double[][] m, double[][] c){
        NDOFS=k.length;
        num_mdofs++;
        id=num_mdofs;
        K=new double[NDOFS][NDOFS];
        M=new double[NDOFS][NDOFS];
        C=new double[NDOFS][NDOFS];
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                K[i][j]=k[i][j];
                M[i][j]=m[i][j];
                C[i][j]=c[i][j];
            }
        }
    }
    
    public mdof(double[][] k, double[][] m){
        NDOFS=k.length;
        num_mdofs++;
        id=num_mdofs;
        K=new double[NDOFS][NDOFS];
        M=new double[NDOFS][NDOFS];
        C=new double[NDOFS][NDOFS];
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                K[i][j]=k[i][j];
                M[i][j]=m[i][j];
                C[i][j]=0.0;
            }
        }
    }
    
    public mdof(int ndofs){
        num_mdofs++;
        id=num_mdofs;
        NDOFS=ndofs;
        K=new double[NDOFS][NDOFS];
        M=new double[NDOFS][NDOFS];
        C=new double[NDOFS][NDOFS];
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                K[i][j]=0.0;
                M[i][j]=0.0;
                C[i][j]=0.0;
            }
        }
    }
    
    public void setC(double c[][]){
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                C[i][j]=c[i][j];
            }
        }
    }
    
    public void setC(double a0, double a1){
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                C[i][j]=a0*M[i][j]+a1*K[i][j];
            }
        }
    }
    
    public void setC(double[] xis){
        setC(xis, false);
    }
    
    public void setC(double[] xis, boolean computeC){
        this.xis = new double[NDOFS];
        for(int i=0;i<xis.length;i++)this.xis[i]=xis[i];
        if(computeC){
            if(this.eigen==null){
                eigenAnalysis();
            }
            
            double[][] Phi = new double[NDOFS][NDOFS];
            double[][] Cdiag=new double[NDOFS][NDOFS];
            for(int i=0;i<NDOFS;i++){
                double sqrt_eigmass=Math.sqrt(getEigenMass(i+1));
                for(int j=0;j<NDOFS;j++){
                    Phi[j][i]=(eigen.getEigenvector(eig_index[i])).getEntry(j)/sqrt_eigmass;
                }
                Cdiag[i][i]=2.0*xis[i]*Math.sqrt(eigen.getRealEigenvalue(eig_index[i]));
            }
            RealMatrix PhiMat=new Array2DRowRealMatrix(Phi);
            //LUDecomposition LUD= new LUDecomposition(PhiMat);
            //RealMatrix IPhiMat=LUD.getSolver().getInverse();
            RealMatrix MassMat=new Array2DRowRealMatrix(M);
            RealMatrix CdiagMat=new Array2DRowRealMatrix(Cdiag);
            C=(MassMat.multiply(PhiMat).multiply(CdiagMat).multiply(PhiMat.transpose()).multiply(MassMat)).getData();
        }
    }
    
    public double getEigenMass(int i){
        if(this.eigen==null){
            eigenAnalysis();
        }
        RealMatrix phi = new Array2DRowRealMatrix(getEigenVecArray(i));
        RealMatrix Mmat = new Array2DRowRealMatrix(M);
        return phi.transpose().multiply(Mmat.multiply(phi)).getEntry(0, 0);
    }
    
    public void eigenAnalysis(){
        eigen = new EigenDecomposition(getIMassK());
        eig_index=new Integer[eigen.getRealEigenvalues().length];
        Map<Double, Integer> map = new TreeMap<Double, Integer>();
        for (int i = 0; i < eigen.getRealEigenvalues().length; ++i) {
            map.put(eigen.getRealEigenvalues()[i], i);
        }
        Collection<Integer> indices = map.values();
        eig_index = indices.toArray(new Integer[0]);             
        if(eig_index.length!=eigen.getRealEigenvalues().length){
            System.out.println("WARNING! nDOFs= "+NDOFS+" num of eig_vals= "+eigen.getRealEigenvalues().length+" num indices= "+eig_index.length);
        }
    }
    
    public RealMatrix getIMassK(){
        RealMatrix Kmat=new Array2DRowRealMatrix(K);
        double[][] matrixData2 = new double[NDOFS][NDOFS];
        RealMatrix IM;
        if(this.diagonalMass){
            for(int i=0;i<NDOFS;i++){
                for(int j=0;j<NDOFS;j++){
                    matrixData2[i][j]=Kmat.getEntry(i, j)/M[i][i];
                }
            }
            IM = new Array2DRowRealMatrix(matrixData2);
        }else{
            RealMatrix Mmat=new Array2DRowRealMatrix(M);
            LUDecomposition LUD= new LUDecomposition(Mmat);
            Mmat=LUD.getSolver().getInverse();
            IM = Mmat.multiply(Kmat);
        }
        return IM;
    }
    
    public RealMatrix getIK_M(double mu){
        RealMatrix Mmat=new Array2DRowRealMatrix(M);
        RealMatrix IK;
        RealMatrix Kmat=new Array2DRowRealMatrix(K);
        LUDecomposition LUD= new LUDecomposition(Kmat.add(Mmat.scalarMultiply(-mu)));
        Kmat=LUD.getSolver().getInverse();
        IK = Kmat.multiply(Mmat);
        return IK;
    }
    
    public RealMatrix getIMassC(){
        RealMatrix Kmat=new Array2DRowRealMatrix(C);
        double[][] matrixData2 = new double[NDOFS][NDOFS];
        RealMatrix IM;
        if(this.diagonalMass){
            for(int i=0;i<NDOFS;i++){
                for(int j=0;j<NDOFS;j++){
                    matrixData2[i][j]=Kmat.getEntry(i, j)/M[i][i];
                }
            }
            IM = new Array2DRowRealMatrix(matrixData2);
        }else{
            RealMatrix Mmat=new Array2DRowRealMatrix(M);
            LUDecomposition LUD= new LUDecomposition(Mmat);
            Mmat=LUD.getSolver().getInverse();
            IM = Mmat.multiply(Kmat);
        }
        return IM;
    }
    
    public RealMatrix getIMass(){
        RealMatrix IM;
        double[][] matrixData2 = new double[NDOFS][NDOFS];
        if(this.diagonalMass){
            for(int i=0;i<NDOFS;i++){
                for(int j=0;j<NDOFS;j++){
                    matrixData2[i][j]=0.0;
                    if(i==j)matrixData2[i][j]=1.0/M[i][i];
                }
            }
            IM = new Array2DRowRealMatrix(matrixData2);
        }else{
            RealMatrix Mmat=new Array2DRowRealMatrix(M);
            LUDecomposition LUD= new LUDecomposition(Mmat);
            IM=LUD.getSolver().getInverse();
        }
        return IM;
    }
    
    public double getEigenValue(int i){
        if(eigen==null){eigenAnalysis();}
        return eigen.getRealEigenvalue(eig_index[i-1]);
    }
    
    public double[] getEigenValues(){
        if(eigen==null){eigenAnalysis();}
        double[] eigs = new double[NDOFS];
        for(int i=0;i<eigs.length;i++)eigs[i]=eigen.getRealEigenvalue(eig_index[i]);
        return eigs;
    }
    
    public double[] getEigenVecArray(int i){
        double[] eigv=new double[NDOFS];
        if(eigen==null){eigenAnalysis();}
        RealVector eigvec=eigen.getEigenvector(eig_index[i-1]);
        for(int j=0;j<eigv.length;j++){
            eigv[j]=eigvec.getEntry(j);
        }
        return eigv;
    }
    
    public double[][] getEigenVecArray(){
        double[][] eigv=new double[NDOFS][NDOFS];
        if(eigen==null){eigenAnalysis();}
        for(int i=0;i<this.NDOFS;i++){
            RealVector eigvec=eigen.getEigenvector(eig_index[i-1]);
            for(int j=0;j<eigv.length;j++){
                eigv[i][j]=eigvec.getEntry(j);
            }
        }
        return eigv;
    }
    
    public RealVector getEigenVector(int i){
        if(eigen==null){eigenAnalysis();}
        RealVector eigvec=eigen.getEigenvector(eig_index[i-1]);
        return eigvec;
    }
    
    public double getEigenValueIm(int i){
        if(eigen==null){eigenAnalysis();}
        return eigen.getImagEigenvalue(eig_index[i-1]);
    }
    
    public double[] getEigenValuesIm(){
        if(eigen==null){eigenAnalysis();}
        double[] eigs = new double[NDOFS];
        for(int i=0;i<eigs.length;i++)eigs[i]=eigen.getImagEigenvalue(eig_index[i]);
        return eigs;
    }
    
    public void setDiagonalMass(boolean b){this.diagonalMass=b;}
    
    Complex[][] TransferFunction(double omext){
        // W. W. Smith, Jr. and S. Erdman, ‘A note on the inversion of complex matrices’, 
        // IEEE Trans. AC, AC-19,64 (1974)
        double[][] G = new double[2*NDOFS][2*NDOFS];
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                G[i][j]=K[i][j]-omext*omext*M[i][j]; 
                G[i+NDOFS][j+NDOFS]=G[i][j];
                
                G[i][j+NDOFS]=omext*C[i][j];
                G[i+NDOFS][j]=-G[i][j+NDOFS];
            }
        }
        RealMatrix Gmat=new Array2DRowRealMatrix(G);
        G=null;
        LUDecomposition LUD= new LUDecomposition(Gmat);
        Gmat=LUD.getSolver().getInverse();
        
        // create a 2x2 complex matrix
        Complex[][] TransferFun = new Complex[NDOFS][NDOFS];
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                TransferFun[i][j]=new Complex(Gmat.getEntry(i, j),Gmat.getEntry(i, j+NDOFS));
            }
        }
        return TransferFun;
    }
    
    Complex[][][] TransferFunction(double om1, double om2, int nom){
        Complex[][][] TransferFun = new Complex[NDOFS][NDOFS][nom];
        for(int iom=0;iom<nom;iom++){
            double omext=om1+(om2-om1)*iom/(nom-1);
            double[][] G = new double[2*NDOFS][2*NDOFS];
            for(int i=0;i<NDOFS;i++){
                for(int j=0;j<NDOFS;j++){
                    G[i][j]=K[i][j]-omext*omext*M[i][j]; 
                    G[i+NDOFS][j+NDOFS]=G[i][j];

                    G[i][j+NDOFS]=omext*C[i][j]; 
                    G[i+NDOFS][j]=-G[i][j+NDOFS];
                }
            }
            RealMatrix Gmat=new Array2DRowRealMatrix(G);
            G=null;
            LUDecomposition LUD= new LUDecomposition(Gmat);
            Gmat=LUD.getSolver().getInverse();
            
            for(int i=0;i<NDOFS;i++){
                for(int j=0;j<NDOFS;j++){
                    TransferFun[i][j][iom]=new Complex(Gmat.getEntry(i, j),Gmat.getEntry(i, j+NDOFS));
                }
            }
        }
        return TransferFun;
    }
    
    double[][][] ImpulseFunction(double tot, int Nt){
        double[][][] ImpulseFun = new double[NDOFS][NDOFS][Nt+1];
        u0=new double[NDOFS];
        v0=new double[NDOFS];
        df = new DoubleFunction[NDOFS];
        for(int i=0;i<u0.length;i++){
            u0[i]=0.0;
            v0[i]=0.0;
        }
        
        double h=tot/Nt;
        double amp=getMaxK();
        for(int ni=0;ni<NDOFS;ni++){
            for(int i=0;i<NDOFS;i++)df[i]=new ZeroDoubleFunction();
            LocalPowerDoubleFunction Dirac=new LocalPowerDoubleFunction(); Dirac.setAmplification(amp);
            Dirac.setTol(2.0*h);
            df[ni]= Dirac;
            this.solve(tot, h);
            for(int nj=0;nj<NDOFS;nj++){
                for(int it=0;it<this.Disp(nj+1).length;it++){
                    ImpulseFun[ni][nj][it]=this.Disp(nj+1)[it];
                }
            }
        }
        return ImpulseFun;
    }
    
    public double[] Disp(int w){
        int len = this.u[0].length;
        double[] disp = new double[len];
        for(int i=0;i<len;i++)disp[i]=u[w-1][i];
        return disp;
    }
    
    public double getMaxK(){
        double val = Double.MIN_VALUE;
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<=i;j++){
            if(K[i][j]>val)val=K[i][j];}
        }
        return val;
    }

    @Override
    public int getDimension() {
        return 2*NDOFS;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        RealMatrix IM=this.getIMass();
        RealMatrix IMC=IM.multiply(new Array2DRowRealMatrix(C));
        RealMatrix IMK=IM.multiply(new Array2DRowRealMatrix(K));
        for(int i=0;i<NDOFS;i++){
            yDot[i] = y[NDOFS+i];
            yDot[NDOFS+i] =0.0;
            for(int j=0;j<NDOFS;j++){
                DoubleFunction fex=new ZeroDoubleFunction();
                if(df[i]!=null)fex=df[i];
                yDot[NDOFS+i] += fex.run(t)*IM.getEntry(i, j)-y[NDOFS+j]*IMC.getEntry(i, j)-y[j]*IMK.getEntry(i, j);
            }
        }
    }
    
    public void solve(double tot, double dt){
        ClassicalRungeKuttaIntegrator CRK = new ClassicalRungeKuttaIntegrator(dt);
        ContinuousOutputModel tracker = new ContinuousOutputModel();
        CRK.addStepHandler(tracker);
        int N =(int) Math.round(tot/dt)+1;
        double[] y = new double[2*NDOFS];
        for(int i=0;i<NDOFS;i++){
            y[i]=this.u0[i];
            y[i+NDOFS]=this.v0[i];
        }
        CRK.integrate(this, 0.0, y, tot, y);
        double h=dt;
        u=new double[NDOFS][N];
        v=new double[NDOFS][N];
        a=new double[NDOFS][N];
        for(int it=0;it<N;it++){
            tracker.setInterpolatedTime(it*dt);
            double[] res=tracker.getInterpolatedState();
            double[] ares=tracker.getInterpolatedDerivatives();
            //System.out.println("Solution step = "+it);
            for(int i=0;i<NDOFS;i++){
                u[i][it]=res[i];
                //System.out.println(res[i]);
                v[i][it]=res[i+NDOFS];
                a[i][it]=ares[i+NDOFS];
            }
        }
    }
    public double[] InvVecIter(){
        return InvVecIter(0.0);
    }
    
    public double[] InvVecIter(double mu){
        double[] x= new double[NDOFS];
        for(int i=0;i<NDOFS;i++)x[i]=1.0;
        RealMatrix Xp=new Array2DRowRealMatrix(x);
        RealMatrix Xn;
        RealMatrix IKM=getIK_M(mu);
        RealMatrix Mmat=new Array2DRowRealMatrix(M);
        double lamp=1.0; double lamn = 0;
        double norm=Double.POSITIVE_INFINITY;
        int count=0; int count2=0;
        while(count2<maxIter){
            if(norm>tol){
                count++;
                count2=count;
                Xn=IKM.multiply(Xp);
                double numerator =(Xn.transpose().multiply(Mmat).multiply(Xp)).getEntry(0, 0);
                double denominator=(Xn.transpose().multiply(Mmat).multiply(Xn)).getEntry(0, 0);
                lamn=numerator/denominator;
                norm=Math.abs((lamn-lamp)/lamn);
                double scale=Math.sqrt((Xn.transpose().multiply(Mmat).multiply(Xn)).getEntry(0, 0));
                for(int i=0;i<NDOFS;i++)x[i]=Xn.getEntry(i, 0)/scale;
                Xp=new Array2DRowRealMatrix(x); // this Xp is actually the approximation of eigenvector
                //System.out.print(count+"/"+norm+": "+lamn+", ");
                //for(int i=0;i<NDOFS;i++)System.out.print(x[i]+" ");
                //System.out.println();
            }else{
                count2=this.maxIter+1;
            }
        }
        System.out.println("Eigenvalue = "+(lamn+mu)+" ("+count+"/"+norm+")");
        return Xp.getColumn(0);
    }
}