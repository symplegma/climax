/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package courses.structuraldynamics;

import climax.contraption;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.GeneralPath;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import mathman.LocalPowerDoubleFunction;
import mathman.DoubleFunction;
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
public class sdofSeries  implements FirstOrderDifferentialEquations, contraption{
    protected static int num_sdofSeries=0;
    private final double[] masses;
    private final double[] springs;
    private final double[] dashpots;
    private final double[] lengths; 
    DoubleFunction[] df;
    private double[] u0;
    private double[] v0;
    double[][] u,v,a;
    
    public double x0=0.0,y0=0.0;
    public double totLength;
    public double dx,dy, elen;
    protected int id;
    public Color SelfColor=Color.RED;
    public Color MotionColor=Color.ORANGE;
    
    EigenDecomposition eigen=null;
    Integer[] eig_index;
    
    private boolean Rayleigh=false;
    private double am,ak;
    
    public sdofSeries(double[] masses, double[] springs){
        num_sdofSeries++;
        id=num_sdofSeries;
        int len=springs.length;
        this.masses=new double[len-1];
        this.df=new DoubleFunction[len-1];
        this.springs=new double[len];
        this.dashpots=new double[len];
        this.lengths=new double[len];
        this.u0=new double[len];
        this.v0=new double[len];
        totLength=0.0;
        elen=Double.MAX_VALUE;
        for(int i=0;i<len;i++){
            if(i<len-1){
                this.masses[i]=masses[i];
                this.df[i]=new ZeroDoubleFunction();
            }
            this.springs[i]=springs[i];
            this.dashpots[i]=0.0;
            this.lengths[i]=1.0;
            totLength+=lengths[i];
            if(elen>lengths[i])elen=lengths[i];
            u0[i]=0.0;v0[i]=0.0;
            
        }
        dx = elen/20.0;
        dy=dx;
    }
    
    public sdofSeries(double[] masses, double[] springs, double[] dashpots){
        num_sdofSeries++;
        id=num_sdofSeries;
        int len=springs.length;
        this.masses=new double[len-1];
        this.df=new DoubleFunction[len-1];
        this.springs=new double[len];
        this.dashpots=new double[len];
        this.lengths=new double[len];
        this.u0=new double[len];
        this.v0=new double[len];
        totLength=0.0;
        elen=Double.MAX_VALUE;
        for(int i=0;i<len;i++){
            if(i<len-1){
                this.masses[i]=masses[i];
                this.df[i]=new ZeroDoubleFunction();
            }
            this.springs[i]=springs[i];
            this.dashpots[i]=dashpots[i];
            this.lengths[i]=1.0;
            totLength+=lengths[i];
            if(elen>lengths[i])elen=lengths[i];
            u0[i]=0.0;v0[i]=0.0;
        }
        dx = elen/20.0;
        dy=dx;
    }
    
    public void setLengths(double[] lengths){
        for(int i=0;i<this.lengths.length;i++){
            this.lengths[i]=lengths[i];
        }
    }
    
    public void setInitCond(double[] th0){
        u=new double[th0.length][1];
        this.u0=new double[th0.length]; 
        this.v0=new double[th0.length];
        for(int i=0;i<th0.length;i++){u0[i]=th0[i];u[i][0]=th0[i];v0[i]=0.0;}
    }
    
    public void setInitCond(double[] th0,double[] dth0){
        u=new double[th0.length][1];
        this.u0=new double[th0.length]; for(int i=0;i<th0.length;i++){this.u0[i]=th0[i];u[i][0]=th0[i];}
        this.v0=new double[dth0.length]; for(int i=0;i<dth0.length;i++)this.v0[i]=dth0[i];
    }
    
    public void setRayleigh(double am, double ak){
        this.Rayleigh=true;
        this.am=am;
        this.ak=ak;
    }
    
    @Override
    public int getID() {
        return this.id;
    }
    
    public double totLength(){
        return this.totLength;
    }
    
    // methods relevant to contraption and graphical representation
    @Override
    public void SelfPortrait(Graphics2D g2, double ymin, double ymax, int ye, int ys, double xmin, double xmax, int xe, int xs) {
        Function<Double, Integer> transX = x -> {
            return ((int) ((int) ((-xmin * xe + x * (xe - xs) + xmax * xs) / (xmax - xmin))));
        };
        
        Function<Double, Integer> transY = y -> {
            //int yd;
            return ((int) ((int) ((-ymin * ye + y * (ye - ys) + ymax * ys) / (ymax - ymin))));
            //return yd;
        };
        Color origColor=g2.getColor();
        g2.setColor(this.SelfColor);
        
        // top horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0), transX.apply(x0+dx/2.0),  transY.apply(y0));
        for(int i=0;i<=9;i++){
            double xi=x0-dx/2.0+i*dx/10;
            g2.drawLine(transX.apply(xi), transY.apply(y0), transX.apply(xi+dx/5.0), transY.apply(y0+dy/2.0));
        }
        
        double ystart=y0;
        double yend=y0;
        for(int i=0;i<this.lengths.length;i++){
            yend-=lengths[i];
            
            // top vertical dashpot and spring
            g2.drawLine(transX.apply(x0-dx/4.0), transY.apply(ystart), transX.apply(x0-dx/4.0),  transY.apply( ystart+(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/4.0), transY.apply(ystart), transX.apply(x0+dx/4.0),  transY.apply( ystart+(yend-ystart)/3.0 ));
            
            // dashpot
            g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 ), transX.apply(x0-dx/2.0+dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0+(yend-ystart)/6.0 ));
            g2.drawLine(transX.apply(x0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0+(yend-ystart)/6.0 ));

            g2.drawLine(transX.apply(x0-dx/4.0), transY.apply(ystart+(yend-ystart)/3.0+(yend-ystart)/10.0 ), transX.apply(x0-dx/4.0),  transY.apply(ystart+2.0*(yend-ystart)/3.0));
            g2.drawLine(transX.apply(x0-dx/4.0-dx/10.0), transY.apply(ystart+(yend-ystart)/3.0+(yend-ystart)/10.0 ), transX.apply(x0-dx/4.0+dx/10.0),  transY.apply(ystart+(yend-ystart)/3.0+(yend-ystart)/10.0 ));

            //spring
            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply(ystart+(yend-ystart)/3.0),                       transX.apply(x0+dx/8.0),         transY.apply(ystart+(yend-ystart)/3.0));

            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0),                      transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.1*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( ystart+(yend-ystart)/3.0 + 0.1*(yend-ystart)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 + 0.2*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.2*(yend-ystart)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.3*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( ystart+(yend-ystart)/3.0 + 0.3*(yend-ystart)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 + 0.4*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.4*(yend-ystart)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.5*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( ystart+(yend-ystart)/3.0 + 0.5*(yend-ystart)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ));

            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ),                       transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ));

            // bottom vertical dashpot and spring 
            g2.drawLine(transX.apply(x0-dx/4.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ), transX.apply(x0-dx/4.0),  transY.apply( ystart+(yend-ystart) ));
            g2.drawLine(transX.apply(x0+dx/4.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ), transX.apply(x0+dx/4.0),  transY.apply( ystart+(yend-ystart) ));

            // mass
            if(i<this.masses.length){
                GeneralPath gp = new GeneralPath();
                gp.moveTo(transX.apply(x0-dx), transY.apply(yend+dy));
                gp.lineTo(transX.apply(x0+dx), transY.apply(yend+dy));
                gp.lineTo(transX.apply(x0+dx), transY.apply(yend-dy));
                gp.lineTo(transX.apply(x0-dx), transY.apply(yend-dy));
                gp.closePath();
                g2.fill(gp);
                //g2.draw(gp);
            }
            ystart=yend;
        }
        
        // bottom horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0-totLength), transX.apply(x0+dx/2.0),  transY.apply(y0-totLength));
        for(int i=-1;i<=8;i++){
            double xi=x0-dx/2.0+i*dx/10;
            g2.drawLine(transX.apply(xi), transY.apply(y0-totLength-dy/2.0), transX.apply(xi+dx/5.0), transY.apply(y0-totLength));
        }
        
        g2.setColor(origColor);
    }

    @Override
    public double maximum_X_coordinate() {
        return x0+dx;
    }

    @Override
    public double minimum_X_coordinate() {
        return x0-dx;
    }

    @Override
    public double maximum_Y_coordinate() {
        return y0;
    }

    @Override
    public double minimum_Y_coordinate() {
        return y0-dy-this.totLength;
    }

    @Override
    public void Motion(Graphics2D g2, double ymin, double ymax, int ye, int ys, double xmin, double xmax, int xe, int xs, int step, double scale) {
        Function<Double, Integer> transX = x -> {
            return ((int) ((int) ((-xmin * xe + x * (xe - xs) + xmax * xs) / (xmax - xmin))));
        };
        
        Function<Double, Integer> transY = y -> {
            //int yd;
            return ((int) ((int) ((-ymin * ye + y * (ye - ys) + ymax * ys) / (ymax - ymin))));
            //return yd;
        };
        Color origColor=g2.getColor();
        g2.setColor(this.MotionColor);double[] disp= new double[this.masses.length];
        for(int i=0;i<disp.length;i++)disp[i]=this.u[i][step]*scale;
        
        // top horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0), transX.apply(x0+dx/2.0),  transY.apply(y0));
        for(int i=0;i<=9;i++){
            double xi=x0-dx/2.0+i*dx/10;
            g2.drawLine(transX.apply(xi), transY.apply(y0), transX.apply(xi+dx/5.0), transY.apply(y0+dy/2.0));
        }
        
        double ystart=y0;
        double yend=y0;
        double accumL=0.0;
        for(int i=0;i<this.lengths.length;i++){
            accumL-=lengths[i];
            if(i==0){
                ystart=y0;
                yend=-lengths[0]-disp[0]+dy;
            }
            if(i==(lengths.length-1)){
                ystart=accumL+lengths[i]-disp[i-1]-dy;
                yend=y0-this.totLength;
            }
            if(i>0 && i<(lengths.length-1)){
                ystart=accumL+lengths[i]-disp[i-1]-dy;
                yend=accumL-disp[i]+dy;
            }
            
            // top vertical dashpot and spring
            g2.drawLine(transX.apply(x0-dx/4.0), transY.apply(ystart), transX.apply(x0-dx/4.0),  transY.apply( ystart+(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/4.0), transY.apply(ystart), transX.apply(x0+dx/4.0),  transY.apply( ystart+(yend-ystart)/3.0 ));
            
            // dashpot
            g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0-dx/2.0+dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 ), transX.apply(x0-dx/2.0+dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0+(yend-ystart)/6.0 ));
            g2.drawLine(transX.apply(x0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 ), transX.apply(x0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0+(yend-ystart)/6.0 ));

            g2.drawLine(transX.apply(x0-dx/4.0), transY.apply(ystart+(yend-ystart)/3.0+(yend-ystart)/10.0 ), transX.apply(x0-dx/4.0),  transY.apply(ystart+2.0*(yend-ystart)/3.0));
            g2.drawLine(transX.apply(x0-dx/4.0-dx/10.0), transY.apply(ystart+(yend-ystart)/3.0+(yend-ystart)/10.0 ), transX.apply(x0-dx/4.0+dx/10.0),  transY.apply(ystart+(yend-ystart)/3.0+(yend-ystart)/10.0 ));

            //spring
            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply(ystart+(yend-ystart)/3.0),                       transX.apply(x0+dx/8.0),         transY.apply(ystart+(yend-ystart)/3.0));

            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0),                      transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.1*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( ystart+(yend-ystart)/3.0 + 0.1*(yend-ystart)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 + 0.2*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.2*(yend-ystart)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.3*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( ystart+(yend-ystart)/3.0 + 0.3*(yend-ystart)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 + 0.4*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.4*(yend-ystart)/3.0 ), transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.5*(yend-ystart)/3.0 ));
            g2.drawLine(transX.apply(x0+dx/8.0),        transY.apply( ystart+(yend-ystart)/3.0 + 0.5*(yend-ystart)/3.0 ), transX.apply(x0+dx/2.0-dx/8.0),  transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ));

            g2.drawLine(transX.apply(x0+dx/2.0-dx/8.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ),                       transX.apply(x0+dx/8.0),         transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ));

            // bottom vertical dashpot and spring 
            g2.drawLine(transX.apply(x0-dx/4.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ), transX.apply(x0-dx/4.0),  transY.apply( ystart+(yend-ystart) ));
            g2.drawLine(transX.apply(x0+dx/4.0), transY.apply( ystart+(yend-ystart)/3.0 + 0.6*(yend-ystart)/3.0 ), transX.apply(x0+dx/4.0),  transY.apply( ystart+(yend-ystart) ));

            // mass
            if(i<this.masses.length){
                GeneralPath gp = new GeneralPath();
                gp.moveTo(transX.apply(x0-dx), transY.apply(accumL+dy-disp[i]));
                gp.lineTo(transX.apply(x0+dx), transY.apply(accumL+dy-disp[i]));
                gp.lineTo(transX.apply(x0+dx), transY.apply(accumL-dy-disp[i]));
                gp.lineTo(transX.apply(x0-dx), transY.apply(accumL-dy-disp[i]));
                gp.closePath();
                g2.fill(gp);
                //g2.draw(gp);
            }
            //ystart=yend;
        }
        
        // bottom horizontal
        g2.drawLine(transX.apply(x0-dx/2.0), transY.apply(y0-totLength), transX.apply(x0+dx/2.0),  transY.apply(y0-totLength));
        for(int i=-1;i<=8;i++){
            double xi=x0-dx/2.0+i*dx/10;
            g2.drawLine(transX.apply(xi), transY.apply(y0-totLength-dy/2.0), transX.apply(xi+dx/5.0), transY.apply(y0-totLength));
        }
        
        g2.setColor(origColor);
    }
    
    public void setOrigin(double x0, double y0){
        this.x0=x0; this.y0=y0;
    }
    
    public RealMatrix getIMassC(){
        RealMatrix Cmat=this.getDamping();
        double[][] matrixData2 = new double[this.masses.length][this.masses.length];
        for(int i=0;i<this.masses.length;i++){
            for(int j=0;j<this.masses.length;j++){
                matrixData2[i][j]=Cmat.getEntry(i, j)/masses[i];
            }
        }
        RealMatrix IM = new Array2DRowRealMatrix(matrixData2);
        return IM;
    }
    
    public RealMatrix getIMassK(){
        RealMatrix Cmat=this.getStiffness();
        double[][] matrixData2 = new double[this.masses.length][this.masses.length];
        for(int i=0;i<this.masses.length;i++){
            for(int j=0;j<this.masses.length;j++){
                matrixData2[i][j]=Cmat.getEntry(i, j)/masses[i];
            }
        }
        RealMatrix IM = new Array2DRowRealMatrix(matrixData2);
        return IM;
    }
    
    public RealMatrix getMass(){
        double[][] matrixData2 = new double[this.masses.length][this.masses.length];
        for(int i=0;i<this.masses.length;i++){
            matrixData2[i][i]=masses[i];
        }
        RealMatrix IM = new Array2DRowRealMatrix(matrixData2);
        return IM;
    }
    
    public RealMatrix getStiffness(){
        double[][] matrixData2 = new double[this.masses.length][this.masses.length];
        matrixData2[0][0]=this.springs[0]+this.springs[1];
        if(this.masses.length>1){
            matrixData2[0][1]=-this.springs[1];
            for(int i=1;i<this.masses.length-1;i++){
                matrixData2[i][i]=this.springs[i]+this.springs[i+1];
                matrixData2[i][i-1]=-this.springs[i];
                matrixData2[i][i+1]=-this.springs[i+1];
            }
            matrixData2[masses.length-1][masses.length-1]=this.springs[masses.length-1]+this.springs[masses.length];
            matrixData2[masses.length-1][masses.length-2]=-this.springs[masses.length-1];
        }
        RealMatrix IM = new Array2DRowRealMatrix(matrixData2);
        return IM;
    }
    
    public RealMatrix getDamping(){
        RealMatrix IM = null;
        RealMatrix MassM=this.getMass();
        RealMatrix StifM=this.getStiffness();
        double[][] matrixData2;
        if(this.Rayleigh){
            matrixData2 = new double[this.masses.length][this.masses.length];
            for(int i=0;i<this.masses.length;i++){
                for(int j=0;j<this.masses.length;j++){
                    matrixData2[i][j]=am*MassM.getEntry(i, j)+ak*StifM.getEntry(i, j);
                }
            }
        }else{
            matrixData2 = new double[this.masses.length][this.masses.length];
            matrixData2[0][0]=this.dashpots[0]+this.dashpots[1];
            if(this.masses.length>1){
                matrixData2[0][1]=-this.dashpots[1];
                for(int i=1;i<this.masses.length-1;i++){
                    matrixData2[i][i]=this.dashpots[i]+this.dashpots[i+1];
                    matrixData2[i][i-1]=-this.dashpots[i];
                    matrixData2[i][i+1]=-this.dashpots[i+1];
                }
                matrixData2[masses.length-1][masses.length-1]=this.dashpots[masses.length-1]+this.dashpots[masses.length];
                matrixData2[masses.length-1][masses.length-2]=-this.dashpots[masses.length-1];
            }
        }
        IM = new Array2DRowRealMatrix(matrixData2);
        return IM;
    }
    
    public void solve(double tot, double dt){
        ClassicalRungeKuttaIntegrator CRK = new ClassicalRungeKuttaIntegrator(dt);
        ContinuousOutputModel tracker = new ContinuousOutputModel();
        CRK.addStepHandler(tracker);
        int N =(int) Math.round(tot/dt)+1;
        double[] y = new double[2*masses.length];
        for(int i=0;i<masses.length;i++){
            y[i]=this.u0[i];
            y[i+masses.length]=this.v0[i];
        }
        CRK.integrate(this, 0.0, y, tot, y);
        double h=dt;
        u=new double[masses.length][N];
        v=new double[masses.length][N];
        a=new double[masses.length][N];
        for(int it=0;it<N;it++){
            tracker.setInterpolatedTime(it*dt);
            double[] res=tracker.getInterpolatedState();
            double[] ares=tracker.getInterpolatedDerivatives();
            //System.out.println("Solution step = "+it);
            for(int i=0;i<masses.length;i++){
                u[i][it]=res[i];
                //System.out.println(res[i]);
                v[i][it]=res[i+masses.length];
                a[i][it]=ares[i+masses.length];
            }
        }
    }

    @Override
    public int getDimension() {
        return 2*this.masses.length;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        RealMatrix IMC=this.getIMassC();
        RealMatrix IMK=this.getIMassK();
        for(int i=0;i<masses.length;i++){
            yDot[i] = y[masses.length+i];
            yDot[masses.length+i] =df[i].run(t)/this.masses[i];
            for(int j=0;j<masses.length;j++){
                yDot[masses.length+i] -= y[masses.length+j]*IMC.getEntry(i, j)+y[j]*IMK.getEntry(i, j);
            }
        }
    }
    
    public void setRHS(DoubleFunction df){
        for(int i=0;i<this.df.length;i++)this.df[i]=df;
    }
    
    public void setRHS(DoubleFunction df, int dof){
        this.df[dof-1]=df;
    }
    
    public void setRHS(DoubleFunction[] df){
        for(int i=0;i<df.length;i++){
            this.df[i]=df[i];
        }
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
            System.out.println("WARNING! nDOFs= "+this.masses.length+" num of eig_vals= "+eigen.getRealEigenvalues().length+" num indices= "+eig_index.length);
        }
    }
    
    public double getEigenFrequency(int i){
        if(eigen==null){eigenAnalysis();}
        return Math.sqrt(eigen.getRealEigenvalue(eig_index[i-1]));
    }
    
    public double[] getEigenFrequencies(){
        if(eigen==null){eigenAnalysis();}
        double[] eigs = new double[this.masses.length];
        for(int i=0;i<eigs.length;i++)eigs[i]=Math.sqrt(eigen.getRealEigenvalue(eig_index[i]));
        return eigs;
    }
    
    public double[] getEigenVecArray(int i){
        double[] eigv=new double[this.masses.length];
        if(eigen==null){eigenAnalysis();}
        RealVector eigvec=eigen.getEigenvector(eig_index[i-1]);
        for(int j=0;j<eigv.length;j++){
            eigv[j]=eigvec.getEntry(j);
        }
        return eigv;
    }
    
    public RealVector getEigenVector(int i){
        if(eigen==null){eigenAnalysis();}
        RealVector eigvec=eigen.getEigenvector(eig_index[i-1]);
        return eigvec;
    }
    
    public double[] Disp(int w){
        int len = this.u[0].length;
        double[] disp = new double[len];
        for(int i=0;i<len;i++)disp[i]=u[w-1][i];
        return disp;
    }
    
    Complex[][] TransferFunction(double omext){
        // W. W. Smith, Jr. and S. Erdman, ‘A note on the inversion of complex matrices’, 
        // IEEE Trans. AC, AC-19,64 (1974)
        int NDOFS=this.masses.length;
        double[][] G = new double[2*NDOFS][2*NDOFS];
        RealMatrix K= this.getStiffness();
        RealMatrix M= this.getMass();
        RealMatrix C=this.getDamping();
        for(int i=0;i<NDOFS;i++){
            for(int j=0;j<NDOFS;j++){
                G[i][j]=K.getEntry(i, j)-omext*omext*M.getEntry(i, j); 
                G[i+NDOFS][j+NDOFS]=G[i][j];
                
                G[i][j+NDOFS]=omext*C.getEntry(i, j); 
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
        int NDOFS=this.masses.length;
        Complex[][][] TransferFun = new Complex[NDOFS][NDOFS][nom];
        for(int iom=0;iom<nom;iom++){
            double omext=om1+(om2-om1)*iom/(nom-1);
            double[][] G = new double[2*NDOFS][2*NDOFS];
            RealMatrix K= this.getStiffness();
            RealMatrix M= this.getMass();
            RealMatrix C=this.getDamping();
            for(int i=0;i<NDOFS;i++){
                for(int j=0;j<NDOFS;j++){
                    G[i][j]=K.getEntry(i, j)-omext*omext*M.getEntry(i, j); 
                    G[i+NDOFS][j+NDOFS]=G[i][j];

                    G[i][j+NDOFS]=omext*C.getEntry(i, j); 
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
        double[][][] ImpulseFun = new double[this.masses.length][this.masses.length][Nt+1];
        double[] temp_u0=null,temp_v0=null;
        if(u0!=null){
            temp_u0 = new double[this.u0.length];
        }
        if(v0!=null){
            temp_v0 = new double[this.v0.length];
        }
        u0=new double[this.masses.length];
        v0=new double[this.masses.length];
        for(int i=0;i<u0.length;i++){
            u0[i]=0.0;
            v0[i]=0.0;
        }
        DoubleFunction[] tempdf = new DoubleFunction[df.length];
        for(int i=0;i<tempdf.length;i++){
            tempdf[i]=df[i];
        }
        double h=tot/Nt;
        double amp=getMaxK();
        for(int ni=0;ni<this.masses.length;ni++){
            for(int i=0;i<this.masses.length;i++)df[i]=new ZeroDoubleFunction();
            LocalPowerDoubleFunction Dirac=new LocalPowerDoubleFunction(); Dirac.setAmplification(amp);
            Dirac.setTol(2.0*h);
            df[ni]= Dirac;
            this.solve(tot, h);
            for(int nj=0;nj<this.masses.length;nj++){
                for(int it=0;it<this.Disp(nj+1).length;it++){
                    ImpulseFun[ni][nj][it]=this.Disp(nj+1)[it];
                }
            }
        }
        
        for(int i=0;i<tempdf.length;i++){
            df[i]=tempdf[i];
        }
        u0=null;
        if(temp_u0!=null){
            u0 = new double[temp_u0.length];
        }
        for(int i=0;i<u0.length;i++){
            u0[i]=temp_u0[i];
        }
        v0=null;
        if(temp_v0!=null){
            v0 = new double[temp_v0.length];
        }
        for(int i=0;i<v0.length;i++){
            v0[i]=temp_v0[i];
        }
        
        return ImpulseFun;
    }
    
    public double getMaxK(){
        double val = Double.MIN_VALUE;
        for(int i=0;i<this.springs.length;i++){if(springs[i]>val)val=springs[i];}
        return val;
    }
}
