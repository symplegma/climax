/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

/**
 *
 * @author pchr
 */
public class ResultPoint extends geom.Point{
    private static int numPointInternals=0;
    private Domain onDomain=null;
    private boolean Boundary;
    private double[][] stress=null;
    
    // constructor
    public ResultPoint(){
    }
    
    public ResultPoint(int id, double[] coords){
        if(id>0){
            ++numPointInternals;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            for(int i=0; i<3; i++){
                if(i<coords.length){this.coordinates[i]=coords[i];}
                else{this.coordinates[i]=0.;}
            }
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public ResultPoint(int id, double x, double y){
        double[] coords = new double[2];
        coords[0]=x; coords[1]=y;
        if(id>0){
            ++numPointInternals;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            for(int i=0; i<3; i++){
                if(i<coords.length){this.coordinates[i]=coords[i];}
                else{this.coordinates[i]=0.;}
            }
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public ResultPoint(int id, double[] coords,Domain inDomain){
        if(id>0){
            ++numPointInternals;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            for(int i=0; i<3; i++){
                if(i<coords.length){this.coordinates[i]=coords[i];}
                else{this.coordinates[i]=0.;}
            }
            this.onDomain=inDomain;
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public ResultPoint(int id, double x, double y, Domain inDomain){
        double[] coords = new double[2];
        coords[0]=x; coords[1]=y;
        if(id>0){
            ++numPointInternals;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            for(int i=0; i<3; i++){
                if(i<coords.length){this.coordinates[i]=coords[i];}
                else{this.coordinates[i]=0.;}
            }
            this.onDomain=inDomain;
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public void setInDomain(Domain inDomain){this.onDomain=inDomain;}
    
    public Domain getDomain(){return this.onDomain;}
    
    public double getu(int wdisp, int step, int wstate){
        if(this.Boundary){
            return this.onDomain.getInternalPointDispOnBoundary(this.id, wdisp, step, wstate);
        }else{
            return this.onDomain.getInternalPointDispInterior(this.id, wdisp, step, wstate);
        }
    }
    
    public double getu(int wdisp, int step){
        return getu(wdisp,step, 0);
    }
    
    public double getu(int wdisp){
        return getu(wdisp,0, 0);
    }
    
    public double getstress(int wcomponent, int step, int wstate){
        if(this.Boundary){
            return this.onDomain.getInternalPointStressOnBoundary(this.id, wcomponent, step, wstate);
        }else{
            return this.onDomain.getInternalPointStress(this.id, wcomponent, step, wstate);
        }
    }
    
    public double getstressAux(int wdisp, int step, int wstate){
        if(this.Boundary){
            return this.onDomain.getInternalPointStressOnBoundaryAux(this.id, wdisp, step, wstate);
        }else{
            return this.onDomain.getInternalPointStressAux(this.id, wdisp, step, wstate);
        }
    }
    
    public double getstress(int wdisp, int step){
        return getstress(wdisp,step, 0);
    }
    
    public double getstress(int wdisp){
        return getstress(wdisp,0, 0);
    }
    
    public double getstrain(int wdisp, int step, int wstate){
        if(this.Boundary){
            return this.onDomain.getInternalPointStrainOnBoundary(this.id, wdisp, step, wstate);
        }else{
            return this.onDomain.getInternalPointStrain(this.id, wdisp, step, wstate);
        }
    }
    
    public double getstrain(int wdisp, int step){
        return getstrain(wdisp,step, 0);
    }
    
    public double getstrain(int wdisp){
        return getstrain(wdisp,0, 0);
    }
    
    public double getEnergyDensity(int step, int wstate){
        if(this.Boundary){
            return this.onDomain.getInternalPointEnergyDensityOnBoundary(this.id, step, wstate);
        }else{
            return this.onDomain.getInternalPointEnergyDensity(this.id, step, wstate);
        }
    }
    
    public double getEnergyDensity(int step){
        return getEnergyDensity(step, 0);
    }
    
    public double getEnergyDensity(){
        return getEnergyDensity(0, 0);
    }
    
    public double getTIEnergyDensityRate(double tau, int step, int wstate){
        // the Time Integral of the Energy Density Rate 
        if(this.Boundary){
            return this.onDomain.getInternalPointTIEnergyDensityRateOnBoundary(this.id, tau, step, wstate);
        }else{
            return this.onDomain.getInternalPointTIEnergyDensityRate(this.id, tau, step, wstate);
        }
    }
    
//    public double getTIEnergyDensityRate(double nu, double tau, int step){
//        // the Time Integral of the Energy Density Rate 
//        return getTIEnergyDensityRate(nu, tau, step, 0);
//    }
    
    public void setOnBoundary(){
        this.Boundary=true;
        //throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public boolean isOnBoundary(){
        return this.Boundary;
    }
}
