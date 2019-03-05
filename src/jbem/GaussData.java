/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import javax.swing.JOptionPane;

/**
 *
 * @author pchr
 */
public class GaussData {
    private int numofGaussPoints;
    private double[] gaussQuadDist;
    private double[] gaussQuadWeight;
    private double lower=-1.;
    private double upper=1.;
    private int type=1; //1: "Gauss-Legendre"
                        //2: "Gauss-Laguere"
    
    // constructor
    public GaussData(int numofGaussPoints){
        this.numofGaussPoints=numofGaussPoints;
        this.gaussQuadDist = new double[numofGaussPoints];
        this.gaussQuadWeight = new double[numofGaussPoints];
    }

    
    // methods
    public void computeGaussPoints(){
        switch(type){
            
            case 1:gaussQuadCoeff();break;
            case 2:gausslaguere();break;
            default: JOptionPane.showMessageDialog(null, "unkown type of Gauss quadrature",
                    "error: program will exit",JOptionPane.ERROR_MESSAGE);
            System.err.println("error (101) in:"+this.getClass().toString()+", method:"+"computeGaussPoints");
            break;
            
        }
    }
    
    /**
     * from http://www.ee.ucl.ac.uk/~mflanaga
     * WRITTEN BY: Dr Michael Thomas Flanagan
     */
    public void gaussQuadCoeff(){
        int n=this.numofGaussPoints;
        double	z=0.0D, z1=0.0D;
        double  pp=0.0D, p1=0.0D, p2=0.0D, p3=0.0D;

        double 	eps = 3e-11;	// set required precision
        double	x1 = lower;	// lower limit
        double	x2 = upper;	// upper limit

        //  Calculate roots
        // Roots are symmetrical - only half calculated
        int m  = (n+1)/2;
        double	xm = 0.5D*(x2+x1);
        double	xl = 0.5D*(x2-x1);

        // Loop for  each root
        for(int i=1; i<=m; i++){
            // Approximation of ith root
            z = Math.cos(Math.PI*(i-0.25D)/(n+0.5D));

            // Refinement on above using Newton's method
            do{
                p1 = 1.0D;
                p2 = 0.0D;

                // Legendre polynomial (p1, evaluated at z, p2 is polynomial of
                //  one order lower) recurrence relationsip
                for(int j=1; j<=n; j++){
                        p3 = p2;
                        p2 = p1;
                        p1= ((2.0D*j - 1.0D)*z*p2 - (j - 1.0D)*p3)/j;
                }
                pp = n*(z*p1 - p2)/(z*z - 1.0D);    // Derivative of p1
                z1 = z;
                z = z1 - p1/pp;
            } while(Math.abs(z - z1) > eps);

            gaussQuadDist[i-1] = xm - xl*z;                     // Scale root to desired interval
            gaussQuadDist[n-i] = xm + xl*z;                     // Symmetric counterpart
            gaussQuadWeight[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp);	// Compute weight
            gaussQuadWeight[n-i] = gaussQuadWeight[i-1];	// Symmetric counterpart
        }
    }
    
    public void gausslaguere(){
        if(this.numofGaussPoints>8){
            JOptionPane.showMessageDialog(null, "for Gauss Laguerre quadrature up to 8 points data are given",
                    "number of Guass point will be set equal to 8",JOptionPane.WARNING_MESSAGE);
            this.numofGaussPoints=8;
            this.gaussQuadDist = new double[numofGaussPoints];
            this.gaussQuadWeight = new double[numofGaussPoints];
        }
        switch(numofGaussPoints){
            case 1: gaussQuadDist[0]=0.5 ; gaussQuadWeight[0]=1.;
                    break;
            case 2: gaussQuadDist[0]=.112008806 ; gaussQuadWeight[0]=.718539319;
                    gaussQuadDist[1]=.602276908 ; gaussQuadWeight[1]=.281460680;
                    break;
            case 3: gaussQuadDist[0]=.063890793 ; gaussQuadWeight[0]=.513404552;
                    gaussQuadDist[1]=.368997063 ; gaussQuadWeight[1]=.391980041;
                    gaussQuadDist[2]=.766880303 ; gaussQuadWeight[2]=.0946154065;
                    break;
            case 4: gaussQuadDist[0]=.0414484801 ; gaussQuadWeight[0]=.3834640680;
                    gaussQuadDist[1]=.245274914 ; gaussQuadWeight[1]=.386875317;
                    gaussQuadDist[2]=.556165453 ; gaussQuadWeight[2]=.190435126;
                    gaussQuadDist[3]=.848982394 ; gaussQuadWeight[3]=.0392254871;
                    break;
            case 5: gaussQuadDist[0]=.0291344721 ; gaussQuadWeight[0]=.297893471;
                    gaussQuadDist[1]=.173977213 ; gaussQuadWeight[1]=.349776226;
                    gaussQuadDist[2]=.411702520; gaussQuadWeight[2]=.234488290;
                    gaussQuadDist[3]=.677314174 ; gaussQuadWeight[3]=.0989304595;
                    gaussQuadDist[4]=.894771361 ; gaussQuadWeight[4]=.0189115521;
                    break;
            case 6: gaussQuadDist[0]=.0216340058; gaussQuadWeight[0]=.238763662;
                    gaussQuadDist[1]=.129583391 ; gaussQuadWeight[1]=.308286573;
                    gaussQuadDist[2]=.314020449 ; gaussQuadWeight[2]=.245317426;
                    gaussQuadDist[3]=.538657217 ; gaussQuadWeight[3]=.142008756;
                    gaussQuadDist[4]=.756915337 ; gaussQuadWeight[4]=.0554546223;
                    gaussQuadDist[5]=.922668851 ; gaussQuadWeight[5]=.0101689586;
                    break;
            case 7: gaussQuadDist[0]=.0167193554; gaussQuadWeight[0]=.196169389;
                    gaussQuadDist[1]=.100185677 ; gaussQuadWeight[1]=.270302644;
                    gaussQuadDist[2]=.246294246 ; gaussQuadWeight[2]=.239681873;
                    gaussQuadDist[3]=.433463493 ; gaussQuadWeight[3]=.165775774;
                    gaussQuadDist[4]=.632350988 ; gaussQuadWeight[4]=.0889432271;
                    gaussQuadDist[5]=.811118626 ; gaussQuadWeight[5]=.0331943043;
                    gaussQuadDist[6]=.940848166 ; gaussQuadWeight[6]=.00593278701;
                    break;
            case 8: gaussQuadDist[0]=.0133202441 ; gaussQuadWeight[0]=.164416604;
                    gaussQuadDist[1]=.0797504290 ; gaussQuadWeight[1]=.237525610;
                    gaussQuadDist[2]=.197871029 ; gaussQuadWeight[2]=.226841984;
                    gaussQuadDist[3]=.354153994 ; gaussQuadWeight[3]=.175754079;
                    gaussQuadDist[4]=.529458575 ; gaussQuadWeight[4]=.112924030;
                    gaussQuadDist[5]=.701814529 ; gaussQuadWeight[5]=.0578722107;
                    gaussQuadDist[6]=.849379320 ; gaussQuadWeight[6]=.0209790737;
                    gaussQuadDist[7]=.953326450 ; gaussQuadWeight[7]=.00368640710;
                    break;
        }
        
    }
    
    public void print(){
        for(int i=0; i<numofGaussPoints; i++){
            System.out.print(gaussQuadDist[i]);
            System.out.print(" ");
            System.out.print(gaussQuadWeight[i]);
            System.out.println();
        }
    }
    
    public void setLimits(double lower, double upper){
        this.lower=lower;
        this.upper=upper;
    }
    
    public double getGaussCoordinate(int w){
        return this.gaussQuadDist[w];
    }
    
    public double getGaussWeight(int w){
        return this.gaussQuadWeight[w];
    }
    
    public void setGaussType(int w){
        this.type=w;
    }

}
