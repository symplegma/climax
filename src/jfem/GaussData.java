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
    
    // constructor
    public GaussData(int numofGaussPoints){
        this.numofGaussPoints=numofGaussPoints;
        this.gaussQuadDist = new double[numofGaussPoints];
        this.gaussQuadWeight = new double[numofGaussPoints];
        computeGaussPoints();
    }

    
    // methods
    private void computeGaussPoints(){
        gaussQuadCoeff();
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

}
