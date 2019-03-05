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

/* 
Boundary Elements in Dynamics
J. Dominguez
SUBROUTINE BESSEL(Z,FK0,FK1,FK2,FKOM,FK1M,FK2M)
С
С
С	THIS SUBROUTINE COMPUTES THE MOO[F°E3 3SSSEL FUNCTIONS OF THE SECCMO
С	KINO «NO ORDERS ZERO, ONE AND TUO CfltO, FK1.FK2) «NO COMPLEX ARCUME4T
С	Z USING A POL INCH IUM
С
С	THE COEFFICIENTS OF THE POLINOMIUM HAVE 3EEN OBTAINED FROM THE
С	TERMS OF THE SERIES EXPANSION
С
С	TUO DIFFEERENT SERIES HAVE BEEN USE). ONE FOR HOOULUSCZ).LT.5
С	ANO THE OTHER FOR MOOUIUS(Z).QT.5
С
С	FOR SMALL ARGUMENTS THE PARTS OF THE MOOIFIED BESSEL FUNCTIONS
С	TENDING TO ZERO WITH THE ARGUMENT ARE COMPUTED SEPARATEDLY FROM
С	THE BESSEL FUNCTIONS AND STORED IN FKOM.FKIM ANO FK2M. THESE
С	PARTS ARE ONLY USED BY THE QUADRATIC ELEMENT CODES IN ORDER
С	ТО COMPUTE THE DIFFERENCE BETWEEN THE DYNAMIC ANO THE STATIC
С	FUNDAMENTAL SOLUTION TRACTIONS
*/
package mathman;

import org.apache.commons.math3.complex.Complex;

/**
 *
 * @author pchr
 */
public class BesselK {
    
    protected static final double[] CI0P = {
   		1.00000000E+00, 
                0.62500000E+01,
                0.97656250E+01,
                0.67816840E+01, 
                0.26490953E+01, 
                0.66227383E+00,
                0.11497810E+00, 
                0.14665573E-01, 
                0.14321849E-02,
                0.11050809E-03, 
                0.69067559E-05, 
                0.35675392E-06,
                0.15484111E-07, 
                0.57263725E-09, 
                0.18260116E-10
   		};
    
    protected static final double[] CI1P = {
   		0.50000000E+00, 
                0.15625000E+01, 
                0.16276042E+01,
                0.84771050E+00, 
                0.26490953E+00, 
                0.55189486E-01,
                0.82127211E-02, 
                0.91659834E-03, 
                0.79565828E-04,
                0.55254047E-05, 
                0.31394345E-06, 
                0.14864747E-07,
                0.59554274E-09, 
                0.20451330E-10, 
                0.60867054E-12
   		};
    
    protected static final double[] CK0P = {
   		-0.57721564E+00, 
                0.26424021E+01, 
                0.90115658E+01,
                0.85185931E+01, 
                0.39898493E+01, 
                0.11299171E+01,
                0.21532918E+00, 
                0.29560538E-01, 
                0.30657944E-02,
                0.24883689E-03, 
                0.16242981E-04, 
                0.87142913E-06,
                0.39112788E-07, 
                0.14905278E-08, 
                0.48833881E-10
   		};
    
    protected static final double[] CK1P = {
   		1.00000000E+00, 
                0.96519581E+00,
                -0.26280638E+02,
                -0.44329875E+02, 
                -0.29269699E+02,
                -0.10636897E+02,
                -0.24689720E+01, 
                -0.39918196E+00, 
                -0.47620526E-01,
                -0.43685559E-02,
                -0.31795287E-03,
                -0.18814687E-04,
                -0.92322279E-06,
                -0.38181087E-07,
                -0.13490886E-08
   		};
    
    protected static final double[] CK0G = {
   		0.12533141E+01,
                -0.31332853E-01, 
                0.35249460E-02,
                -0.73436375E-03, 
                0.22489890E-03,
                -0.91084054E-04,
                0.45921544E-04,
                -0.27716932E-04, 
                0.19488468E-04,
                -0.15644909E-04, 
                0.14119530E-04,
                -0.14151620E-04
   		};
    
    protected static final double[] CK1G = {
   		0.12533141E+01, 
                0.93998560E-01,
                -0.58749100E-02,
                0.10281093E-02,
                -0.28915573E-03, 
                0.11132496E-03,
                -0.54270916E-04, 
                0.31981075E-04,
                -0.22086930E-04,
                0.17485486E-04,
                -0.15605797E-04, 
                0.15499393E-04
   		};
    
    static public Complex[] BesselKs(Complex Z){
        // FK0,FK1,FK2,FK0M,FK1M,FK2M
        Complex[] Ks= new Complex[6];
        Complex FI0M,FI1M,TempC;
        double R=Z.abs();
        double RR=5.0;
        if(R<RR){
            Complex ZZ=Z.divide(new Complex(RR,0.0));
            ZZ=ZZ.multiply(ZZ);
            int N=(int)Math.ceil(Math.max(2,Math.min(15.0,24.177/Math.log(RR/R))));
            
            TempC=SER(CI0P,ZZ,2,N);
            FI0M=new Complex(TempC.getReal(),TempC.getImaginary());
            TempC=SER(CI1P,ZZ,2,N);
            TempC=TempC.multiply(Z);
            FI1M=new Complex(TempC.getReal(),TempC.getImaginary());
            
            TempC=SER(CK0P,ZZ,2,N);
            TempC=TempC.add(FI0M.multiply(-1.0).multiply((Z.divide(2.0)).log()));
            Ks[3]=new Complex(TempC.getReal(),TempC.getImaginary());
            TempC=((Z.divide(2.0)).log().multiply(-1.0).add(new Complex(CK0P[0],0.0)));
            TempC=Ks[3].add(TempC);
            Ks[0]=new Complex(TempC.getReal(),TempC.getImaginary());
            
            TempC=SER(CK1P,ZZ,2,N);
            TempC=TempC.add(FI1M.multiply(Z).multiply((Z.divide(2.0)).log()));
            TempC=TempC.divide(Z);
            TempC=TempC.add((new Complex(CK0P[0]+0.5,0.0)).multiply(Z.divide(2.0)));
            Ks[4]=new Complex(TempC.getReal(),TempC.getImaginary());
            TempC=((new Complex(1.0,0.0)).divide(Z)).add(Ks[4]);
            TempC=TempC.add(Z.multiply(Z.divide(2.0).log()).divide(2.0));
            TempC=TempC.add(Z.multiply((new Complex(CK0P[0]+0.5,0.0)).multiply(-1.0)).divide(2.0));
            Ks[1]=new Complex(TempC.getReal(),TempC.getImaginary());
            
            TempC=Ks[4].multiply(new Complex(2.0,0.0)).divide(Z);
            TempC=TempC.add(Ks[3]);
            Ks[5]=new Complex(TempC.getReal(),TempC.getImaginary());
            TempC=Ks[1].multiply(new Complex(2.0,0.0)).divide(Z);
            TempC=TempC.add(Ks[0]);
            Ks[2]=new Complex(TempC.getReal(),TempC.getImaginary());
        }else{
            Complex ZZ=(new Complex(RR,0.0)).divide(Z);
            
            int N=(int)Math.ceil(Math.max(2,Math.min(12.0,12.088/Math.log(R/RR))));
            
            TempC=SER(CK0G,ZZ,2,N).add(new Complex(CK0G[0],0.0));
            TempC=TempC.multiply( Z.multiply(-1.0).exp() );
            TempC=TempC.divide(Z.sqrt());
            Ks[0]=new Complex(TempC.getReal(),TempC.getImaginary());
            TempC=(new Complex(CK0P[0],0.0)).multiply(-1.0);
            TempC=TempC.add(Z.divide(2.0).log());
            TempC=TempC.add(Ks[0]);
            Ks[3]=new Complex(TempC.getReal(),TempC.getImaginary());
            
            //FK1=(CMPLX(CK1G(1),0.)+SER(CK1G,ZZ,2,N))*CEXP(-Z)/CSQRT(Z)
            TempC=SER(CK1G,ZZ,2,N).add(new Complex(CK1G[0],0.0));
            TempC=TempC.multiply(Z.multiply(-1.0).exp());
            TempC=TempC.divide(Z.sqrt());
            Ks[1]=new Complex(TempC.getReal(),TempC.getImaginary());
            //FK1M=FK1-(1.,0.)/Z-Z*(CLOG(Z/(2.,0.))-CMPLX(CK0P(1)+0.5,0.))/2
            TempC=(new Complex(CK0P[0]+0.5,0.0)).multiply(-1.0);
            TempC=TempC.add((Z.divide(2.0)).log());
            TempC=TempC.multiply(Z).multiply(-0.5);
            TempC=TempC.add((new Complex(-1.0,0.0)).divide(Z));
            TempC=TempC.add(Ks[1]);
            Ks[4]=new Complex(TempC.getReal(),TempC.getImaginary());
            
            //FK2=FK0+FK1*(2.,0.)/Z
            TempC=Ks[1].multiply((new Complex(2.0,0.0)).divide(Z));
            TempC=TempC.add(Ks[0]);
            Ks[2]=new Complex(TempC.getReal(),TempC.getImaginary());
            //FK2M=FK0M+FK1M*(2.,0.)/Z
            TempC=Ks[4].multiply((new Complex(2.0,0.0)).divide(Z));
            TempC=TempC.add(Ks[3]);
            Ks[5]=new Complex(TempC.getReal(),TempC.getImaginary());
        }
        return Ks;
    }
    
    static private Complex SER(double[] S, Complex ZZ, int N1, int N2){
        /*
        THIS FUNCTION COMPUTES THE VALUE OF A SERIES OF N2-N1+1
        TERMS CONSISTING OF INCREASING POWERS OF ZZ TIMES THE
        COEFFICIENTS IN ARRAY S
        */
        Complex val,ZZZ;
        val = new Complex(0.0,0.0);
        ZZZ = new Complex(1.0,0.0);
        for(int i=N1;i<=N2;i++){
            ZZZ=ZZZ.multiply(ZZ);
            val=val.add(ZZZ.multiply(S[i-1]));
        }
        return val;
    }
    
    static public Complex INBESS(Complex Z){
        /*
        THIS FUNCTION COMPUTES THE INTEGRAL OF THE MODIFIED BESSEL FUNCTION
        OF ZERO ORDER BETWEEN 0 AND Z USING A POLINOMIUM.
        
        THE COEFFICIENTS OF THE POLINOMIUM HAVE BEEN OBTAINED FROM THE 
        TERMS OF THE SERIES EXPANSION
        */
        Complex FINKO = null;
        
        double S0 = 0.42278434;
        double[] S1 = {
            0.90734120E+01, 0.72756425E+02, 0.25901019E+03,
            0.52398212E+03, 0.68598066E+03, 0.62976149E+03,
            0.42826262E+03, 0.22451590E+03, 0.93540077E+02,
            0.31723213E+02, 0.89292258E+01, 0.21196877E+01,
            0.43013488E+00, 0.75474799E-01, 0.11565621E-01,
            0.15612002E-02, 0.18705605E-03, 0.20027866E-04,
            0.19277885E-05, 0.16772276E-06
            
        };
        double[] S2 = {
            0.12000000E+02, 0.64800000E+02, 0.18514286E+03,
            0.32400000E+03, 0.38173091E+03, 0.32300308E+03,
            0.20566727E+03, 0.10207750E+03, 0.40592223E+02,
            0.13221467E+02, 0.35916023E+01, 0.82606852E+00,
            0.16293265E+00, 0.27862515E-01, 0.41703893E-02,
            0.55091790E-03, 0.64704940E-04, 0.68008195E-05,
            0.64341868E-06, 0.55082916E-07
        };
        double R=Z.abs();
        double RR=12.0;
        int N;
        if(R-RR<0.0){
            N=(int)Math.ceil(Math.max(2,Math.min(20.0,19.572/Math.log(RR/R))));
        }else{
            N=20;
            System.err.println("INBESS: Results may be inaccurate");
        }
        Complex ZZ=Z.divide(new Complex(RR,0.0));
        ZZ=ZZ.multiply(ZZ);
        
        FINKO=(Z.divide(new Complex(2.,0.))).negate().log();
        FINKO=FINKO.add(SER(S2, ZZ, 1, N));
        FINKO=FINKO.add(SER(S1, ZZ, 1, N));
        FINKO=FINKO.add(new Complex(S0,0.0));
        FINKO=FINKO.multiply(Z);
        return FINKO;
    }
}