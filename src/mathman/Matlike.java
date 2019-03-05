/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mathman;

import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.complex.Complex;

/**
 *
 * @author pchr
 */
public final class Matlike {
    
    public static int nextpow2(double A){
        int nextN=0;
        while(Math.pow(2, nextN)<Math.abs(A))nextN++;
        return nextN;
    }
    
    public static int nextpow2(int A){
        int nextN=0;
        while(Math.pow(2, nextN)<Math.abs(A))nextN++;
        return nextN;
    }
    
    public static double[] linspace(double min, double max, int points) {  
        double[] d = new double[points];  
        for (int i = 0; i < points; i++){  
            d[i] = min + i * (max - min) / (points - 1);  
        }  
        return d;  
    }
    
    public static double[] linspace(double min, int points, double dh) {  
        double[] d = new double[points];  
        for (int i = 0; i < points; i++){  
            d[i] = min + i * dh;  
        }  
        return d;  
    }
    
    public static double[] linspace(int points, double dh) {  
        double[] d = new double[points];  
        for (int i = 0; i < points; i++){  
            d[i] =  i * dh;  
        }  
        return d;  
    }
    
    public static double[] linspace(double min, double max) {  
        return linspace(min, max, 100);
    } 
    
    public static double random(double min, double max){
        return min +Math.random()*(max-min);
    }
    
    /**
     * Convolves sequences a and b.  The resulting convolution has
     * length a.length+b.length-1.
     */
    public static double[] conv(double[] a, double[] b)
    {
        double[] y = new double[a.length+b.length-1];

        // make sure that a is the shorter sequence
        if(a.length > b.length)
        {
            double[] tmp = a;
            a = b;
            b = tmp;
        }

        for(int lag = 0; lag < y.length; lag++)
        {
            y[lag] = 0;

            // where do the two signals overlap?
            int start = 0;
            // we can't go past the left end of (time reversed) a
            if(lag > a.length-1) 
                start = lag-a.length+1;

            int end = lag;
            // we can't go past the right end of b
            if(end > b.length-1)
                end = b.length-1;

            //System.out.println("lag = " + lag +": "+ start+" to " + end);
            for(int n = start; n <= end; n++)
            {
                //System.out.println("  ai = " + (lag-n) + ", bi = " + n); 
                y[lag] += b[n]*a[lag-n];
            }
        }

        return(y);
    }

    /**
     * Computes the cross correlation between sequences a and b.
     */
    public static double[] xcorr(double[] a, double[] b)
    {
        int len = a.length;
        if(b.length > a.length)
            len = b.length;

        return xcorr(a, b, len-1);

        // // reverse b in time
        // double[] brev = new double[b.length];
        // for(int x = 0; x < b.length; x++)
        //     brev[x] = b[b.length-x-1];
        // 
        // return conv(a, brev);
    }

    /**
     * Computes the auto correlation of a.
     */
    public static double[] xcorr(double[] a)
    {
        return xcorr(a, a);
    }

    /**
     * Computes the cross correlation between sequences a and b.
     * maxlag is the maximum lag to
     */
    public static double[] xcorr(double[] a, double[] b, int maxlag)
    {
        double[] y = new double[2*maxlag+1];
        Arrays.fill(y, 0);
        
        for(int lag = b.length-1, idx = maxlag-b.length+1; 
            lag > -a.length; lag--, idx++)
        {
            if(idx < 0)
                continue;
            
            if(idx >= y.length)
                break;

            // where do the two signals overlap?
            int start = 0;
            // we can't start past the left end of b
            if(lag < 0) 
            {
                //System.out.println("b");
                start = -lag;
            }

            int end = a.length-1;
            // we can't go past the right end of b
            if(end > b.length-lag-1)
            {
                end = b.length-lag-1;
                //System.out.println("a "+end);
            }

            //System.out.println("lag = " + lag +": "+ start+" to " + end+"   idx = "+idx);
            for(int n = start; n <= end; n++)
            {
                //System.out.println("  bi = " + (lag+n) + ", ai = " + n); 
                y[idx] += a[n]*b[lag+n];
            }
            //System.out.println(y[idx]);
        }

        return(y);
    }
    
    public static double[] pad(double[] f){
        int Nold=f.length;
        Double myDouble = Math.ceil(Math.pow(2, nextpow2(Nold)));
        int N = myDouble.intValue();
        double[] fn=new double[N];
        for(int i=0;i<Nold;i++){
            fn[i]=f[i];
        }
        for(int i=Nold;i<N;i++){
            fn[i]=0.0;
        }
        return fn;
    }
    
    public static double[] pad(ArrayList x){
        double[] f=new double[x.size()];
        for(int i=0;i<x.size();i++){
            f[i]=(double)x.get(i);
        }
        return pad(f);
    }
    
    public static Complex[] pad(Complex[] f){
        int Nold=f.length;
        Double myDouble = Math.ceil(Math.pow(2, nextpow2(Nold)));
        int N = myDouble.intValue();
        Complex[] fn=new Complex[N];
        for(int i=0;i<Nold;i++){
            fn[i]=new Complex(f[i].getReal(),f[i].getImaginary());
        }
        for(int i=Nold;i<N;i++){
            fn[i]=new Complex(0.0,0.0);
        }
        return fn;
    }
}
