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

/**
 *
 * @author pchr
 */
public class CrossSection {
    protected static int numberOfCrossSection = 0;
    protected int id;
    protected double A=0.;
    protected double As2=0.;
    protected double As3=0.;
    protected double J1=0.;
    protected double J2=0.;
    protected double J3=0.;
    
    protected double thickness=0.;
    protected double angle=0.0;
    
    
    // constructor
    public CrossSection(){}
    
    public CrossSection(int id){
        this.id=id;
        ++numberOfCrossSection;
    }
    
    public CrossSection(int id, double width, double height ){
        this.id=id;
        ++numberOfCrossSection;
        double w=width;
	double h=height;
	A=w*h;
	As2=5./6.*A;
	As3=5./6.*A;
	J1=h*w*(w*w+h*h)/12.;
	J2=h*w*w*w/12.;
	J3=w*h*h*h/12.;
    }
    
    public CrossSection(int id, double A){
        this.id=id;
        ++numberOfCrossSection;
	this.A=A;
    }
    
    public CrossSection(int id, double A, double As2, double As3, double J1, double J2, double J3 ){
        this.id=id;
        ++numberOfCrossSection;
	this.A=A;
	this.As2=As2;
	this.As3=As3;
	this.J1=J1;
	this.J2=J2;
	this.J3=J3;
    }
    
    public CrossSection(int id, double A, double As2, double J3){
        this.id=id;
        ++numberOfCrossSection;
	this.A=A;
	this.As2=As2;
	this.J3=J3;
    }
    
    public CrossSection(int id, double[] X, double[] Y){
        /*
        The geometric property constant J1 is the torsional stiffness constant for cross section. 
        For circular sections J1 is equal to the polar moment of inertia (J2 + J3). 
        For non-circular sections the torsional stiffness constant is not equal to the polar moment of 
        inertia (see torsion of non-circular section in a solid mechanics reference). 
        If no value is entered for J1, JBEM will compute the torsional stiffness constant as J2 + J3 
        which is correct for circular sections but not correct for non-circular sections.
        */
        this.id=id;
        ++numberOfCrossSection;
        int n=X.length;
        this.A=0.;
        this.J2=0.0;
        this.J3=0.0;
        this.J1=0.0;
        for(int i=0;i<(n-1);i++){
            A+=X[i]*Y[i+1]-X[i+1]*Y[i];
//            J3+=(Y[i]*Y[i]+Y[i]*Y[i+1]+Y[i+1]*Y[i+1])*(X[i]*Y[i+1]-X[i+1]*Y[i]);
//            J2+=(X[i]*X[i]+X[i]*X[i+1]+X[i+1]*X[i+1])*(X[i]*Y[i+1]-X[i+1]*Y[i]);
            J2+=((Y[i]+Y[i+1])*(Y[i]+Y[i+1])-Y[i]*Y[i+1])*(X[i]*Y[i+1]-X[i+1]*Y[i]);
            J3+=((X[i]+X[i+1])*(X[i]+X[i+1])-X[i]*X[i+1])*(X[i]*Y[i+1]-X[i+1]*Y[i]);
        }
        A+=X[n-1]*Y[0]-X[0]*Y[n-1];
//        J3+=(Y[n-1]*Y[n-1]+Y[n-1]*Y[0]+Y[0]*Y[0])*(X[n-1]*Y[0]-X[0]*Y[n-1]);
//        J2+=(X[n-1]*X[n-1]+X[n-1]*X[0]+X[0]*X[0])*(X[n-1]*Y[0]-X[0]*Y[n-1]);
        J2+=((Y[n-1]+Y[0])*(Y[n-1]+Y[0])-Y[n-1]*Y[0])*(X[n-1]*Y[0]-X[0]*Y[n-1]);
        J3+=((X[n-1]+X[0])*(X[n-1]+X[0])-X[n-1]*X[0])*(X[n-1]*Y[0]-X[0]*Y[n-1]);
        A=A/2.;
        J2=J2/12.;
        J3=J3/12.;
        J2=J2-this.getA()*this.getXs(X,Y)*this.getXs(X,Y);
        J3=J3-this.getA()*this.getYs(X,Y)*this.getYs(X,Y);
        
        J1=J2+J3;
    }
    
    private double getSy(double[] X, double[] Y){
        int n=X.length;
        double Ys=0.;
        for(int i=0;i<(n-1);i++){
            Ys+=(X[i]*Y[i+1]-X[i+1]*Y[i])*(Y[i]+Y[i+1]);
        }
        Ys+=(X[n-1]*Y[0]-X[0]*Y[n-1])*(Y[n-1]+Y[0]);
        return Ys/6.;
    }
    
    private double getSx(double[] X, double[] Y){
        int n=X.length;
        double Xs=0.;
        for(int i=0;i<(n-1);i++){
            Xs+=(X[i]*Y[i+1]-X[i+1]*Y[i])*(X[i]+X[i+1]);
        }
        Xs+=(X[n-1]*Y[0]-X[0]*Y[n-1])*(X[n-1]+X[0]);
        return Xs/6.;
    }
    
    public double getYs(double[] X, double[] Y){
        double Ys=0;
        Ys=this.getSx(X,Y)/this.getA();
        return Ys;
    }
    
    
    
    public double getXs(double[] X, double[] Y){
        double Xs=0;
        Xs=this.getSy(X,Y)/this.getA();
        return Xs;
    }
    
    // methods
    public double getA() {
        return A;
    }
    
    public double getAs2() {
        return As2;
    }
    
    public double getAs3() {
        return As3;
    }
    
    public double getJ1() {
        return J1;
    }
    
    public double getJ2() {
        return J2;
    }
    
    public double getJ3() {
        return J3;
    }
    
    public int getID() {
        return id;
    }
    
    public double getThickness() {
        return this.thickness;
    }
    
    public void setThickness(double thikness) {
        this.thickness=thikness;
    }
    
    public void setAngle(double theta){this.angle=theta;}
    
    public double getAngle(){return this.angle;}
    
    public void setJ1(double J1value){
        // user defined torsional stiffness
        this.J1=J1value;
    }
    
    public double DA_height(){
        return 0.0;
    }
    
    public double DJ3_height(){
        return 0.0;
    }
}
