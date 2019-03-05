/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sympmesher;

import geom.Point;

/**
 *
 * @author pchr
 */
public class Patch {
    private final int id;
    private final Curve3P C1,C2,C3,C4;
    private boolean r1,r2,r3,r4;
    private double ratioXsi=1.0;
    private double ratioEta=1.0;
    private int meshXsi=1;
    private int meshEta=1;
    private Area2D theArea2D;
    
    public Patch(int id, Curve3P C1, Curve3P C2, Curve3P C3, Curve3P C4){
        this.id=id;
        this.C1=C1;
        this.C2=C2;
        this.C3=C3;
        this.C4=C4;
    }
    
    public Patch(int id, Curve3P C1, Curve3P C2, Curve3P C3, Curve3P C4, int meshx, int meshy){
        this.id=id;
        this.meshXsi=meshx;
        this.meshEta=meshy;
        this.C1=C1; this.C1.putPatch(this,1);
        this.C2=C2; this.C2.putPatch(this,2);
        this.C3=C3; this.C3.putPatch(this,3);
        this.C4=C4; this.C4.putPatch(this,4);
    }
    
    public Patch(int id, Curve3P C1, Curve3P C2, Curve3P C3, Curve3P C4, int meshx, int meshy, double rx, double ry){
        this.id=id;
        this.meshXsi=meshx;
        this.meshEta=meshy;
        this.ratioXsi=rx;
        this.ratioEta=ry;
        this.C1=C1; this.C1.putPatch(this,1);
        this.C2=C2; this.C2.putPatch(this,2);
        this.C3=C3; this.C3.putPatch(this,3);
        this.C4=C4; this.C4.putPatch(this,4);
    }
    
    public Patch(int id, Curve3P C1, Curve3P C2, Curve3P C3, Curve3P C4, double rx, double ry){
        this.id=id;
        this.ratioXsi=rx;
        this.ratioEta=ry;
        this.C1=C1; this.C1.putPatch(this,1);
        this.C2=C2; this.C2.putPatch(this,2);
        this.C3=C3; this.C3.putPatch(this,3);
        this.C4=C4; this.C4.putPatch(this,4);
    }
    
    private void setReversedCurve(int wC){
        switch(wC){
            case 1: r1=true; break;
            case 2: r2=true; break;
            case 3: r3=true; break;
            case 4: r4=true; break;
        }
    }
    
    private boolean getReversedCurve(int wC){
        boolean var=false;
        switch(wC){
            case 1: var=r1; break;
            case 2: var=r2; break;
            case 3: var=r3; break;
            case 4: var=r4; break;
        }
        return var;
    }
    
    public void setArea2D(Area2D area){this.theArea2D=area;}
    
    public double getratioXsi(){return this.ratioXsi;}
    
    public double getratioEta(){return this.ratioEta;}
    
    public int getmeshXsi(){return meshXsi;}
    
    public int getmeshEta(){return meshEta;}

    public int getID(){
        return this.id;
    }
    
    public void meshit(){
        int[][] nodes =new int[this.meshEta+1][this.meshXsi+1];
        // fo rthe geometric progression
        // https://en.wikipedia.org/wiki/Geometric_progression
        double alphaXsi = 0, alphaEta = 0;

        for(int i=1; i<=this.meshXsi; i++){
            alphaXsi+=Math.pow(this.ratioXsi, i-1);
        }
        alphaXsi=2./alphaXsi;

        for(int i=1; i<=this.meshEta; i++){
            alphaEta+=Math.pow(this.ratioEta, i-1);
        }
        alphaEta=2./alphaEta;
        
        double xsi; double eta=-1.0;
        for(int j=1; j<=meshEta+1; j++){
            xsi=-1.0;
            for(int i=1; i<=meshXsi+1; i++){
                double x,y,z;
                if(this.getReversedCurve(1)){
                    x=SF(2,xsi,eta)*C1.getPointStart().X()+SF(5,xsi,eta)*C1.getPointMid().X()+SF(1,xsi,eta)*C1.getPointEnd().X();
                    y=SF(2,xsi,eta)*C1.getPointStart().Y()+SF(5,xsi,eta)*C1.getPointMid().Y()+SF(1,xsi,eta)*C1.getPointEnd().Y();
                    z=SF(2,xsi,eta)*C1.getPointStart().Z()+SF(5,xsi,eta)*C1.getPointMid().Z()+SF(1,xsi,eta)*C1.getPointEnd().Z();
                }else{
                    x=SF(1,xsi,eta)*C1.getPointStart().X()+SF(5,xsi,eta)*C1.getPointMid().X()+SF(2,xsi,eta)*C1.getPointEnd().X();
                    y=SF(1,xsi,eta)*C1.getPointStart().Y()+SF(5,xsi,eta)*C1.getPointMid().Y()+SF(2,xsi,eta)*C1.getPointEnd().Y();
                    z=SF(1,xsi,eta)*C1.getPointStart().Z()+SF(5,xsi,eta)*C1.getPointMid().Z()+SF(2,xsi,eta)*C1.getPointEnd().Z();
                }
                if(this.getReversedCurve(2)){
                    x+=SF(6,xsi,eta)*C2.getPointMid().X()+SF(3,xsi,eta)*C2.getPointStart().X();
                    y+=SF(6,xsi,eta)*C2.getPointMid().Y()+SF(3,xsi,eta)*C2.getPointStart().Y();
                    z+=SF(6,xsi,eta)*C2.getPointMid().Z()+SF(3,xsi,eta)*C2.getPointStart().Z();
                }else{
                    x+=SF(6,xsi,eta)*C2.getPointMid().X()+SF(3,xsi,eta)*C2.getPointEnd().X();
                    y+=SF(6,xsi,eta)*C2.getPointMid().Y()+SF(3,xsi,eta)*C2.getPointEnd().Y();
                    z+=SF(6,xsi,eta)*C2.getPointMid().Z()+SF(3,xsi,eta)*C2.getPointEnd().Z();
                }
                if(this.getReversedCurve(3)){
                    x+=SF(7,xsi,eta)*C3.getPointMid().X()+SF(4,xsi,eta)*C3.getPointStart().X();
                    y+=SF(7,xsi,eta)*C3.getPointMid().Y()+SF(4,xsi,eta)*C3.getPointStart().Y();
                    z+=SF(7,xsi,eta)*C3.getPointMid().Z()+SF(4,xsi,eta)*C3.getPointStart().Z();
                }else{
                    x+=SF(7,xsi,eta)*C3.getPointMid().X()+SF(4,xsi,eta)*C3.getPointEnd().X();
                    y+=SF(7,xsi,eta)*C3.getPointMid().Y()+SF(4,xsi,eta)*C3.getPointEnd().Y();
                    z+=SF(7,xsi,eta)*C3.getPointMid().Z()+SF(4,xsi,eta)*C3.getPointEnd().Z();
                }
                x+=SF(8,xsi,eta)*C4.getPointMid().X();
                y+=SF(8,xsi,eta)*C4.getPointMid().Y();
                z+=SF(8,xsi,eta)*C4.getPointMid().Z();
                
                int existedID =this.theArea2D.ExistBoundaryLocation(x, y, z);
                if(existedID==0){
                    Point aP= new Point(x,y, z);
                    this.theArea2D.putPoint(aP);
                    if(j==1 || j==(meshEta+1) || i==1 || i==(meshXsi+1))theArea2D.putBoundaryPoint(aP);
                    nodes[j-1][i-1]=aP.getID();
                }else{
                    nodes[j-1][i-1]=existedID;
                }
                xsi+=alphaXsi*Math.pow(ratioXsi, i-1);
            }
            eta+=alphaEta*Math.pow(ratioEta, j-1);
        }
        for(int j=1; j<meshEta+1; j++){
            for(int i=1; i<meshXsi+1; i++){
                theArea2D.putquads(new quad(nodes[j-1][i-1],nodes[j-1][i],nodes[j][i],nodes[j][i-1]));
            }
        }
    }
    
    public void print(){
        System.out.println(this.id+": "+this.C1.getID()+","+this.C2.getID()+","+this.C3.getID()+","+this.C4.getID()+". (rx,ry)="+" ("+ratioXsi+","+ratioEta+")"+", (mx,my)="+" ("+meshXsi+","+meshEta+")");
    }
    
    private double SF(int wSF, double xsi, double eta){
        double N;
        switch(wSF){
            case 1: N=0.25*(1-xsi)*(1-eta)
                    -0.5*(0.5*(1-xsi*xsi)*(1-eta))
                    -0.5*(0.5*(1-xsi)*(1-eta*eta));
            break;
            case 2: N=0.25*(1+xsi)*(1-eta)
                    -0.5*(0.5*(1-xsi*xsi)*(1-eta))
                    -0.5*(0.5*(1+xsi)*(1-eta*eta));
            break;
            case 3: N=0.25*(1+xsi)*(1+eta)
                    -0.5*(0.5*(1+xsi)*(1-eta*eta))
                    -0.5*(0.5*(1-xsi*xsi)*(1+eta));
            break;
            case 4: N=0.25*(1-xsi)*(1+eta)
                    -0.5*(0.5*(1-xsi*xsi)*(1+eta))
                    -0.5*(0.5*(1-xsi)*(1-eta*eta));
            break;
            case 5: N=0.5*(1-xsi*xsi)*(1-eta);
            break;
            case 6: N=0.5*(1+xsi)*(1-eta*eta);
            break;
            case 7: N=0.5*(1-xsi*xsi)*(1+eta);
            break;
            case 8: N=0.5*(1-xsi)*(1-eta*eta);
            break;
            default:N=0. ;
            break;
        }
        return N;
    }
}
