/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sympmesher;

/**
 *
 * @author pchr
 */
public class quad{
    int n1,n2,n3,n4;
    public quad(int n1, int n2, int n3, int n4){
        this.n1=n1; this.n2=n2; this.n3=n3; this.n4=n4;
    }
    public int getN1(){return n1;}
    public int getN2(){return n2;}
    public int getN3(){return n3;}
    public int getN4(){return n4;}
}
