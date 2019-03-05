/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jfem;

import java.util.Iterator;
import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class TwoDAcousticQuad extends Quad{
    private static int numberOfTwoDAcousticQuad = 0;
    // constructor
    public TwoDAcousticQuad(){}
    
    public TwoDAcousticQuad(int id, Node node1, Node node2, Node node3, Node node4,
                             Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfTwoDAcousticQuad;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        this.putNode(node4);
        int[] dofs = new int[1];
        dofs[0]=7;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        node4.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(4, node4.getID());
        this.dof_per_node=1;
        ndofs = 4;
    }

    public TwoDAcousticQuad(int id, Node node1, Node node2, Node node3, Node node4){
        this.id=id;
        ++numberOfTwoDAcousticQuad;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        this.putNode(node4);
        int[] dofs = new int[1];
        dofs[0]=7;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        node4.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(4, node4.getID());
        this.dof_per_node=1;
        ndofs = 4;
    }

    @Override
    AbstractMatrix getEmat() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public AbstractMatrix getK() {
        double x1,x2,w1,w2;
        double h;
        MyMatrix = new AbstractMatrix(4,4);
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                MyMatrix.putVal(i, j, 0.0);
            }
        }
        double  cmat2=((AcousticMaterial)this.ElementMaterial).getMatConstant();cmat2*=cmat2;

        for(int i=1; i<=GaussPoints; i++){
            for(int j=1; j<=GaussPoints; j++){
                x1=this.getGaussCoord(i);
                x2=this.getGaussCoord(j);
                w1=this.getGaussWeight(i);
                w2=this.getGaussWeight(j);
                h=this.theCrossSection.getThickness();
                double detJ=getDetJ(x1, x2);
                double J11=this.getJ(x1, x2).get(0,0); double J12=this.getJ(x1, x2).get(0,1);
                double J21=this.getJ(x1, x2).get(1,0); double J22=this.getJ(x1, x2).get(1,1);
                detJ=J11*J22-J12*J21;
                double N1x=ShapeFunction_xsi(1,x1,x2); double N1y=ShapeFunction_eta(1,x1,x2);
                double N2x=ShapeFunction_xsi(2,x1,x2); double N2y=ShapeFunction_eta(2,x1,x2);
                double N3x=ShapeFunction_xsi(3,x1,x2); double N3y=ShapeFunction_eta(3,x1,x2);
                double N4x=ShapeFunction_xsi(4,x1,x2); double N4y=ShapeFunction_eta(4,x1,x2);
                double coef=cmat2*h*detJ*w1*w2;
                
                MyMatrix.set(0, 0, MyMatrix.get(0, 0)+coef*((N1x*J22/detJ-N1y*J12/detJ)*(N1x*J22/detJ-N1y*J12/detJ)+(N1y*J11/detJ-N1x*J21/detJ)*(N1y*J11/detJ-N1x*J21/detJ)));
                MyMatrix.set(0, 1, MyMatrix.get(0, 1)+coef*((N1x*J22/detJ-N1y*J12/detJ)*(N2x*J22/detJ-N2y*J12/detJ)+(N1y*J11/detJ-N1x*J21/detJ)*(N2y*J11/detJ-N2x*J21/detJ)));
                MyMatrix.set(0, 2, MyMatrix.get(0, 2)+coef*((N1x*J22/detJ-N1y*J12/detJ)*(N3x*J22/detJ-N3y*J12/detJ)+(N1y*J11/detJ-N1x*J21/detJ)*(N3y*J11/detJ-N3x*J21/detJ)));
                MyMatrix.set(0, 3, MyMatrix.get(0, 3)+coef*((N1x*J22/detJ-N1y*J12/detJ)*(N4x*J22/detJ-N4y*J12/detJ)+(N1y*J11/detJ-N1x*J21/detJ)*(N4y*J11/detJ-N4x*J21/detJ)));
                
                MyMatrix.set(1, 1, MyMatrix.get(1, 1)+coef*((N2x*J22/detJ-N2y*J12/detJ)*(N2x*J22/detJ-N2y*J12/detJ)+(N2y*J11/detJ-N2x*J21/detJ)*(N2y*J11/detJ-N2x*J21/detJ)));
                MyMatrix.set(1, 2, MyMatrix.get(1, 2)+coef*((N2x*J22/detJ-N2y*J12/detJ)*(N3x*J22/detJ-N3y*J12/detJ)+(N2y*J11/detJ-N2x*J21/detJ)*(N3y*J11/detJ-N3x*J21/detJ)));
                MyMatrix.set(1, 3, MyMatrix.get(1, 3)+coef*((N2x*J22/detJ-N2y*J12/detJ)*(N4x*J22/detJ-N4y*J12/detJ)+(N2y*J11/detJ-N2x*J21/detJ)*(N4y*J11/detJ-N4x*J21/detJ)));
                
                MyMatrix.set(2, 2, MyMatrix.get(2, 2)+coef*((N3x*J22/detJ-N3y*J12/detJ)*(N3x*J22/detJ-N3y*J12/detJ)+(N3y*J11/detJ-N3x*J21/detJ)*(N3y*J11/detJ-N3x*J21/detJ)));
                MyMatrix.set(2, 3, MyMatrix.get(2, 3)+coef*((N3x*J22/detJ-N3y*J12/detJ)*(N4x*J22/detJ-N4y*J12/detJ)+(N3y*J11/detJ-N3x*J21/detJ)*(N4y*J11/detJ-N4x*J21/detJ)));
                
                MyMatrix.set(3, 3, MyMatrix.get(3, 3)+coef*((N4x*J22/detJ-N4y*J12/detJ)*(N4x*J22/detJ-N4y*J12/detJ)+(N4y*J11/detJ-N4x*J21/detJ)*(N4y*J11/detJ-N4x*J21/detJ)));
            }
        }
        MyMatrix.set(1,0,MyMatrix.get(0, 1));
        MyMatrix.set(2,0,MyMatrix.get(0, 2)); MyMatrix.set(2,1,MyMatrix.get(1, 2));
        MyMatrix.set(3,0,MyMatrix.get(0, 3)); MyMatrix.set(3,1,MyMatrix.get(1, 3)); MyMatrix.set(3,2,MyMatrix.get(2, 3));
        return MyMatrix;
    }

    @Override
    public int[] getFtable(){
        int[] theEFTable = new int[ndofs];
        int j=0;
        for(Iterator<Node> it=this.ElementNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theEFTable[j]=theNode.getFtable()[6];
            j++;
        }
        return theEFTable;
    }
    
    @Override
    protected AbstractMatrix getM(int consistent) {
        MyMatrix = new AbstractMatrix(ndofs,ndofs);
        double h=this.theCrossSection.getThickness();
        MyMatrix.init();
        switch(consistent){
            case 0:
                // that is for consistent
                for(int i=0; i<ndofs; i++){
                    for(int j=0; j<ndofs; j++){
                        MyMatrix.putVal(i, j, 0.0);
                    }
                }
                for(int i=1; i<=GaussPoints; i++){
                    for(int j=1; j<=GaussPoints; j++){
                        double x1=this.getGaussCoord(i);
                        double x2=this.getGaussCoord(j);
                        double w1=this.getGaussWeight(i);
                        double w2=this.getGaussWeight(j);
                        double detJ=getDetJ(x1, x2);
                        double N1=ShapeFunction(1,x1,x2); 
                        double N2=ShapeFunction(2,x1,x2); 
                        double N3=ShapeFunction(3,x1,x2); 
                        double N4=ShapeFunction(4,x1,x2);
                        double coef=h*detJ*w1*w2;
                        MyMatrix.set(0, 0, MyMatrix.get(0, 0)+coef*(N1*N1));
                        MyMatrix.set(0, 1, MyMatrix.get(0, 1)+coef*(N1*N2));
                        MyMatrix.set(0, 2, MyMatrix.get(0, 2)+coef*(N1*N3));
                        MyMatrix.set(0, 3, MyMatrix.get(0, 3)+coef*(N1*N4));
                        
                        MyMatrix.set(1, 1, MyMatrix.get(1, 1)+coef*(N2*N2));
                        MyMatrix.set(1, 2, MyMatrix.get(1, 2)+coef*(N2*N3));
                        MyMatrix.set(1, 3, MyMatrix.get(1, 3)+coef*(N2*N4));
                        
                        MyMatrix.set(2, 2, MyMatrix.get(2, 2)+coef*(N3*N3));
                        MyMatrix.set(2, 3, MyMatrix.get(2, 3)+coef*(N3*N4));
                        
                        MyMatrix.set(3, 3, MyMatrix.get(3, 3)+coef*(N4*N4));
                    }
                }
                MyMatrix.set(1,0,MyMatrix.get(0, 1));
                MyMatrix.set(2,0,MyMatrix.get(0, 2)); MyMatrix.set(2,1,MyMatrix.get(1, 2));
                MyMatrix.set(3,0,MyMatrix.get(0, 3)); MyMatrix.set(3,1,MyMatrix.get(1, 3)); MyMatrix.set(3,2,MyMatrix.get(2, 3));
                break;
            default:
                double A=getArea();
                for (int i = 0; i < ndofs; i++) {
                    MyMatrix.putVal(i, i,A*h/4.0);
                }
                break;
        }
//        MyMatrix=this.getT().transpose().times(MyMatrix);
//        return MyMatrix.times(this.getT());
        return MyMatrix;
    }
}
