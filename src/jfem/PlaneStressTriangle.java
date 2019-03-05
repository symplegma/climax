/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jfem;

import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class PlaneStressTriangle extends Triangle{    
    private static int numberOfPlaneStressTriangle = 0;
    
    // constructor
    public PlaneStressTriangle(int id, Node node1, Node node2, Node node3,
                             Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfPlaneStressTriangle;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        this.dof_per_node=2;
        ndofs = 6;
    }
    
    public PlaneStressTriangle(int id, Node node1, Node node2, Node node3){
        this.id=id;
        ++numberOfPlaneStressTriangle;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        this.dof_per_node=2;
    }

    @Override
    public boolean ElementPlastified() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix getK() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix getM() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    int getElementNdofs() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    int[] getFtable() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix getF() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    void clear() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getM_v() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getM_a() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getM_u() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getK_u() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getK_v() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getK_a() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    void commit() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getuKu() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getvMv() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getBVNorm(double coef) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix getT() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getVolume() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
