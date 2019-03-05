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

import jmat.AbstractMatrix;
import static jfem.Quad.MyMatrix;

/**
 *
 * @author pchr
 */
public class Shell extends Element{
    private static AbstractMatrix MyMatrix = new AbstractMatrix(24,24);
    private static AbstractMatrix Bm = new AbstractMatrix(3,12,0.);
    private static AbstractMatrix Bb = new AbstractMatrix(3,12,0.);
    private static AbstractMatrix Bs = new AbstractMatrix(2,16,0.);
    private static AbstractMatrix K0 = new AbstractMatrix(12,12,0.);
    private static AbstractMatrix T = new AbstractMatrix(24,24,0.);
    private double[] detJ = new double[5];
    private double[][][] shp = new double[8][3][5];
    private final static int GaussPoints=2; // per direction of integration
    private static int numberOfShell = 0;
    
    public Shell(int id, Node node1, Node node2, Node node3, Node node4,
                             Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfShell;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        this.putNode(node4);
        int[] dofs = new int[6];
        dofs[0]=1; dofs[1]=2; ; dofs[2]=3; dofs[3]=4; dofs[4]=5; ; dofs[5]=6;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        node4.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(4, node4.getID());
        ndofs = 24;
        dimension=2;
    }
    
    public double ShapeFunction(int wsf, double xsi, double eta){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean ElementPlastified() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AbstractMatrix getK() {
        for(int i=0; i<ndofs; i++){
            for(int j=0; j<ndofs; j++){
                MyMatrix.putVal(i, j, 0.0);
            }
        }
        double E =((ElasticMaterial)this.ElementMaterial).getElasticity();
	double nu=((ElasticMaterial)this.ElementMaterial).getPoisson();
	double th=this.theCrossSection.getThickness();
        
        return MyMatrix;
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
    AbstractMatrix getK_v() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    AbstractMatrix getK_a() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public AbstractMatrix getT() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public double getVolume(){
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public double[] getNormal(int nid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
