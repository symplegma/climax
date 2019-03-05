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

/**
 *
 * @author pchr
 */
public class PlaneStressQuad extends Quad{
    private static int numberOfPlaneStressQuad = 0;
    // constructor
    public PlaneStressQuad(int id, Node node1, Node node2, Node node3, Node node4,
                             Material Mater, CrossSection Sect){
        this.id=id;
        ++numberOfPlaneStressQuad;
        this.ElementMaterial = Mater;
        this.theCrossSection = Sect;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        this.putNode(node4);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        node4.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(4, node4.getID());
        this.dof_per_node=2;
        ndofs = 8;
    }
    
    public PlaneStressQuad(int id, Node node1, Node node2, Node node3, Node node4){
        this.id=id;
        ++numberOfPlaneStressQuad;
        this.putNode(node1);
        this.putNode(node2);
        this.putNode(node3);
        this.putNode(node4);
        int[] dofs = new int[2];
        dofs[0]=1; dofs[1]=2;
        node1.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(1, node1.getID());
        node2.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(2, node2.getID());
        node3.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(3, node3.getID());
        node4.setNdofs_ofNode(dofs);this.theNodesHierarchy.put(4, node4.getID());
        this.dof_per_node=2;
    }
    
    public AbstractMatrix getEmat(){
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                Matrix_3_3.putVal(i, j, 0.0);
            }
            double E=((ElasticMaterial) this.ElementMaterial).getElasticity();
            double v=((ElasticMaterial) this.ElementMaterial).getPoisson();
            double coef=E/(1.-v*v);
            Matrix_3_3.putVal(0, 0,  coef );
            Matrix_3_3.putVal(0, 1,  coef*v);
            Matrix_3_3.putVal(1, 0,  coef*v);
            Matrix_3_3.putVal(1, 1,  coef);
            Matrix_3_3.putVal(2, 2,  0.5*(1.-v)*coef );
        }
        return Quad.Matrix_3_3;
    }
    
    @Override
    public double getStressVonMisses(double[] at, int LC,int step){
        double[] Stress = this.getStress(at, LC, step);
        return Math.sqrt(0.5*(
                (Stress[0]-Stress[1])*(Stress[0]-Stress[1])
                +Stress[1]*Stress[1]+Stress[0]*Stress[0]
                +6.0*Stress[2]*Stress[2])
        );
    }
}
