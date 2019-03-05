/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class TransientHeatTransfer2D extends TimeFundamentalSolution{
    
    // constructor
    public TransientHeatTransfer2D(){
        this.ndofs=1;
        RespectiveStatic = new SteadyStatePotential2DFS();
        SpaceDimension=2;
    }

    @Override
    public AbstractMatrix get_pdif_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_v_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int get_u_DOFs() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int get_v_DOFs() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int get_p_DOFs() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_u_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_p_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_s_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AbstractMatrix get_r_fund() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
