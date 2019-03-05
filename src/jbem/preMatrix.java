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
public class preMatrix {
    private int m_minus_n;
    private int pos;
    private AbstractMatrix theMatrix;
    
    // constructor
    public preMatrix(int m_minus_n, int pos, AbstractMatrix theMatrix){        
        this.m_minus_n=m_minus_n;
        this.pos=pos;
        int row=theMatrix.getColumnDimension();
        int col=theMatrix.getRowDimension();
        this.theMatrix = new AbstractMatrix(col,row);
        this.theMatrix.init();
        this.theMatrix=this.theMatrix.plus(theMatrix);
    }
    
    // methods
    public int get_m_minus_n(){
        return this.m_minus_n;
    }
    
    public int get_pos(){
        return this.pos;
    }
    
    public AbstractMatrix get_Matrix(){
        return this.theMatrix;
    }

}
