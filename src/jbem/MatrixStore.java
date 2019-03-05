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
abstract public class MatrixStore {
    protected static int numOfMatrices=0;
    
    // constructor
    public MatrixStore(){}
    
    // methods
    abstract public void putMatrix(int m,int n, int pos,AbstractMatrix theMatrix);
    
    abstract public boolean isStored(int m, int n, int pos);
    
    abstract public AbstractMatrix getMatrix(int m, int n, int pos);

}
