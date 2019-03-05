/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
abstract public class AbstractSOE {
    
    
    abstract public AbstractMatrix getB();
    abstract public AbstractMatrix getA();
    abstract public AbstractMatrix getX();
    
    abstract public void setB(int row, double value);
    abstract public void setA(int row, int clmn, double value);
    abstract public void addA(int row, int clmn, double value);
    abstract public void zeroA();
    abstract public void addB(int row, double value);
    
    abstract public void init(int order);
    abstract public void init(int order, int steps);
    
    abstract public void solve();
    
    abstract public void zeroB();
}
