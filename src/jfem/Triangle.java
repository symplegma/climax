/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jfem;

/**
 *
 * @author pchr
 */
abstract public class Triangle extends Element{
    
// default constructor
    public Triangle(){
        ndofs = 6;
        dimension=2;
    }
    
}
