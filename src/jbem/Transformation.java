/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
abstract public class Transformation {
    protected int[] theNodes;
    
    public Transformation(){}
    
    public void setNodes(int[] theNodes){
        this.theNodes=theNodes;
    }

}
