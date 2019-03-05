/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
abstract public class TransformationPlane extends Transformation{
    
    public TransformationPlane(){}
    
    abstract public double getJacobian(double xsi, double eta);
    abstract protected double getShapeFunction(int wSF, double xsi, double eta);
    abstract protected double getShapeFunction_xsi(int wSF, double xsi, double eta);
    abstract protected double getShapeFunction_eta(int wSF, double xsi, double eta);
    abstract public double getXsi(double xsi_, double eta_);
    abstract public double getEta(double xsi_, double eta_);

}
