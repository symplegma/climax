/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mathman;

/**
 *
 * @author pchr
 */
public class ConstantDoubleFunction implements DoubleFunction{
    double constVal=0.0;
    
    public ConstantDoubleFunction(double constv){
        this.constVal=constv;
    }
    
    @Override
    public double run(double x) {
        return constVal;
    }

    @Override
    public double run(double x, double y) {
        return constVal;
    }
}
