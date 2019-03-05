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
public class ZeroDoubleFunction implements DoubleFunction{

    @Override
    public double run(double x) {
        return 0.0;
    }

    @Override
    public double run(double x, double y) {
        return 0.0;
    }
    
}
