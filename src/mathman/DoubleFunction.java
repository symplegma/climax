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
public interface DoubleFunction {
    // http://www.javaworld.com/article/2092260/java-se/java-programming-with-lambda-expressions.html
    public double run(double x);
    
    public double run(double x, double y);
    
}
