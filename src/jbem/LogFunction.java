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
public interface LogFunction {
    AbstractMatrix getuLogPart();
    AbstractMatrix getuNonLogPart(double val);
    AbstractMatrix getvLogPart();
    AbstractMatrix getvNonLogPart();
}
