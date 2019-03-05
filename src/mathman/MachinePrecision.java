/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package mathman;

/**
 *
 * @author pchr
 */
public final class MachinePrecision {
    
    public static double getMachinePrecision(){
        double m_Epsilon=1.0;
	while(1.0+m_Epsilon > 1.0){
	    m_Epsilon /= 2.0;	    
	}
	m_Epsilon *= 2.0;
        return m_Epsilon;
    }
    
    public static double getMachineZero(){
        double m_Zero;
	m_Zero = Math.sqrt(getMachinePrecision());
        return m_Zero;
    }
}
