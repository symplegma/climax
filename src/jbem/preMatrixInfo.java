/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class preMatrixInfo {
    private String fileName;
    private int m_minus_n;
    private int pos;
    
    // constructor
    public preMatrixInfo(int m_minus_n, int pos){
        this.m_minus_n=m_minus_n;
        this.pos=pos;
        fileName="mat"+m_minus_n+"_"+pos+".spm";
    }
    
    
    public int get_m_minus_n(){
        return this.m_minus_n;
    }
    
    public int get_pos(){
        return this.pos;
    }
    
    public String get_Name(){
        return this.fileName;
    }

}
