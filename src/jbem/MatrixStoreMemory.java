/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author pchr
 */
public class MatrixStoreMemory extends MatrixStore{
    private Map<Integer,preMatrix> theMatrices = new HashMap<Integer,preMatrix>();
    
    // constructor
    public MatrixStoreMemory(){}
    
    // methods
    public void putMatrix(int m,int n, int pos,AbstractMatrix theMatrix){
        ++numOfMatrices;
        preMatrix apreMatrix = new preMatrix(m-n,pos,theMatrix);
        this.theMatrices.put(numOfMatrices, apreMatrix);
    }
    
    public boolean isStored(int m, int n, int pos){
        boolean check=false;
        int mmn=m-n;
        for(Iterator<preMatrix>it=this.theMatrices.values().iterator();it.hasNext();){
            preMatrix apreMatrix=it.next();
            int the_m=apreMatrix.get_m_minus_n();
            int the_pos=apreMatrix.get_pos();
            if( (the_m==mmn)&&(the_pos==pos) ){
                check=true;
            }
            if(check!=true)break;
        }
        return check;
    }
    
    public AbstractMatrix getMatrix(int m, int n, int pos){
        boolean check=false;
        AbstractMatrix Mat = null;
        int mmn=m-n;
        for(Iterator<preMatrix>it=this.theMatrices.values().iterator();it.hasNext();){
            preMatrix apreMatrix=it.next();
            int the_m=apreMatrix.get_m_minus_n();
            int the_pos=apreMatrix.get_pos();
            if( (the_m==mmn)&&(the_pos==pos) ){
                Mat=apreMatrix.get_Matrix();
                check=true;
            }
            if(check!=true)break;
        }
        return Mat;
    }

}
