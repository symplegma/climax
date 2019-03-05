/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gendomain;

/**
 *
 * @author pchr
 */
public class DOFGroup {
    private int[] FTable;
    private DOFDescription dofDescription;
    
    public enum DOFDescription {
        DISPLACEMENTS, STRESS, STRAIN, VELOCITY
    }
    
    public DOFGroup(int n, DOFDescription theDescription){
        FTable = new int[n];
        for(int i=0;i<n;i++){
            FTable[i]=0;
        }
        dofDescription=theDescription;
    }
}
