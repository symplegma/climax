/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

/**
 *
 * @author pchr
 */
public class MinDOF {
    private int nodeID;
    private int dof_sequence;
    private int DomainID;
    private int mdof;
    private static int numDOFS=0;

    public MinDOF(int DomainID, int nodeID, int dof_sequence){
        ++ numDOFS;
        this.DomainID=DomainID;
        this.mdof = numDOFS;
        this.nodeID=nodeID;
        this.dof_sequence=dof_sequence;
    }

    public MinDOF(int DomainID, int nodeID){
        ++ numDOFS;
        this.DomainID=DomainID;
        this.mdof = numDOFS;
        this.nodeID=nodeID;
    }


    public int getMDOF(){return this.mdof;}

    public int getNodeID(){return this.nodeID;}

    public int getDomainID(){return this.DomainID;}

    public int getDOF_sequence(){return this.dof_sequence;}

    public static int getNUMofDOFS(){return numDOFS;}

}
