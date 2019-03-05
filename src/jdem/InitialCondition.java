package jdem;

import gendomain.Node;

/**
 *
 * @author pchr
 */
public class InitialCondition {
    protected int id;
    protected static int numberOfConstraints = 0;
    protected double InitialValue=0.;
    protected Node wnode;
    protected int wdof; // local dof of wnode
    protected int type; // 0 fo disp, 1 for vel
    
    // constructor    
    public InitialCondition(double InitialValue, Node wnode, int wdof){
        this.id=++numberOfConstraints;
        this.InitialValue=InitialValue;
        this.wnode=wnode;
        this.wdof=wdof;
    }
   
    
    // methods
    public double getInitialValue(){
        return InitialValue;
    }
    
    public Node getNode(){
        return wnode;
    }
    
    public int getDOF(){
        return wdof;
    }
    
    public int getID(){
        return id;
    }
    
    public void setType(int t){
        this.type=t;
    }
    
    public int getType(){
        return type;
    }
}
