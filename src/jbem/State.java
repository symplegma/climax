/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

/**
 *
 * @author pchr
 */
public class State {
    private int id;
    private String name;
    private int responseID;
    private double eigval;
    
    public State(int id, String name){
        this.id=id;
        this.name=name;
    }
    
    public int getID(){return this.id;}
    
    public String getName(){return this.name;}
    
    public void setResponseID(int responseID){this.responseID=responseID;}
    
    public int getResponseID(){return this.responseID;}
    
    public void setEigenValue(int responseID){this.responseID=responseID;}
    
    public int getEigenValue(){return this.responseID;}
}
