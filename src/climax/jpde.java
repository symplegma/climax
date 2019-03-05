/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package climax;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;

/**
 *
 * @author pchr
 */
public class jpde {
    public Mesh mesh;
    public Vector[] um;
    public int numSteps;
    
    public jpde(){}
    
    public jpde(Mesh m){this.mesh=m;}
    
    public void setMesh(Mesh m){this.mesh=m;}
    
    public void initSol(int steps){
        numSteps=steps;
        um=new Vector[steps];
    }

    public void setSol(Vector u, int atStep){
        um[atStep]=u;
    }

    public Vector getSol(int atStep){
        return this.um[atStep];
    }

    public int getNumSteps(){return this.numSteps;}
}
