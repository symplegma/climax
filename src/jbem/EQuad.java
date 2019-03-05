/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import java.util.Iterator;

/**
 *
 * @author pchr
 */
abstract public class EQuad extends EPlane{
    
    // constructor
    public EQuad(){};
    
    public double[] getNormal(int nid){
        double xsi = 0., eta=0.;
        switch(this.getHierOfNode(nid)){
            case 1: xsi= -1.; eta= -1.; break;
            case 2: xsi=  1.; eta= -1.; break;
            case 3: xsi=  1.; eta=  1.; break;
            case 4: xsi= -1.; eta=  1.; break;
            case 5: xsi=  0.; eta= -1.; break;
            case 6: xsi=  1.; eta=  0.; break;
            case 7: xsi=  0.; eta=  1.; break;
            case 8: xsi= -1.; eta=  0.; break;
            case 9: xsi=  0.; eta=  0.; break;
            default: System.err.println("does not exist node with id= "+nid+" on element with id= "+id); System.exit(nid); break;
        }
        return this.getNormal(xsi, eta);
    }
    
    public double getInternalPointDisp(Domain theDomain, ResultPoint theInternalPoint, int wdisp,int step,int wstate){
        double val=0.;
        SpaceQuadIntegrator quadSI = (SpaceQuadIntegrator) this.theSIntegrator;
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
            val+=quadSI.IntegrateInternalPointDisp(this,theNode, theDomain, theInternalPoint, wdisp, step, wstate);
        }
        return val;
    }
    
    @Override
    public double getInternalPointStress(Domain theDomain, ResultPoint theInternalPoint, int ws1,int step,int wstate){
        double val=0.;
        SpaceQuadIntegrator quadSI = (SpaceQuadIntegrator) this.theSIntegrator;
        for(Iterator<Node> nt=this.getNodes().values().iterator(); nt.hasNext();){
            Node theNode = nt.next();
            val+=quadSI.IntegrateInternalPointStress(this,theNode, theDomain, theInternalPoint, ws1, step, wstate);
        }
        return val;
    }
}
