/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jbem;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class InterfaceZeroDimension {
    protected int id;
    protected Map<Integer,Node> theNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,ENode> theElements = new TreeMap<Integer,ENode>();
    protected Material theMaterial;
    protected Spring NormalSpring;
    protected Spring TangentSpring;
    protected EnergyDissipator theDissipator;
    protected Map<Integer,InterfaceNode> theInterfaceNodes = new TreeMap<Integer,InterfaceNode>();
    protected Map<Integer,IENode> theInterfaceElements = new TreeMap<Integer,IENode>();
    protected Map<Integer,MinimizationConstraintEquation> theMinimizationEquations = new TreeMap<Integer,MinimizationConstraintEquation>();
    
    protected boolean zconstant = false;
    protected boolean sconstant = false;
    protected boolean pconstant = false;

    protected boolean zcontiuous = true;
    protected boolean scontiuous = true;
    protected boolean pcontiuous = true;

    protected boolean tangential = false;
    protected boolean normal = false;

    protected boolean damage = false;
    protected boolean slip = false;
    protected boolean plast = false;

    protected int sumDOFS;
    protected int utDOFS;
    protected int unDOFS;
    protected int zDOFS;
    protected int sDOFS;
    protected int pDOFS;
    private int beginDOF=0;
    
    public InterfaceZeroDimension(){}

    public InterfaceZeroDimension(int id){
        this.id=id;
    }
    
    public int getID(){
        return this.id;
    }

    public void addMaterial(Material aMaterial){
        this.theMaterial=aMaterial;
    }


    public void putNode(Node aNode){
        this.theNodes.put(aNode.getID(), aNode);
    }

    public void putElement(ENode anElement){
        this.theElements.put(anElement.getID(), anElement);
    }

    public Element getElement(int id){
        return this.theElements.get(id);
    }

    public Node getNode(int id){
        return this.theNodes.get(id);
    }

    public Map<Integer,ENode> getElements(){
        return this.theElements;
    }

    public Map<Integer,IENode> getInterfaceElements(){
        return this.theInterfaceElements;
    }

    public Map<Integer,Node> getNodes(){
        return this.theNodes;
    }

    public Map<Integer,InterfaceNode> getInterfaceNodes(){
        return this.theInterfaceNodes;
    }

    public void setSpaceIntegrator(SpaceIntegrator SI){
        for(Iterator<IENode> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            IENode theElement = it.next();
            theElement.setSpaceIntegrator(SI);
        }
    }
}
