/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gendomain;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import mater.Material;

/**
 *
 * @author pchr
 */
public class Element {
    protected int id;
    private static int numberOfElements = 0;
    //protected Vector ElementNodes;
    protected Map<Integer,Node> ElementNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,Integer> theNodesHierarchy = new TreeMap<Integer,Integer>();
    protected int numNodes=0;
    protected Material ElementMaterial;
    
    // constructors
    public Element(){}
        
    // methods
    public int getnumberOfElements() {
        return numberOfElements;
    }
    
    public int getID() {
        return id;
    }
    
    protected void putNode(Node aNode) {
        ElementNodes.put(++numNodes, aNode);
        aNode.putElement(this);
    }

    public int getNumNodes(){
        return this.ElementNodes.size();
    }

    public void print(){
        System.out.print("elem "+this.id+": ");
        for(Iterator<Node>it=this.ElementNodes.values().iterator();it.hasNext();){
            System.out.print(it.next().getID()+" ");
        }
//        System.out.println(((EBeam2d)this).getL());
        System.out.println(this.ElementMaterial.getID());
    }

    public Node getNodeHierarchy(int nodeHier) {
        return this.ElementNodes.get(nodeHier);
    }
    
    public int getHierOfNode(int theNode_id){
        int seq=0;
        int i=0;
        for(Iterator<Integer> it=this.theNodesHierarchy.values().iterator(); it.hasNext();){
            ++i;
            if(it.next()==theNode_id){seq=i;}
        }
        return seq;
    }

    public void setMaterial(Material Mater){
        this.ElementMaterial = Mater;
    }
    
     public Material getMaterial(){return this.ElementMaterial;}
}
