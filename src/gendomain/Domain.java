/*******************************************************************************
* Climax.                                                                      *
* Copyright (C) 2009-2017 C.G. Panagiotopoulos [http://www.symplegma.org]      *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): C.G. Panagiotopoulos (pchr76@gmail.com)
// *****************************************************************************
package gendomain;

import delaunay.Pnt;
import delaunay.Triangulation;
import geom.Point;
import geom.Shape;
import geom.Triangle;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import mater.ElasticMat;
import mater.Material;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;

/**
 *
 * @author pchr
 */
public class Domain {
    private int id;
    private static int numberOfDomains = 0;
    protected Map<Integer,Node> theNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,Element> theElements = new TreeMap<Integer,Element>();
    protected Map<Integer,Material> theMaterials = new TreeMap<Integer,Material>();
    protected Map<Integer,Shape> theShapes = new TreeMap<Integer,Shape>();
    protected Map<Integer,ConstraintEquation> theConstraintEquations = new TreeMap<Integer,ConstraintEquation>();
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Z_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Z_coordinate=Double.POSITIVE_INFINITY;
    private double Lx;
    private double Ly;
    private int Nx;
    private int Ny;
    private int dofcounter=0;
    private boolean isTriangulized=false;
    public Color ColorNodesUND=Color.green;
    public Color ColorNodesDEF=Color.yellow;
    public boolean fillNodes=true;
    
    public Domain(){
        id = ++numberOfDomains;
    }
    
    public void clsNumberOfDomains(){numberOfDomains=0;}
    
    public Domain(double Lx, double Ly, int Nx, int Ny){
        id = ++numberOfDomains;
        this.Lx=Lx;
        this.Ly=Ly;
        this.Nx=Nx;
        this.Ny=Ny;
    }
    
    public Domain(double L, int N){
        id = ++numberOfDomains;
        this.Lx=L;
        this.Ly=L;
        this.Nx=N;
        this.Ny=N;
    }
    
    public Domain(double Lx, double Ly, int N){
        id = ++numberOfDomains;
        this.Lx=Lx;
        this.Ly=Ly;
        this.Nx=N;
        this.Ny=N;
    }
    
    public Domain(double L, int Nx, int Ny){
        id = ++numberOfDomains;
        this.Lx=L;
        this.Ly=L;
        this.Nx=Nx;
        this.Ny=Ny;
    }
    
    public int getID(){
        return this.id;
    }
    
    public void setLxy(double Lx, double Ly){this.Lx=Lx; this.Ly=Ly;}
    public void setNxy(int Nx, int Ny){this.Nx=Nx; this.Ny=Ny;}
    
    public double getLx(){return this.Lx;}
    
    public double getLy(){return this.Ly;}
    
    public void putNode(Node aNode){
        aNode.setuEFTable(dofcounter+1);
        dofcounter+=aNode.getNumDOFS();
        if(aNode.getCoordinates()[0]>maximum_X_coordinate)maximum_X_coordinate=aNode.getCoordinates()[0];
        if(aNode.getCoordinates()[0]<minimum_X_coordinate)minimum_X_coordinate=aNode.getCoordinates()[0];
        if(aNode.getCoordinates().length>1){
            if(aNode.getCoordinates()[1]>maximum_Y_coordinate)maximum_Y_coordinate=aNode.getCoordinates()[1];
            if(aNode.getCoordinates()[1]<minimum_Y_coordinate)minimum_Y_coordinate=aNode.getCoordinates()[1];
        }
        if(aNode.getCoordinates().length>2){
            if(aNode.getCoordinates()[2]>maximum_Z_coordinate)maximum_Z_coordinate=aNode.getCoordinates()[2];
            if(aNode.getCoordinates()[2]<minimum_Z_coordinate)minimum_Z_coordinate=aNode.getCoordinates()[2];
        }
        this.theNodes.put(aNode.getID(), aNode);
    }
    
    public Node getNode(int id){
        return this.theNodes.get(id);
    }

    public Map<Integer,Node> getNodes(){
        return this.theNodes;
    }
    
    public void putElement(Element aElement){
        this.theElements.put(aElement.getID(), aElement);
    }
    
    public Element getElement(int id){
        return this.theElements.get(id);
    }

    public Map<Integer,Element> getElements(){
        return this.theElements;
    }
    
    public Map<Integer,Shape> getShapes(){
        return this.theShapes;
    }
    
    public void setShapes(Map<Integer,Shape> someShapes){
        this.theShapes=someShapes;
        this.isTriangulized=true;
    }
    
    public void putMaterial(Material aMaterial){
        this.theMaterials.put(aMaterial.getID(), aMaterial);
    }
    
    public Material getMaterial(int id){
        return this.theMaterials.get(id);
    }

    public Map<Integer,Material> getMaterials(){
        return this.theMaterials;
    }
    
    public double getVp_max(){
        double Vp_max=Double.NEGATIVE_INFINITY;
        for (Element theElement : this.theElements.values()) {
            if(Vp_max<((ElasticMat)theElement.getMaterial()).getVp())Vp_max=((ElasticMat)theElement.getMaterial()).getVp();
        }
        return Vp_max;
    }
    
    public double getVs_max(){
        double Vs_max=Double.NEGATIVE_INFINITY;
        for (Element theElement : this.theElements.values()) {
            if(Vs_max<((ElasticMat)theElement.getMaterial()).getVs())Vs_max=((ElasticMat)theElement.getMaterial()).getVs();
        }
        return Vs_max;
    }
    
    public double getVp_min(){
        double Vp_min=Double.NEGATIVE_INFINITY;
        for (Element theElement : this.theElements.values()) {
            if(Vp_min>((ElasticMat)theElement.getMaterial()).getVp())Vp_min=((ElasticMat)theElement.getMaterial()).getVp();
        }
        return Vp_min;
    }
    
    public double getVs_min(){
        double Vs_min=Double.NEGATIVE_INFINITY;
        for (Element theElement : this.theElements.values()) {
            if(Vs_min>((ElasticMat)theElement.getMaterial()).getVs())Vs_min=((ElasticMat)theElement.getMaterial()).getVs();
        }
        return Vs_min;
    }
    
    public int find(double x){
        return find(x, 0.0);
    }
    
    public int find(double x, double y){
        int nodeID= -1;
        double dist=Double.POSITIVE_INFINITY;
        for (Node theNode : this.theNodes.values()) {
            if(dist>theNode.getDist(new Point(x,y,0.0))){
                dist=theNode.getDist(new Point(x,y,0.0));
                nodeID=theNode.getID();
            }
        }
        return nodeID;
    }
    
    public int find(double[] array){
        double z = 0.0,y = 0.0,x=array[0];
        if(array.length>1)y=array[1];
        if(array.length>2)z=array[2];
        int nodeID= -1;
        double dist=Double.POSITIVE_INFINITY;
        for (Node theNode : this.theNodes.values()) {
            if(dist>theNode.getDist(new Point(x,y,z))){
                dist=theNode.getDist(new Point(x,y,z));
                nodeID=theNode.getID();
            }
        }
        return nodeID;
    }
    
    public void printNodesGeom(){
        System.out.println("Nodes of domain with id= "+this.id);
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            aNode.print();
        }
    }
   
    public double maximum_X_coordinate(){return this.maximum_X_coordinate;}
    
    public double maximum_Y_coordinate(){return this.maximum_Y_coordinate;}
    
    public double maximum_Z_coordinate(){return this.maximum_Z_coordinate;}
    
    public double minimum_X_coordinate(){return this.minimum_X_coordinate;}
    
    public double minimum_Y_coordinate(){return this.minimum_Y_coordinate;}
    
    public double minimum_Z_coordinate(){return this.minimum_Z_coordinate;}
    
    public int getNumNodes(){return this.theNodes.size();}
    
    public int getNumElements(){return this.theElements.size();}
    
    public int getNx(){return this.Nx;}
    
    public int getNy(){return this.Ny;}
    
    public void updateNodes(double[] vals, int step){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            double[] nvals=new double[aNode.getNumDOFS()]; 
            for(int i=0;i<nvals.length;i++){
                nvals[i]=vals[(aNode.getID()-1)*nvals.length+i]; 
            }
            aNode.updateVals(nvals, step);
        }
    }
    
    public void updateNodes(double[] X){updateNodes(X, 0);}
    
    public void updateNodes(Array2DRowRealMatrix vals, int step){
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            double[] nvals=new double[aNode.getNumDOFS()]; 
            for(int i=0;i<nvals.length;i++){
                nvals[i]=vals.getEntry((aNode.getID()-1)*nvals.length+i, 0); 
            }
            aNode.updateVals(nvals, step);
        }
    }
    
    public void updateNodes(Array2DRowRealMatrix vals){
        updateNodes(vals, 0);
    }
    
    public void updateNodesDOF(double[] vals, int step, int dof){
        // dof is 1 or 2 or N <= ndofs_ofNode
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            aNode.updateValDOF(vals[(aNode.getID()-1)], step, dof);
        }
    }
    
    public void Triangulate(){
        //delaunay.Triangle.moreInfo = true;
        double mx,my,m;
        mx=Math.max(Math.abs(minimum_X_coordinate()), maximum_X_coordinate());
        my=Math.max(Math.abs(minimum_X_coordinate()), maximum_X_coordinate());
        m=Math.max(mx,my);
        delaunay.Triangle tri =
            new delaunay.Triangle(   new Pnt(-m*10000.0,m*10000.0), 
                            new Pnt(m*10000.0,m*10000.0), 
                            new Pnt(0,-m*10000.0));
        //System.out.println(tri);
        Triangulation dt = new Triangulation(tri);
        boolean repeat=false;
//        System.out.println("DelaunayTriangulation created: " + dt);
        for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
            Node aNode=it.next();
            Pnt aPnt = new Pnt(aNode.getCoords()[0],aNode.getCoords()[1]);
            aPnt.setIndex(aNode.getID());
//            System.out.println("Point id = "+aNode.getID());
            try{
                dt.delaunayPlace(aPnt);
            }catch(java.util.NoSuchElementException e){
                repeat=true;
            }
        }
        int indtr=0;indtr=10; 
        while(repeat && indtr<3){ // never used since indtr=10>3
            indtr++;
            repeat=false;
            System.out.println("Delaunay triangulation to be repeated: "+indtr);
            for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
                Node aNode=it.next();
                Pnt aPnt = new Pnt(aNode.getCoords()[0],aNode.getCoords()[1]);
                aPnt.setIndex(aNode.getID());
    //            System.out.println("Point id = "+aNode.getID());
                try{
                    dt.delaunayPlace(aPnt);
                }catch(java.util.NoSuchElementException e){
                    repeat=true;
                }
            }
        }
        indtr=0;
        for (delaunay.Triangle triangle : dt) {
            if (Collections.disjoint(tri, triangle)) {
                Node[] nds = new Node[3];
                int count=0;
//                System.out.println("triangle: "+triangle.getID());
                for(Iterator<Pnt>it=triangle.iterator();it.hasNext();){
                    Pnt aPnt=it.next();
//                    System.out.print(aPnt.getIndex()+" ");
                    nds[count++] = this.theNodes.get(aPnt.getIndex());
                }
//                System.out.println();
                Triangle aTrg = new Triangle(++indtr,nds);
                theShapes.put(aTrg.getID(), aTrg);
            }
        }
        isTriangulized=true;
    }
    
    public boolean isTriangulized(){return isTriangulized;}
    
    public void putConstraintEquation(ConstraintEquation aConstraintEquation){
        this.theConstraintEquations.put(aConstraintEquation.getID(), aConstraintEquation);
    }
    
    public ConstraintEquation getConstraintEquation(int id){
        return this.theConstraintEquations.get(id);
    }

    public Map<Integer,ConstraintEquation> getConstraintEquations(){
        return this.theConstraintEquations;
    }
    
    public void setNumDOFs_Steps(int n, int nsteps){
        for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
                Node aNode=it.next();
                aNode.setNumDOFs_Steps(n, nsteps);
            }
    }
    
    public void setNumDOFs(int n){
        for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
                Node aNode=it.next();
                aNode.setNumDOFs_Steps(n, 1);
            }
    }
    
    public ArrayList<Double> getX(){
        ArrayList<Double> list = new ArrayList<Double>();
        for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
                Node aNode=it.next();
                list.add(aNode.X());
        }
        return list;
    }
    
    public ArrayList<Double> getY(){
        ArrayList<Double> list = new ArrayList<Double>();
        for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
                Node aNode=it.next();
                list.add(aNode.Y());
        }
        return list;
    }
    
    public ArrayList<Double> getZ(){
        ArrayList<Double> list = new ArrayList<Double>();
        for(Iterator<gendomain.Node>it=getNodes().values().iterator();it.hasNext();){
                Node aNode=it.next();
                list.add(aNode.Z());
        }
        return list;
    }
}
