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

import geom.Point;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Node extends Point{
    private double[][] dofValues; // dofValues[dof][step]
    private double[][] dofTDValues; // dofValues[dof][step]
    private double[][] dof2TDValues; // dofValues[dof][step]
    private static int numberOfNodes = 0;
    private Map<Integer,Element> ConnectedElements = new TreeMap<Integer,Element>();
    private Map<DOFGroup.DOFDescription,Element> DOFGroups = new TreeMap<DOFGroup.DOFDescription,Element>();
    private int[] uEFTable = new int[6];
    private boolean boundary=false;
    private static int ndofs=0;
    private int ndofs_ofNode=0;
    private int nSteps;
    
     public Node(){};
     
     public Node(int id,double v[]){
        dofValues=new double[1][1];
        this.coordinates = new double[v.length];
        for (int i=0;i<coordinates.length; i++){
            this.coordinates[i] = v[i];
        }
        this.id =id;
        for(int i=0; i<6; i++){
            uEFTable[i]=0;
        }
        ++numberOfNodes;
     }
     
     public Node(int id,double x){
        dofValues=new double[1][1];
        this.coordinates = new double[2];
        this.coordinates[0] = x;
        this.coordinates[1] = 0;
        this.id =id;
        for(int i=0; i<6; i++){
            uEFTable[i]=0;
        }
        ++numberOfNodes;
     }
     
     public Node(int id,double x, double y){
        dofValues=new double[1][1];
        this.coordinates = new double[2];
        this.coordinates[0] = x;
        this.coordinates[1] = y;
        this.id =id;
        for(int i=0; i<6; i++){
            uEFTable[i]=0;
        }
        ++numberOfNodes;
     }
     
     public Node(int id,double x, double y, double z){
        dofValues=new double[1][1];
        this.coordinates = new double[3];
        this.coordinates[0] = x;
        this.coordinates[1] = y;
        this.coordinates[2] = z;
        this.id =id;
        for(int i=0; i<6; i++){
            uEFTable[i]=0;
        }
        ++numberOfNodes;
     }
     
     public void setBoundary(){this.boundary=true;}
     
     public boolean getBoundary(){return this.boundary;}
     
     public void putElement(Element aElement) {
        this.ConnectedElements.put(aElement.getID(), aElement);
    }
     
     public double[] getCoords(){
        return this.coordinates;
    }
     
    public void updateVals(double[] x, int step){
        for(int i=0;i<dofValues.length;i++){
            dofValues[i][step]=x[i];
        }
    }
    
    public void updateValDOF(double val, int step, int dof){
        dofValues[dof][step]=val;
    }
     
    public void accumulateVals(double[] x, int step){
        for(int i=0;i<dofValues.length;i++){
            dofValues[i][step]=+x[i];
        }
    }
    
    public void updateVals(double x, int step){
        dofValues[0][step]=x;
    }
     
    public void updateTDVals(double[] x, int step){
        for(int i=0;i<dofTDValues.length;i++){
            dofTDValues[i][step]=x[i];
        }
    }
     
    public void update2TDVals(double[] x, int step){
        for(int i=0;i<dof2TDValues.length;i++){
            dof2TDValues[i][step]=x[i];
        }
    }
     
     public void setNumDOFs(int n){setNumDOFs_Steps(n,1);}
     
     public void setNumDOFs_Steps(int n, int nsteps){
         this.dofValues=new double[n][nsteps];
         this.nSteps=nsteps;
     }
     
     public void setNumDOFsTD_Steps(int n, int nsteps){
         this.dofValues=new double[n][nsteps];
         this.dofTDValues=new double[n][nsteps];
         this.dof2TDValues=new double[n][nsteps];
     }
     
     public void initInitialialConditions(double val){
         for(int i=0;i<dofValues.length;i++)dofValues[i][0]=val;
     }
     
     public void initInitialialConditions(double val, int dof){
         dofValues[dof][0]=val;
     }
     
     public void initInitialialConditionsTD(double val, int dof){
         dofTDValues[dof][0]=val;
     }
     
     public int getNumDOFS(){return this.dofValues.length;}
     
     public int getNumSteps(){return this.nSteps;}
     
     public void setuEFTable(int n){
         // this was originally here
         // however after DEM implementation I have also set
         // setNdofs_ofNode(int[] dofs)
        this.uEFTable=new int[dofValues.length];
        for(int i=0; i<dofValues.length; i++){
            this.uEFTable[i]=n+i;
        }
    }
     
     public void setNdofs_ofNode(int[] dofs) {
         // right now I can not see the difference between ndofs and ndofs_ofNode
        for(int i=0; i<dofs.length; i++ ){
            if(this.uEFTable[dofs[i]-1]==0){
                this.uEFTable[dofs[i]-1]=++ndofs;
                ++ndofs_ofNode;
            }
        }
    }
     
     public int[] getuEFTable(){
        return this.uEFTable;
    }
     
    public double getDisps(int step,int wdof){
        return this.dofValues[wdof][step];
    }
     
    public double getSolution(int step,int wdof){
        return this.dofValues[wdof][step];
    }
    
    public double getSolution(int wdof){
        return getSolution(0,wdof);
    }
    
    public double getSolution(){
        return getSolution(0,0);
    }
     
    public double getSolutionTD(int step,int wdof){
        return this.dofTDValues[wdof][step];
    }
     
    public double getSolution2TD(int step,int wdof){
        return this.dof2TDValues[wdof][step];
    }
        
    public double[] getSolutionHistory(int wdof){
        double[] retv=new double[this.nSteps];
        System.arraycopy(dofValues[wdof], 0, retv, 0, this.nSteps);
        return retv;
    }
    
     public static Node parsePoint(Point aPoint){
         Node pNode=new  Node(aPoint.getID(), aPoint.getCoordinates());
         return pNode;
     }
     
     
//     public void putElement(Element aElement) {
//         this.ConnectedElements.put(aElement.getID(), aElement);
//     }
}
