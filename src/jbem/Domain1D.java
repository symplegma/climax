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

package jbem;

import java.util.Iterator;
//import java.util.Map;

/**
 *
 * @author pchr
 */
public class Domain1D extends Domain{
    int idFnode;
    int idSnode;
    // constructor
    public Domain1D(int id){
        this.id=id;
    }
    
    
    @Override
    public void putNode(Node aNode){
        int size=this.theNodes.size();
        if(size>1){
            System.err.println("One Dimensional Domain can not have more than 2 nodes");
            System.exit(111);
        }
        if(size==0){idFnode=aNode.getID();}else{idSnode=aNode.getID();}
        this.theNodes.put(aNode.getID(), aNode);
    }
    
    public double getL(){
        double L=0.;
        double[] coord1=new double[0];
        double[] coord2=new double[0];
        int i=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            for(Iterator<Node> nit=theElement.getNodes().values().iterator(); nit.hasNext();){
                Node node_of_elem = nit.next();
                if(i==0){coord1=node_of_elem.getCoordinates();}
                if(i==1){coord2=node_of_elem.getCoordinates();}
                ++i;
            }
        }
        int s1=coord1.length;
        int s2=coord2.length;
        double[] coord13=new double[3];
        double[] coord23=new double[3];
        for(i=0; i<3; i++){
            if(i<s1){coord13[i]=coord1[i];}else{coord13[i]=0.;}
            if(i<s2){coord23[i]=coord2[i];}else{coord23[i]=0.;}
        }
        L=Math.sqrt(
                (coord23[0]-coord13[0])*(coord23[0]-coord13[0])
               +(coord23[1]-coord13[1])*(coord23[1]-coord13[1])
               +(coord23[2]-coord13[2])*(coord23[2]-coord13[2]));
        return L;
    }
    
    public double getCOS(){
        double c=0;
        c=(theNodes.get(idSnode).getCoordinates()[0]-theNodes.get(idFnode).getCoordinates()[0])/this.getL();
        return c;
    }
    
    public double getSIN(){
        double c=0;
        c=(theNodes.get(idSnode).getCoordinates()[1]-theNodes.get(idFnode).getCoordinates()[1])/this.getL();
        return c;
    }
    
    public int getNodeHier(int hier){
        int node_id=0;
        if(hier==1){node_id=this.idFnode;}else{node_id=this.idSnode;}
        return node_id;
    }

}
