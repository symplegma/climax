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

package jfem;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import static jfem.Tracker.numberOfTrackers;

/**
 *
 * @author pchr
 */
public class TrackerVTK extends Tracker{
    // exei graftei elafros proxeira
    // isxyei mono gia 4komva stoixeia
    // ENTELOS EXIDIKEYMENO
    int[] map;
    private int variable=0; // 0: nrg, 1:disp, 2: velc
    
    // constructor
    public TrackerVTK(Analysis theAnalysis,double dt, String Prefix, int var){
        this.dt=dt;
        this.variable=var;
        try {
            ++numberOfTrackers;
            this.theAnalysis = theAnalysis;
            this.id = numberOfTrackers;

//            String filename = "r_" + theAnalysis.getID() +".vtk";
            filename = "r_"+Prefix+'_'+theAnalysis.getID()+'_'+variable+".vtk";
            TrackerFile = new PrintWriter(new FileWriter(filename));
            int numN=theAnalysis.getDomain().getNumNodes();
            map = new int[numN];
            TrackerFile.println("# vtk DataFile Version 3.0");
            TrackerFile.println("jfem vtk output");
            TrackerFile.println("ASCII");
            TrackerFile.println("DATASET UNSTRUCTURED_GRID");
            TrackerFile.println();
            TrackerFile.println("POINTS "+numN+" float");
            int count=0;
            for(Iterator<Node> it=theAnalysis.getDomain().getNodes().values().iterator(); it.hasNext();){
                Node aNode=it.next();
                map[count]=aNode.getID();
                TrackerFile.println(aNode.getCoords()[0]+" "+
                                aNode.getCoords()[1]+" "+
                                aNode.getCoords()[2]);
                count+=1;
            }
            TrackerFile.println();
            numN=theAnalysis.getDomain().getNumQuadElement();
            int nn=numN*5;
            int numTruss = theAnalysis.getDomain().getNumLineElement();
            int nnn=numTruss*3;
            TrackerFile.println("CELLS "+(numN+numTruss)+" "+(nn+nnn));
            int n1,n2,n3,n4;
            int m1 = 0,m2 = 0,m3 = 0,m4 = 0;
            for(Iterator<Element> it=theAnalysis.getDomain().getElements().values().iterator(); it.hasNext();){
                Element aElement=it.next();
                if(aElement.getClass()==PlaneStressQuad.class || aElement.getClass()==PlaneStrainQuad.class){
                    n1=aElement.getNodeHierarchy(1).getID();
                    n2=aElement.getNodeHierarchy(2).getID();
                    n3=aElement.getNodeHierarchy(3).getID();
                    n4=aElement.getNodeHierarchy(4).getID();

                    for(int i=0;i<map.length;i++){
                        if(n1==map[i]){
                            m1=i;
                            break;
                        }
                    }

                    for(int i=0;i<map.length;i++){
                        if(n2==map[i]){
                            m2=i;
                            break;
                        }
                    }

                    for(int i=0;i<map.length;i++){
                        if(n3==map[i]){
                            m3=i;
                            break;
                        }
                    }

                    for(int i=0;i<map.length;i++){
                        if(n4==map[i]){
                            m4=i;
                            break;
                        }
                    }
                    TrackerFile.println(4+" "+m1+" "+m2+" "+m3+" "+m4);
                }
            }
            for(Iterator<Element> it=theAnalysis.getDomain().getElements().values().iterator(); it.hasNext();){
                Element aElement=it.next();
                if(aElement.getClass()==Truss2d.class || aElement.getClass()==EBeam2d.class || aElement.getClass()==EBeam3d.class){
                    n1=aElement.getNodeHierarchy(1).getID();
                    n2=aElement.getNodeHierarchy(2).getID();

                    for(int i=0;i<map.length;i++){
                        if(n1==map[i]){
                            m1=i;
                            break;
                        }
                    }

                    for(int i=0;i<map.length;i++){
                        if(n2==map[i]){
                            m2=i;
                            break;
                        }
                    }
                    TrackerFile.println(2+" "+m1+" "+m2);
                }
            }
            TrackerFile.println();
            TrackerFile.println("CELL_TYPES"+" "+(numN+numTruss));
            for(int i=1;i<=numN;i++){
                TrackerFile.println(9);
            }
            for(int i=1;i<=numTruss;i++){
                TrackerFile.println(3);
            }
            TrackerFile.println();
            TrackerFile.println("POINT_DATA"+" "+theAnalysis.getDomain().getNumNodes());
            TrackerFile.println();
        } catch (IOException ex) {
            Logger.getLogger(Tracker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Override
    public void write(int step) {
        double[] data=new double[3];
//        
        switch(variable){
            case 1: 
                TrackerFile.println("VECTORS u"+step+" float");
                for(int i=0;i<map.length;i++){
                    data[0]=theAnalysis.getDomain().getNode(map[i]).getDisp(1);
                    data[1]=theAnalysis.getDomain().getNode(map[i]).getDisp(2);
                    data[2]=theAnalysis.getDomain().getNode(map[i]).getDisp(3);
//                    data[2]=theAnalysis.getDomain().getNode(map[i]).getDisp(6);
                    TrackerFile.println(data[0]+" "+data[1]+" "+data[2]);
                }
                TrackerFile.println();
                break;
            case 2: 
                TrackerFile.println("VECTORS u"+step+" float");
                for(int i=0;i<map.length;i++){
                    data[0]=theAnalysis.getDomain().getNode(map[i]).getVelc(1);
                    data[1]=theAnalysis.getDomain().getNode(map[i]).getVelc(2);
                    data[2]=theAnalysis.getDomain().getNode(map[i]).getVelc(3);
//                    data[2]=theAnalysis.getDomain().getNode(map[i]).getDisp(6);
                    TrackerFile.println(data[0]+" "+data[1]+" "+data[2]);
                }
                TrackerFile.println();
                break;
            default: 
            TrackerFile.println("VECTORS nrg"+step+" float");
            for(int i=0;i<map.length;i++){
                data[0]=0.0;
                data[1]=0.0;
                for(Iterator<Element> it=theAnalysis.getDomain().getNode(map[i]).getElements().values().iterator(); it.hasNext();){
                    Element theElement = it.next();
                    data[0]+=theElement.getuKu()/2.;
                    data[1]+=theElement.getvMv()/2.;
                }
                data[2]=data[0]+data[1];
//                TrackerFile.println(data[0]/theAnalysis.getDomain().getNode(map[i]).getNumElement()+" "+0.0+" "+0.0);
                TrackerFile.println(data[0]+" "+data[1]+" "+data[2]);
            }
            TrackerFile.println();
                break;
        }
    }
}