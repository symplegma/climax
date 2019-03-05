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
public class Tracker2DLineElems extends Tracker{
    
    // constructor
    public Tracker2DLineElems(Analysis theAnalysis,double dt, String Prefix){
        this.dt=dt;
        try {
            ++numberOfTrackers;
            this.theAnalysis = theAnalysis;
            this.id = numberOfTrackers;

//            String filename = "r_" + theAnalysis.getID() +".vtk";
            filename = "r_"+Prefix+'_'+theAnalysis.getID() +".tdl";
            TrackerFile = new PrintWriter(new FileWriter(filename));
        } catch (IOException ex) {
            Logger.getLogger(Tracker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Override
    public void write(int step) {
        TrackerFile.println("Dispalcements u"+step);
        
//        for(Iterator<Node> it=theAnalysis.getDomain().getNodes().values().iterator(); it.hasNext();){
//            Node aNode=it.next();
//            TrackerFile.println(aNode.getID()+" "+
//                    aNode.getCoords()[0]+" "+
//                    aNode.getCoords()[1]+" "+
//                    aNode.getDisp(1)+" "+
//                    aNode.getDisp(2)+" "+
//                    aNode.getDisp(6));
//        }
        
        for(Iterator<Element> et=theAnalysis.getDomain().getElements().values().iterator(); et.hasNext();){
            Element anElement=et.next();
            double ElastE=anElement.getuKu()/2.;
            double kinE=anElement.getvMv()/2.;
            for(Iterator<Node> it=anElement.ElementNodes.values().iterator(); it.hasNext();){
                Node aNode=it.next();
                TrackerFile.println(aNode.getID()+" "+
                        aNode.getCoords()[0]+" "+
                        aNode.getCoords()[1]+" "+
                        aNode.getDisp(1)+" "+
                        aNode.getDisp(2)+" "+
                        aNode.getDisp(6)+" "+
                        ElastE+" "+
                        kinE+" "+
                        aNode.getVelc(1)+" "+
                        aNode.getVelc(2)+" "+
                        aNode.getVelc(6));
            }
        }
        
        TrackerFile.println();
    }
    
}
