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

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import static jfem.Tracker.numberOfTrackers;

/**
 *
 * @author pchr
 */
public class TrackerNode extends Tracker{
    
    // constructor
    public TrackerNode(){}
    
    public TrackerNode(Analysis theAnalysis, Node theNode, int wdof, int resp){
        try {
            ++numberOfTrackers;
            this.theAnalysis = theAnalysis;
            this.id = numberOfTrackers;
            this.theNode = theNode;
            this.wdof = wdof;
            this.resp = resp;
            filename = "n_" + theAnalysis.getID() + "_" + theNode.getID() + "_" + wdof + "_" + resp + ".jtn";
            TrackerFile = new PrintWriter(new FileWriter(filename));
        } catch (IOException ex) {
            Logger.getLogger(Tracker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public TrackerNode(Analysis theAnalysis, Node theNode, int wdof, int resp, String Prefix){
        ++numberOfTrackers;
        this.theAnalysis = theAnalysis;
        this.id = numberOfTrackers;
        this.theNode = theNode;
        this.wdof = wdof;
        this.resp = resp;
        filename = "n_"+Prefix+'_'+theAnalysis.getID() + "_" + theNode.getID() + "_" + wdof + "_" + resp + ".jtn";
        this.f=new double[Tracker.numsteps];
    }
    
    public void write(int step){
        double data;
        switch(resp){
            case 2:data=theNode.getVelc(wdof); break;
            case 3:data=theNode.getAccl(wdof); break;
            default:data=theNode.getDisp(wdof); break;
        }
        TrackerFile.print(step*dt);
        TrackerFile.print(" ");
        TrackerFile.println(data);
    }
}
