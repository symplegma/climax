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
import jmat.AbstractMatrix;
import static jfem.Tracker.numberOfTrackers;

/**
 *
 * @author pchr
 */
public class TrackerElement extends Tracker{
    private Element theElement;
    
    // constructor
    public TrackerElement(Analysis theAnalysis, Element theElement){
        try {
            ++numberOfTrackers;
            this.theAnalysis = theAnalysis;
            this.id = numberOfTrackers;
            this.theElement=theElement;

            filename = "e_" + theAnalysis.getID() + "_" +theElement.getID() +".jte";
            TrackerFile = new PrintWriter(new FileWriter(filename));
        } catch (IOException ex) {
            Logger.getLogger(Tracker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public TrackerElement(Analysis theAnalysis, Element theElement, String Prefix){
        try {
            ++numberOfTrackers;
            this.theAnalysis = theAnalysis;
            this.id = numberOfTrackers;
            this.theElement=theElement;

            filename = "e_"+Prefix+'_'+theAnalysis.getID() + "_" + theElement.getID() +".jte";
            TrackerFile = new PrintWriter(new FileWriter(filename));
        } catch (IOException ex) {
            Logger.getLogger(Tracker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void write(int step){
        AbstractMatrix data;
        data = theElement.getF();
        int ndof=data.getRowDimension();
        TrackerFile.print(step*dt);
        TrackerFile.print(" ");
        for(int i=0; i<ndof; i++){
            TrackerFile.print(data.get(i, 0));
            TrackerFile.print(" ");
        }
        TrackerFile.print(theElement.getuKu()/2.+" "+theElement.getvMv()/2.);
        TrackerFile.println();
    }

}
