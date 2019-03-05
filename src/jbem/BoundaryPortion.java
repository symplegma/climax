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
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class BoundaryPortion {
    private static int numBP=0;
    private int id;
    private Map<Integer,Element> theElements = new TreeMap<Integer,Element>();
    private Domain onDomain;
    private int TrackNodeID;
    private boolean trackNodeDisp = false;
    private boolean trackNodeTract = false;
    private String nameOfBP;

    public BoundaryPortion(Domain onDomain){
        this.id=++numBP;
        this.onDomain=onDomain;
    }
    
    public BoundaryPortion(Domain onDomain, String thename){
        this.id=++numBP;
        this.onDomain=onDomain;
        this.nameOfBP=thename;
    }

    public void putElement(Element anElement){
        this.theElements.put(anElement.getID(), anElement);
    }

    public int getID(){
        return this.id;
    }

    public Map<Integer,Element> getElements(){
        return this.theElements;
    }
    
    public void printElements(){
        System.out.println("Boundary Portion with id/Name : "+this.id+"/"+this.nameOfBP);
        System.out.println("Boundary Portion TrackNodeID : "+this.TrackNodeID);
        System.out.println("Boundary Portion trackNodeDisp : "+this.trackNodeDisp);
        System.out.println("Boundary Portion trackNodeTract : "+this.trackNodeTract);
        System.out.println("Elements");
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element anElement = it.next();
            System.out.println(anElement.getID());
        }
        System.out.println("____________________________________");
    }

    public Domain getDomain(){return this.onDomain;}
    
    public void setTrackNodeID(int idn){this.TrackNodeID=idn;}
    
    public int getTrackNodeID(){return this.TrackNodeID;}
    
    public void setTrackNodeDisp(boolean t){this.trackNodeDisp=t;}
    
    public void setTrackNodeTract(boolean t){this.trackNodeTract=t;}
    
    public boolean getTrackNodeDisp(){return this.trackNodeDisp;}
    
    public boolean getTrackNodeTract(){return this.trackNodeTract;}
    
    public String getName(){return this.nameOfBP;}
}
