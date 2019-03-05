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
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class LoadCase {
    private Map<Integer,Load> theLoads = new TreeMap<Integer,Load>();
    //private Vector Loads = new Vector(0);
    private int id;
    private static int numberOfLoadCases = 0;
    private int numberOfLoads = 0;
    private int numOfIncrements=1;
   // private static Vector idsOfLoadCases = new Vector(0);
    
    // constructor
    public LoadCase(int id){
        this.id=id;
        ++numberOfLoadCases;
    }
    
    public LoadCase(int id, int nIncr){
        this.id=id;
        this.numOfIncrements=nIncr;
        ++numberOfLoadCases;
    }
    
    // methods
    public int getID() {
        return id;
    }
    
    public int getNumOfIncrements() {
        return numOfIncrements;
    }
    
    public int getnumberOfLoadCases() {
        return numberOfLoadCases;
    }
    
    public Map getLoads() {
        return this.theLoads;
    }
    
    public void putLoad(Load aLoad){
        this.theLoads.put(numberOfLoads, aLoad);
        ++numberOfLoads;
    }
    
    public void setNumOfIncrements(int nIncr) {
        numOfIncrements=nIncr;
    }

}
