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

/**
 *
 * @author pchr
 */
  abstract public class Material {
    protected int id;
    protected static int numberOfMaterials = 0;
    protected int type;
    
    // constructor
    public Material(){
    }
    
    public Material(int id){
        this.id =id;
        ++numberOfMaterials;
    }
    
    
    // methods
    public int getnumberOfMaterials() {
        return numberOfMaterials;
    }
    
    public int getID() {
        return id;
    }
    
    public double getType(){
        return this.type;
    }



    
    
}
