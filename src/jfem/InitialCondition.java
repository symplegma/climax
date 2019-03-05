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
abstract public class InitialCondition {
    protected int id;
    protected static int numberOfConstraints = 0;
    protected double InitialValue=0.;
    protected Node wnode;
    protected int wdof; // local dof of wnode
    protected int type; // 0 fo disp, 1 for vel
    
    // constructor    
    public InitialCondition(){
    }
   
    
    // methods
    public double getInitialValue(){
        return InitialValue;
    }
    
    public Node getNode(){
        return wnode;
    }
    
    public int getDOF(){
        return wdof;
    }
    
    public int getID(){
        return id;
    }
    
    public int getType(){
        return type;
    }

}
