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
public class Load {
    protected static int numberOfLoads = 0;
    protected double LoadValue;
    protected int wnode;
    protected int wdof; // local dof of wnode
    protected int timeStep=0;
    
    // constructor
    public Load(){}
    
    public Load(double LoadValue, int wnode, int wdof){
        this.LoadValue=LoadValue;
        this.wnode=wnode;
        this.wdof=wdof;
    }
    
    public Load(double LoadValue, int wnode, int wdof,int timeStep){
        this.LoadValue=LoadValue;
        this.wnode=wnode;
        this.wdof=wdof;
        this.timeStep=timeStep;
    }
    
    // methods
    public double getLoadValue() {
        return LoadValue;
    }
    
    public int getwnode() {
        return wnode;
    }
    
    public int getwdof() {
        return wdof;
    }
    
    public int getTime() {
        return timeStep;
    }


}
