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

/**
 *
 * @author pchr
 */
public class Combo {
    private int id;
    private int[] states;
    private int responseID;
    private double[] coefs;
    
    public Combo(int id){
        this.id=id;
    }
    
    public void setComboStates(int[] states){
        this.states = new int[states.length];
        this.coefs = new double[states.length];
        System.arraycopy(states, 0, this.states, 0, states.length);
        for(int i=0;i<coefs.length;i++)coefs[i]=1.0;
    }
    
    public void setComboStates(int[] states, double[] coefs){
        this.states = new int[states.length];
        this.coefs = new double[coefs.length];
        System.arraycopy(states, 0, this.states, 0, states.length);
        System.arraycopy(states, 0, this.coefs, 0, coefs.length);
    }
    
    public int[] getComboStates(){return this.states;}
    
    public double[] getComboCoefs(){return this.coefs;}
    
    public int getID(){return this.id;}
    
    public void setResponseID(int responseID){this.responseID=responseID;}
    
    public int getResponseID(){return this.responseID;}
    
}
