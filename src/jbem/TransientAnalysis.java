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

/**
 *
 * @author pchr
 */
abstract public class TransientAnalysis extends Analysis{
    int currentStep;
    
    // constructor
    public TransientAnalysis(){
        currentStep=0;
    }
    
    public int getCurrentStep(){return this.currentStep;}
    
    public void printDistances(){
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();

            System.out.println("Node distances matrix");
            System.out.println("=====================");
            theDomain.getLmat().print(12, 6);
        }
        
    }
    
    public void printDistances(double vel, double dt){
        for(Iterator<Domain> dit=this.theDomains.values().iterator(); dit.hasNext();){
            Domain theDomain = dit.next();

            System.out.println("Node distances matrix ! L/(c*h)");
            System.out.println("===============================");
            theDomain.getLmat().times(1./(dt*vel)).print(12, 6);
        }
    }

}
