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

/**
 *
 * @author pchr
 */
abstract public class Tracker {
    protected static double dt=1.;
    protected static int numsteps;
    protected Analysis theAnalysis;
    protected int id;
    protected static int numberOfTrackers = 0;
    protected PrintWriter TrackerFile;
    protected String filename;
    protected int print_step=1;
    protected double[] f;
    protected int colm=1;
    protected int count=0;
    
    protected Node theNode;
    protected int wdof;
    protected int resp;
//    protected static String prefix;
    // constructor
    public Tracker(){}
    
    
    //methods
    public int getID(){
        return this.id;
    }
    
    public void close(){
//        System.out.println("CLOSE TRACKER: "+this.id+" on (nodeid,dof,resp)=("+theNode.getID()+","+wdof+","+resp+")");
        this.TrackerFile.close();
    }
    
    public void setDT(double h){Tracker.dt=h;}
    
    public void setNumSteps(int n){Tracker.numsteps=n;}
    
//    public void setPrefix(String Prefix){Tracker.prefix=Prefix;}
    
    public void setStep(int step){this.print_step=step;}
    
    public int getStep(){return this.print_step;}
    
    public void setColumnsFFT(int c){colm=c;}
    
    abstract public void write(int step);
    
    public void write(double dstep) {
        int step=(int) Math.round(dstep/dt);
//        int step=(int) Math.floor(dstep/dt);
        write(step);
    }
    
    public boolean isDouble( String input ){
        boolean theb=false;
        try{
             Double.parseDouble( input );
            theb= true;
        }catch(Exception e){
            theb= false;
        }
        return theb;
    }
    
    public String getFileName(){return this.filename;}
}
