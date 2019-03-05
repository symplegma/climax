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
public class FundamentalSolutionData {
    private double[] r;
    private double[] outwardnormal;
    private Material theMaterial;
    private double timestep;
    private int mStep;
    private int nStep;
    private int utimedist;
    private int ptimedist;
    private int vtimedist;
    private int timepos;
    private double currentFrequency;
    
    // constructor
    public FundamentalSolutionData(){}
    
    // methods
    public void setR(double[] r_){
        int l=r_.length;
        r=new double[l];
        for(int i=0; i<l; i++){
            r[i]=r_[i];
        }
    }
    
    public void setOutwardNormal(double[] r_){
        int l=r_.length;
        outwardnormal=new double[l];
        System.arraycopy(r_, 0, outwardnormal, 0, l);
    }
    
    public void setMaterial(Material aMaterial){
        this.theMaterial=aMaterial;
    }
    
    public void setTimeStep(double time){
        this.timestep=time;
    }
    
    public void setMstep(int m){
        this.mStep=m;
    }
    
    public void setuTimeDist(int dist){
        this.utimedist=dist;
    }
    
    public void setvTimeDist(int dist){
        this.vtimedist=dist;
    }
    
    public void setpTimeDist(int dist){
        this.ptimedist=dist;
    }
    
    public void setTimePos(int pos){
        this.timepos=pos;
    }
    
    public void setNstep(int n){
        this.nStep=n;
    }
    
    public double[] getR(){
        return this.r;
    }
    
    public double[] getDR(){
        int len=this.r.length;
        double[] dr = new double[len];
        for(int i=0; i<len; i++){
            dr[i]=r[i]/this.getAbsR();
        }
        return dr;
    }
    
    public double[] getOutwardNormal(){
        return this.outwardnormal;
    }
    
    public Material getMaterial(){
        return this.theMaterial;
    }
    
    public double getAbsR(){
        double AbsR=0;
        for (int i=0; i<r.length; i++){
            AbsR+=r[i]*r[i];
        }
        AbsR=Math.sqrt(AbsR);
        return AbsR;
    }
    
    public double getTimeStep(){
        return this.timestep;
    }
    
    public int getMstep(){
        return this.mStep;
    }
    
    public int getNstep(){
        return this.nStep;
    }
    
    public int getuTimeDist(){
        return this.utimedist;
    }
    
    public int getvTimeDist(){
        return this.vtimedist;
    }
    
    public int getpTimeDist(){
        return this.ptimedist;
    }
    
    public int getTimePos(){
        return this.timepos;
    }
    
    public void setCurrentFrequency(double freq){this.currentFrequency=freq;}
    
    public double getCurrentFrequency(){return this.currentFrequency;}

}
