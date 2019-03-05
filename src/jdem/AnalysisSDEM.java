/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import gendomain.Node;
import geom.Circle;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 *
 * @author pchr
 */
public class AnalysisSDEM extends Analysis{
    // SDEM stands for Stress-based Discrete Element Method 
    
    //1. Egholm, D. L. (2007), A new strategy for discrete element numerical models: 1. Theory, J. Geophys. Res., 112. 
    //2. Egholm, D. L., Sandiford, M., Clausen, O. R., Nielsen, S. B. (2007), A new strategy for discrete element numerical models: 2. Sandbox applications, J. Geophys. Res., 112. 
    
    // Users guide to DISC: Stress-based Discrete Element Method software, Anders Damsgaard Christensen (https://adamsgaard.dk/)
    public AnalysisSDEM(DEMdomain dom){
        this.theDomain=dom;
    }

    @Override
    public int analyse(){
        int tag=0;
        this.theDomain.DomainInitialConditions();
       
        for(int istep=1; istep<=numSteps; istep++){
            this.step=istep;

            // 1. Predictor
            for(Iterator<Particle>it=theDomain.getParticles().values().iterator();it.hasNext();){
                 Particle P=it.next();
                 Node Nd=P.getCentroid();
                 double[] vals= new double[6];
                 vals[0]=Nd.getSolutionTD(step-1, 0)+timeStep*P.getFx_o()/P.get_m();
                 vals[1]=Nd.getSolutionTD(step-1, 1)+timeStep*P.getFy_o()/P.get_m();
                 vals[5]=Nd.getSolutionTD(step-1, 5)+timeStep*P.getMz_o()/P.get_Jm();
                 Nd.updateTDVals(vals, step);

                 vals[0]=Nd.getSolution(step-1, 0)+Nd.getSolutionTD(step-1, 0)*timeStep+0.5*timeStep*timeStep*P.getFx_o()/P.get_m();
                 vals[1]=Nd.getSolution(step-1, 1)+Nd.getSolutionTD(step-1, 1)*timeStep+0.5*timeStep*timeStep*P.getFy_o()/P.get_m();
                 vals[5]=Nd.getSolution(step-1, 5)+Nd.getSolutionTD(step-1, 5)*timeStep+0.5*timeStep*timeStep*P.getMz_o()/P.get_Jm();
                 
                 vals[0]=Nd.getSolution(step-1, 0)+Nd.getSolutionTD(step-1, 0)*timeStep+0.5*timeStep*timeStep*P.getFx_o()/P.get_m();
                 vals[1]=Nd.getSolution(step-1, 1)+Nd.getSolutionTD(step-1, 1)*timeStep+0.5*timeStep*timeStep*P.getFy_o()/P.get_m();
                 vals[5]=Nd.getSolution(step-1, 5)+Nd.getSolutionTD(step-1, 5)*timeStep+0.5*timeStep*timeStep*P.getMz_o()/P.get_Jm();
                 Nd.accumulateVals(vals, step);
            }

            // 2. Contact Searching
            List<Map.Entry<Integer,Integer>> OverlapList= new ArrayList<>();
            System.out.println("Time step: "+istep);
            for (City theCity : theDomain.getCities().values()) {
                //System.out.println("City: "+theCity.getID());
                for (Block theBlock : theCity.getBlocks().values()) {
                    //System.out.println("Block: "+theBlock.getID());
                    for(Iterator<Particle>it=theBlock.getParticles().values().iterator();it.hasNext();){
                        Particle Pi=it.next();
                        for(Iterator<Particle>jt=theBlock.getParticles().values().iterator();jt.hasNext();){
                             Particle Pj=jt.next();
                             double xi,yi,xj,yj;
                             xi=Pi.getCentroid().X()+Pi.getCentroid().getDisps(istep, 0);
                             yi=Pi.getCentroid().Y()+Pi.getCentroid().getDisps(istep, 1);
                             xj=Pj.getCentroid().X()+Pj.getCentroid().getDisps(istep, 0);
                             yj=Pj.getCentroid().Y()+Pj.getCentroid().getDisps(istep, 1);
                             double dij=Math.sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))-((Circle)Pi.getShape()).getR()-((Circle)Pj.getShape()).getR();
                             if( (Pi.getID()!=Pj.getID()) && (dij<0.0) ){
                                 Map.Entry<Integer,Integer> pair=new AbstractMap.SimpleEntry<>(Pi.getID(), Pj.getID());
                                 Map.Entry<Integer,Integer> pair_inv=new AbstractMap.SimpleEntry<>(Pj.getID(), Pi.getID());
                                 if(!OverlapList.contains(pair)){
                                     if(!OverlapList.contains(pair_inv)){
                                         OverlapList.add(pair);
                                     }
                                 }

                             }
                        }
                    }
                }
            }
//            for (Map.Entry<Integer,Integer> temp : OverlapList) {
//                System.out.println(temp.getKey()+" - "+temp.getValue());
//            }
            
            // 3. Strain Rates
            
            // 4. Stress tensors
            
            // 5. Forces
            for (Map.Entry<Integer,Integer> temp : OverlapList) {
                Particle Pi=theDomain.getParticle(temp.getKey());
                Particle Pj=theDomain.getParticle(temp.getValue());
                double xi,yi,xj,yj;
                xi=Pi.getCentroid().X()+Pi.getCentroid().getDisps(istep, 0);
                yi=Pi.getCentroid().Y()+Pi.getCentroid().getDisps(istep, 1);
                xj=Pj.getCentroid().X()+Pj.getCentroid().getDisps(istep, 0);
                yj=Pj.getCentroid().Y()+Pj.getCentroid().getDisps(istep, 1);
                double d=Math.sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
                double Ri=((Circle)Pi.getShape()).getR();
                double Rj=((Circle)Pj.getShape()).getR();
                double dij=-d+Ri+Rj;
                double nx,ny,xc,yc;
                nx=(xj-xi)/d; ny=(yj-yi)/d;
                xc=xi+(Ri-0.5*dij)*nx; yc=yi+(Ri-0.5*dij)*ny;
                System.out.println(temp.getKey()+" - "+temp.getValue()+": "+dij+" at [xc,yc]=["+xc+","+yc+"]");
                double kn=100.0; // shall be computed by particle's material
                theDomain.getParticle(temp.getKey()).addFx_o(-kn*nx*dij);
                theDomain.getParticle(temp.getKey()).addFy_o(-kn*ny*dij);
                theDomain.getParticle(temp.getValue()).addFx_o(kn*nx*dij);
                theDomain.getParticle(temp.getValue()).addFy_o(kn*ny*dij);
                
//                xi=Pi.getCentroid().X()+Pi.getCentroid().getSolutionTD(istep, 0);
//                yi=Pi.getCentroid().Y()+Pi.getCentroid().getSolutionTD(istep, 1);
//                xj=Pj.getCentroid().X()+Pj.getCentroid().getSolutionTD(istep, 0);
//                yj=Pj.getCentroid().Y()+Pj.getCentroid().getSolutionTD(istep, 1);
//                d=-Math.sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
//                kn=1.0;
//                theDomain.getParticle(temp.getKey()).addFx_o(-kn*nx*d);
//                theDomain.getParticle(temp.getKey()).addFy_o(-kn*ny*d);
//                theDomain.getParticle(temp.getValue()).addFx_o(kn*nx*d);
//                theDomain.getParticle(temp.getValue()).addFy_o(kn*ny*d);
            }
//            // 6. Corrector
//            for(Iterator<Particle>it=theDomain.getParticles().values().iterator();it.hasNext();){
//                 Particle P=it.next();
//                 P.setFx_o(P.getFx_o());
//                 P.setFy_o(P.getFy_o());
//                 P.setMz_o(P.getMz_o());
//            }
       }
       return tag;
    }
    
}
