/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import geom.Point;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author pchr
 */
public class InterfaceNode extends Point{
    private static int numNodes=0;
    private int numOfNodeIElements=0;
    private double[] un;
    private double[] un_pos;
    private double[] un_neg;
    private double[] ut;
    private double[] ut_pos;
    private double[] ut_neg;
    private double unh;
    private double uth;
    private double[][] z;
    private double[][] s;
    private double[][] p;
    private int unEFTable;
    private int un_posEFTable;
    private int un_negEFTable;
    private int utEFTable;
    private int ut_posEFTable;
    private int ut_negEFTable;
    private int[] zEFTable;
    private int[] sEFTable;
    private int[] s_posEFTable;
    private int[] s_negEFTable;
    private int[] pEFTable;
    private int[] p_posEFTable;
    private int[] p_negEFTable;
    private Node BaseNode=null;
    private InterfaceNode TwinPoint=null;
    private Element onElement = null;
    private double[] normal;
    private Map<Integer,InterfaceElement> connectedIElements = new HashMap<Integer,InterfaceElement>();
    private Map<Integer,Integer> connectedIElementsHierarchy = new HashMap<Integer,Integer>();
    private boolean mainNode=false;
    private boolean inequation=false;
    
    private double[][] s_pos;
    private double[][] p_pos;
    private double[][] s_neg;
    private double[][] p_neg;
    
    private static boolean PlasticIncrement=true;

    // constructor
    public InterfaceNode(){}
    
    public InterfaceNode(double[] coords){
        ++numNodes;
        this.id=numNodes;
        //this.coordinates = new double[coords.length];
        this.coordinates = new double[3];
        for(int i=0; i<3; i++){
            if(i<coords.length){this.coordinates[i]=coords[i];}
            else{this.coordinates[i]=0.;}
        }
    }

    public InterfaceNode(double[] coords, Node baseNode){
        ++numNodes;
        this.id=numNodes;
        this.BaseNode=baseNode;
        //this.coordinates = new double[coords.length];
        this.coordinates = new double[3];
        for(int i=0; i<3; i++){
            if(i<coords.length){this.coordinates[i]=coords[i];}
            else{this.coordinates[i]=0.;}
        }
    }

    public InterfaceNode(double[] coords, Element onElement){
        ++numNodes;
        this.id=numNodes;
        this.onElement=onElement;
        //this.coordinates = new double[coords.length];
        this.coordinates = new double[3];
        for(int i=0; i<3; i++){
            if(i<coords.length){this.coordinates[i]=coords[i];}
            else{this.coordinates[i]=0.;}
        }
    }

    // methods
    public void init_ut(int steps){
        this.ut=new double[steps];
    }
    
    public void init_ut_pos(int steps){
        this.ut_pos=new double[steps];
    }
    
    public void init_ut_neg(int steps){
        this.ut_neg=new double[steps];
    }

    public void init_un(int steps){
        this.un=new double[steps];
    }
    
    public void init_un_pos(int steps){
        this.un_pos=new double[steps];
    }
    
    public void init_un_neg(int steps){
        this.un_neg=new double[steps];
    }

    public void init_z(int steps, int num){
        this.z=new double[steps][num];
        this.zEFTable = new int[num];
        for(int i=0;i<num;i++)for(int j=0;j<steps;j++){z[j][i]=1.0;}
    }

    public void init_s(int steps, int num){
        //this.s=new double[steps][num];
        //this.sEFTable = new int[num];
        this.init_s_pos(steps, num);
        this.init_s_neg(steps, num);
    }
    
    public void init_s_pos(int steps, int num){
        this.s_pos=new double[steps][num];
        this.s_posEFTable = new int[num];
    }
    
    public void init_s_neg(int steps, int num){
        this.s_neg=new double[steps][num];
        this.s_negEFTable = new int[num];
    }
    
    public void init_p(int steps, int num){
        //this.p=new double[steps][num];
        //this.pEFTable = new int[num];
        this.init_p_pos(steps, num);
        this.init_p_neg(steps, num);
    }
    
    public void init_p_pos(int steps, int num){
        this.p_pos=new double[steps][num];
        this.p_posEFTable = new int[num];
    }
    
    public void init_p_neg(int steps, int num){
        this.p_neg=new double[steps][num];
        this.p_negEFTable = new int[num];
    }

//    public double[] getut(){
//        double[] ret = null;
//        if(this.utEFTable>0){
//            ret = new double[ut.length];
//            for(int i=1;i<ret.length;i++){
//                ret[i] =ut[i];
//            }
//        }else{
//            if(this.ut_neg.length==this.ut_pos.length){
//                ret = new double[ut_neg.length];
//                for(int i=1;i<ret.length;i++){
//                    ret[i] =ret[i-1] + ut_pos[i]-ut_neg[i];
//                }
//            }else{
//                System.err.println("Error in gets of InterfaceNode.getut . STOP");
//                System.exit(this.id);
//            }
//        }
//        return ret;
//    }
//    
//    public double[] getut_pos(){
//        return this.ut_pos;
//    }
//    
//    public double[] getut_neg(){
//        return this.ut_neg;
//    }
    
    public double[] getut(){
        double[] vals = null;
        double L1,L2;
        int n1,n2;
        if(this.utEFTable>0){
            vals=new double[this.ut.length];
            for(int i=0;i<ut.length;i++){
                vals[i]=0.;
            }

            if(this.utEFTable!=0){
                System.arraycopy(ut, 0, vals, 0, this.ut.length);
            }else{
                if(this.isMain()){
                    L1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getLength();
                    if(this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID()==this.id){
                        n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                    }else{
                        n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                    }

                    L2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getLength();
                    if(this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID()==this.id){
                        n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                    }else{
                        n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                    }
                }else{
                    L1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getLength();
                    if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getTwin().getID()==this.id){
                        n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                    }else{
                        n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                    }

                    L2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getLength();
                    if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getTwin().getID()==this.id){
                        n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                    }else{
                        n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                    }
                }

                for(int i=0;i<un.length;i++){
                    if(this.isMain()){
                        vals[i]=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getNodes().get(n1).getut()[i]*L2/(L1+L2)+
                            this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getNodes().get(n2).getut()[i]*L1/(L1+L2);
                    }else{
                        vals[i]=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getNodes().get(n1).getTwin().getut()[i]*L2/(L1+L2)+
                            this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getNodes().get(n2).getTwin().getut()[i]*L1/(L1+L2);
                    }

                }
            }
        }else{
            if(this.ut_neg.length==this.ut_pos.length){
                vals = new double[ut_neg.length];
                for(int i=1;i<vals.length;i++){
                    vals[i] =vals[i-1] + ut_pos[i]-ut_neg[i];
                }
            }else{
                System.err.println("Error in gets of InterfaceNode.getut . STOP");
                System.exit(this.id);
            }
        }
        return vals;
    }
    
    public double[] getut_pos(){
        double[] vals;
        double L1,L2;
        int n1,n2;
        vals=new double[this.ut_pos.length];
        for(int i=0;i<ut_pos.length;i++){
            vals[i]=0.;
        }
        
        if(this.ut_posEFTable!=0){
            System.arraycopy(ut_pos, 0, vals, 0, this.ut_pos.length);
        }else{
            if(this.isMain()){
                L1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getLength();
                if(this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID()==this.id){
                    n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                }else{
                    n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                }
                
                L2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getLength();
                if(this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID()==this.id){
                    n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                }else{
                    n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                }
            }else{
                L1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getLength();
                if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getTwin().getID()==this.id){
                    n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                }else{
                    n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                }
                
                L2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getLength();
                if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getTwin().getID()==this.id){
                    n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                }else{
                    n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                }
            }
            
            for(int i=0;i<un.length;i++){
                if(this.isMain()){
                    vals[i]=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getNodes().get(n1).getut_pos()[i]*L2/(L1+L2)+
                        this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getNodes().get(n2).getut_pos()[i]*L1/(L1+L2);
                }else{
                    vals[i]=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getNodes().get(n1).getTwin().getut_pos()[i]*L2/(L1+L2)+
                        this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getNodes().get(n2).getTwin().getut_pos()[i]*L1/(L1+L2);
                }
                
            }
        }
        return vals;
    }
    
    public double[] getut_neg(){
        double[] vals;
        double L1,L2;
        int n1,n2;
        vals=new double[this.ut_neg.length];
        for(int i=0;i<ut_neg.length;i++){
            vals[i]=0.;
        }
        
        if(this.ut_negEFTable!=0){
            System.arraycopy(ut_neg, 0, vals, 0, this.ut_neg.length);
        }else{
            if(this.isMain()){
                L1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getLength();
                if(this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID()==this.id){
                    n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                }else{
                    n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                }
                
                L2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getLength();
                if(this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID()==this.id){
                    n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                }else{
                    n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                }
            }else{
                L1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getLength();
                if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getTwin().getID()==this.id){
                    n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                }else{
                    n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                }
                
                L2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getLength();
                if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getTwin().getID()==this.id){
                    n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                }else{
                    n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                }
            }
            
            for(int i=0;i<un.length;i++){
                if(this.isMain()){
                    vals[i]=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getNodes().get(n1).getut_neg()[i]*L2/(L1+L2)+
                        this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getNodes().get(n2).getut_neg()[i]*L1/(L1+L2);
                }else{
                    vals[i]=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getNodes().get(n1).getTwin().getut_neg()[i]*L2/(L1+L2)+
                        this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getNodes().get(n2).getTwin().getut_neg()[i]*L1/(L1+L2);
                }
                
            }
        }
        return vals;
    }

    public double[] getun(){
        double[] vals;
        double L1,L2;
        int n1,n2;
        vals=new double[this.un.length];
        for(int i=0;i<un.length;i++){
            vals[i]=0.;
        }
        if(this.unEFTable!=0){
            System.arraycopy(un, 0, vals, 0, this.un.length);
        }else{
            boolean haselems=false;
            if(this.connectedIElements.size()>0)haselems=true;
            if(haselems){
                L1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getLength();
                if(this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID()==this.id){
                    n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                }else{
                    n1=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                }
                
                L2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getLength();
                if(this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID()==this.id){
                    n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                }else{
                    n2=this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                }
            }else{
                L1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getLength();
                if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getTwin().getID()==this.id){
                    n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(2).getID();
                }else{
                    n1=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getINodeHier(1).getID();
                }
                
                L2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getLength();
                if(this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getTwin().getID()==this.id){
                    n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(2).getID();
                }else{
                    n2=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getINodeHier(1).getID();
                }
            }
            
            for(int i=0;i<un.length;i++){
                if(haselems){
                    vals[i]=this.connectedIElements.get(this.getConnectedIElementsIds()[0]).getNodes().get(n1).getun()[i]*L2/(L1+L2)+
                        this.connectedIElements.get(this.getConnectedIElementsIds()[1]).getNodes().get(n2).getun()[i]*L1/(L1+L2);
                }else{
                    vals[i]=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getNodes().get(n1).getTwin().getun()[i]*L2/(L1+L2)+
                        this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getNodes().get(n2).getTwin().getun()[i]*L1/(L1+L2);
                    //vals[i]=this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[0]).getNodes().get(n1).getun()[i]*L2/(L1+L2)+
                    //    this.getTwin().connectedIElements.get(this.getTwin().getConnectedIElementsIds()[1]).getNodes().get(n2).getun()[i]*L1/(L1+L2);
                }
                
            }
        }
        return vals;
    }
    
    public double[] getun_pos(){
        return this.un_pos;
    }
    
    public double[] getun_neg(){
        return this.un_neg;
    }
    
    public double getuth(){
        return this.uth;
    }

    public double getunh(){
        return this.unh;
    }

    public double[] getz(int elemid){
        double[] ret = null;
        if(z!=null){
        int pos = this.connectedIElementsHierarchy.get(elemid);
        if(this.zEFTable.length==1) pos=1;
        ret = new double[z.length];
        for(int i=0;i<ret.length;i++){
            ret[i] = z[i][pos-1];
        }}
        return ret;
    }
    
    public double getdrz(int elemid, int itime,double kn,double kt,double aI){
        double retv=0.;
        if(z!=null){
            int pos = this.connectedIElementsHierarchy.get(elemid);
            if(this.zEFTable.length==1) pos=1;
            for(int i=1;i<=itime;i++){
                retv+=aI*(z[i][pos-1]-z[i-1][pos-1])
                        -getDamageDrivingForce(elemid,i-1,kn,kt)*(z[i][pos-1]-z[i-1][pos-1]);
            }
        }
        
        return retv;
    }
    
    public double getDamageDrivingForce(int elemid, int itime, double kn,double kt){
        double retv=0.;
        double un,ut,sc;
        if(z!=null){
            un=getun()[itime];
            ut=getut()[itime];
            sc=0.;
            if(s_pos!=null){
                sc=gets(elemid)[itime];
            }
            retv=0.5*(un*un*kn+(ut-sc)*(ut-sc)*kt);
        }
        return retv;
    }
    
    public double getdrs(int elemid, int itime,double kt,double kh,double sy){
        double retv=0.;
        for(int i=1;i<=itime;i++){
        double sc=gets(elemid)[i];
        double sp=gets(elemid)[i-1];
            retv+=sy*Math.abs(sc-sp)
                    -getSlipDrivingForce(elemid,i-1,kt,kh)*(sc-sp);
        }
        
        return retv;
    }
    
    public double getSlipDrivingForce(int elemid, int itime, double kt,double kh){
        double retv=0.;
        double z,ut,s;
        z=getz(elemid)[itime];
        s=gets(elemid)[itime];
        ut=getut()[itime];
        retv=z*kt*(s-ut)+kh*s;
        return retv;
    }

    public double[] gets(int elemid){
        double[] ret = null; 
        
        if(this.s_posEFTable.length==this.s_negEFTable.length){
            int pos = this.connectedIElementsHierarchy.get(elemid);
            if(this.s_posEFTable.length==1) pos=1;
            ret = new double[s_pos.length];
            ret[0]=s_pos[0][pos-1]-s_neg[0][pos-1];
            for(int i=1;i<ret.length;i++){
                if(InterfaceNode.PlasticIncrement){
                    ret[i] =ret[i-1] + s_pos[i][pos-1]-s_neg[i][pos-1];
                }else{
                    ret[i] = s_pos[i][pos-1]-s_neg[i][pos-1];
//                    if(getz(elemid)!=null){if(Math.abs(getz(elemid)[i])<=1.e-18){ret[i]=ret[i-1];}}
                    if(Math.abs(ret[i])<=1.e-16){
                        int k=i;
                        while(k>=1 && Math.abs(ret[i])<=1.e-16){
                            ret[i]= s_pos[k-1][pos-1]-s_neg[k-1][pos-1];
                            k-=1;
                        }
//                        ret[i] = s_pos[i-1][pos-1]-s_neg[i-1][pos-1];
                    }
                }
            }
        }else{
            System.err.println("Error in gets of InterfaceNode. STOP");
            System.exit(elemid);
        }
        return ret;
    }
    
    public double[] gets_pos(int elemid){
        int pos = this.connectedIElementsHierarchy.get(elemid);
        if(this.s_posEFTable.length==1) pos=1;
        double[] ret = new double[s_pos.length];
        for(int i=0;i<ret.length;i++){
            ret[i] = s_pos[i][pos-1];
        }
        return ret;
    }
    
    public double[] gets_neg(int elemid){
        int pos = this.connectedIElementsHierarchy.get(elemid);
        if(this.s_negEFTable.length==1) pos=1;
        double[] ret = new double[s_neg.length];
        for(int i=0;i<ret.length;i++){
            ret[i] = s_neg[i][pos-1];
        }
        return ret;
    }
    
    public double[] getp(int elemid){
        double[] ret = null; 
        
        if(this.p_posEFTable.length==this.p_negEFTable.length){
            int pos = this.connectedIElementsHierarchy.get(elemid);
            if(this.p_posEFTable.length==1) pos=1;
            ret = new double[p_pos.length];
            ret[0]=p_pos[0][pos-1]-p_neg[0][pos-1];
            for(int i=1;i<ret.length;i++){
                if(InterfaceNode.PlasticIncrement){ret[i] =ret[i-1] + p_pos[i][pos-1]-p_neg[i][pos-1];}
                else{ret[i] = p_pos[i][pos-1]-p_neg[i][pos-1];}
            }
        }else{
            System.err.println("Error in gets of InterfaceNode. STOP");
            System.exit(elemid);
        }
        
        
        return ret;
    }
    
    public double[] getp_pos(int elemid){
        int pos = this.connectedIElementsHierarchy.get(elemid);
        if(this.p_posEFTable.length==1) pos=1;
        double[] ret = new double[p_pos.length];
        for(int i=0;i<ret.length;i++){
            ret[i] = p_pos[i][pos-1];
        }
        return ret;
    }
    
    public double[] getp_neg(int elemid){
        int pos = this.connectedIElementsHierarchy.get(elemid);
        if(this.p_negEFTable.length==1) pos=1;
        double[] ret = new double[p_neg.length];
        for(int i=0;i<ret.length;i++){
            ret[i] = p_neg[i][pos-1];
        }
        return ret;
    }
    
    public void print(){
        DecimalFormat Places = new DecimalFormat("0.000000000000000000");
        String vv;
        System.out.print(id+" "+this.mainNode+" "+this.inequation+" ");
        for(int i=0; i<this.coordinates.length; i++){
            vv= Places.format(this.coordinates[i]);
            System.out.print(vv+" ");
        }
        if(this.BaseNode!=null){System.out.print(this.BaseNode.getID()+" ");}
        else{System.out.print(0+" ");}
        if(this.TwinPoint!=null){System.out.print(this.TwinPoint.getID()+" ");}
        else{System.out.print(0+" ");}
        if(this.onElement!=null){System.out.print(this.onElement.getID()+" ");}
        else{System.out.print(0+" ");}
        System.out.print(this.getGap()+" ");
        if(normal!=null){
            for(int i=0; i<this.normal.length; i++){
                vv= Places.format(this.normal[i]);
                System.out.print(vv+" ");
            }
        }
        if(this.unEFTable!=0)System.out.print(" un freedoms: "+unEFTable);
        if(this.un_posEFTable!=0)System.out.print(" un{+} freedoms: "+un_posEFTable);
        if(this.un_negEFTable!=0)System.out.print(" un{-} freedoms: "+un_negEFTable);
        if(this.utEFTable!=0)System.out.print(" ut freedoms: "+utEFTable);
        if(this.ut_posEFTable!=0)System.out.print(" ut{+} freedoms: "+ut_posEFTable);
        if(this.ut_negEFTable!=0)System.out.print(" ut{-} freedoms: "+ut_negEFTable);
        if(zEFTable!=null){System.out.print(" z freedoms: ");for(int i=0; i<this.zEFTable.length;i++)System.out.print(zEFTable[i]+" ");}
        if(sEFTable!=null){System.out.print(" s freedoms: ");
        for(int i=0; i<this.sEFTable.length;i++)System.out.print(sEFTable[i]+" ");}
        if(pEFTable!=null){System.out.print(" p freedoms: ");
        for(int i=0; i<this.pEFTable.length;i++)System.out.print(pEFTable[i]+" ");}
        if(s_posEFTable!=null){System.out.print(" s_pos freedoms: ");
        for(int i=0; i<this.s_posEFTable.length;i++)System.out.print(s_posEFTable[i]+" ");}
        if(s_negEFTable!=null){System.out.print(" s_neg freedoms: ");
        for(int i=0; i<this.s_negEFTable.length;i++)System.out.print(s_negEFTable[i]+" ");}
        if(p_posEFTable!=null){System.out.print(" p_pos freedoms: ");
        for(int i=0; i<this.p_posEFTable.length;i++)System.out.print(p_posEFTable[i]+" ");}
        if(p_negEFTable!=null){System.out.print(" p_neg freedoms: ");
        for(int i=0; i<this.p_negEFTable.length;i++)System.out.print(p_negEFTable[i]+" ");}
        System.out.println();
    }
    
    public void printResponse(){
        this.print();
        System.out.print("step"+" ");
        System.out.print("un"+" ");
        if(this.utEFTable!=0)System.out.print("ut"+" ");
        if(this.BaseNode!=null){
            System.out.print("base("+this.getBase().getID()+")_ux"+" ");
            System.out.print("base("+this.getBase().getID()+")_uy"+" ");
            System.out.print("base("+this.getBase().getID()+")_dux"+" ");
            System.out.print("base("+this.getBase().getID()+")_duy"+" ");
        }
        if(this.getTwin().getBase()!=null){
            System.out.print("tbase("+this.getTwin().getBase().getID()+")_ux"+" ");
            System.out.print("tbase("+this.getTwin().getBase().getID()+")_uy"+" ");
            System.out.print("tbase("+this.getTwin().getBase().getID()+")_dux"+" ");
            System.out.print("tbase("+this.getTwin().getBase().getID()+")_duy"+" ");
        }
        int[] elemsID;
        if(this.BaseNode!=null){
            elemsID=this.getBase().getConnectedElementsIds();
            for(int j=0;j<elemsID.length;j++){
                System.out.print("base_el_"+this.getBase().getConnectedElement(elemsID[j]).getID()+"_px ");
                System.out.print("base_el_"+this.getBase().getConnectedElement(elemsID[j]).getID()+"_py ");
            }
        }
        if(this.getTwin().getBase()!=null){
            elemsID=this.getTwin().getBase().getConnectedElementsIds();
            for(int j=0;j<elemsID.length;j++){
                System.out.print("base_el_"+this.getTwin().getBase().getConnectedElement(elemsID[j]).getID()+"_px ");
                System.out.print("base_el_"+this.getTwin().getBase().getConnectedElement(elemsID[j]).getID()+"_py ");
            }
        }
        
        System.out.println();
        
        
        for(int i=0;i<un.length;i++){
            System.out.print(i+" ");
            System.out.print(un[i]+" ");
            if(this.utEFTable!=0)System.out.print(ut[i]+" ");
            if(this.BaseNode!=null){
                System.out.print(this.getBase().getu()[0][i][3]+" ");
                System.out.print(this.getBase().getu()[1][i][3]+" ");
                System.out.print((this.getBase().getCoordinates()[0]+this.getBase().getu()[0][i][3])+" ");
                System.out.print((this.getBase().getCoordinates()[1]+this.getBase().getu()[1][i][3])+" ");
            }
            if(this.getTwin().getBase()!=null){
                System.out.print(this.getTwin().getBase().getu()[0][i][3]+" ");
                System.out.print(this.getTwin().getBase().getu()[1][i][3]+" ");
                System.out.print((this.getTwin().getBase().getCoordinates()[0]+this.getTwin().getBase().getu()[0][i][3])+" ");
                System.out.print((this.getTwin().getBase().getCoordinates()[1]+this.getTwin().getBase().getu()[1][i][3])+" ");
            }
            
            if(this.BaseNode!=null){
                elemsID=this.getBase().getConnectedElementsIds();
                for(int j=0;j<elemsID.length;j++){
                    System.out.print(this.getBase().getp(this.getBase().getConnectedElement(elemsID[j]), 2, 3)[0][i]+" ");
                    System.out.print(this.getBase().getp(this.getBase().getConnectedElement(elemsID[j]), 2, 3)[1][i]+" ");
                }
            }
            if(this.getTwin().getBase()!=null){
                elemsID=this.getTwin().getBase().getConnectedElementsIds();
                for(int j=0;j<elemsID.length;j++){
                    System.out.print(this.getTwin().getBase().getp(this.getTwin().getBase().getConnectedElement(elemsID[j]), 2, 3)[0][i]+" ");
                    System.out.print(this.getTwin().getBase().getp(this.getTwin().getBase().getConnectedElement(elemsID[j]), 2, 3)[1][i]+" ");
                }
            }
            
            System.out.println();
        }
    }

    public void printunDOFS(){
        System.out.print(id+" ");
        if(this.unEFTable!=0) System.out.print(unEFTable);
        System.out.println();
    }

    public void printzDOFS(){
        System.out.print(id+" ");
        for(int i=0;i<zEFTable.length;i++){
            System.out.print(zEFTable[i]+" ");
        }
        System.out.println();
    }

    public int  getNumINodes(){return InterfaceNode.numNodes;}

    public void setTwinINode(InterfaceNode TwinNode){this.TwinPoint = TwinNode;}

    /*public void setNormal(double[] nor){
        this.normal = new double[nor.length];
        for(int i=0; i<nor.length;i++){
            this.normal[i]=nor[i];
        }
    }*/

    public void setNormal(double[] nor, double coef){
        this.normal = new double[nor.length];
        double arg=0.;
        for(int i=0; i<nor.length;i++){
            arg+=nor[i]*nor[i];
        }
        arg = Math.sqrt(arg);
        for(int i=0; i<normal.length;i++){
            this.normal[i]=coef*nor[i]/arg;
        }
//        if(this.getBaseID()==8 ||  this.getBaseID()==12 ||  this.getBaseID()==14){
//            System.err.println(this.getBaseID()+" "+normal[0]+" "+normal[1]);
//        }
    }

    public double[] getNormal(){
        return this.normal;
    }

    public double[] getGapVector(){
        double[] vec;
        vec = new double[this.coordinates.length];
        if(this.TwinPoint==null){
            for(int i=0; i<this.coordinates.length; i++){
                //vec[i] = this.coordinates[i];
                vec[i] = 0.0;
            }
            vec[1] = this.coordinates[1];
        }else{
            for(int i=0; i<this.coordinates.length; i++){
                vec[i] = this.coordinates[i]-this.TwinPoint.getCoordinates()[i];
            }

        }
        return vec;
    }

    public double[] getIntermediateCoords(){
        double[] vec;
        vec = new double[this.coordinates.length];
        if(this.TwinPoint==null){
            for(int i=0; i<this.coordinates.length; i++){
                vec[i] = this.coordinates[i];
            }
        }else{
            for(int i=0; i<this.coordinates.length; i++){
                vec[i] = (this.coordinates[i]+this.TwinPoint.getCoordinates()[i])/2.;
            }

        }
        return vec;
    }

    public double getGap(){
        double val=0;
        for(int i=0; i<this.getGapVector().length;i++){
            val +=this.getGapVector()[i]*this.getGapVector()[i];
        }
        return Math.sqrt(val);
    }

    public int getTwinID(){
        int tid=0;
        if(this.TwinPoint!=null)tid=this.TwinPoint.getID();
        return tid;
    }

    public InterfaceNode getTwin(){
        return TwinPoint;
    }

    public int getBaseID(){
        int tid=0;
        if(this.BaseNode!=null)tid=this.BaseNode.getID();
        return tid;
    }

    public Node getBase(){
        return BaseNode;
    }

    public void putIElement(InterfaceElement theElem){
        ++numOfNodeIElements;
        this.connectedIElements.put(theElem.getID(), theElem);
        this.connectedIElementsHierarchy.put(theElem.getID(),numOfNodeIElements);
    }

    public int getNumOfConnectedIElements(){
        return this.connectedIElements.size();
    }

    public void setMain(){this.mainNode=true;}

    public boolean isMain(){return this.mainNode;}

    public double getDist(InterfaceNode anotherNode){
        double dist=0.;
        for(int i=0;i<this.coordinates.length;i++){
            dist+=(this.coordinates[i]-anotherNode.coordinates[i])*(this.coordinates[i]-anotherNode.coordinates[i]);
        }
        dist=Math.sqrt(dist);
        return dist;
    }

    public void setunEFTable(int n){
        this.unEFTable=n;
    }
    
    public void setun_posEFTable(int n){
        this.un_posEFTable=n;
    }
    
    public void setun_negEFTable(int n){
        this.un_negEFTable=n;
    }

    public void setutEFTable(int n){
        this.utEFTable=n;
    }
    
    public void setut_posEFTable(int n){
        this.ut_posEFTable=n;
    }
    
    public void setut_negEFTable(int n){
        this.ut_negEFTable=n;
    }

    public void setzEFTable(int ndof, int elemid){
        int pos = connectedIElementsHierarchy.get(elemid);
        zEFTable[pos-1]=ndof;
    }

    public void setzEFTable(int ndof){
        zEFTable[0]=ndof;
    }

//    public void setsEFTable(int ndof, int elemid){
//        int pos = connectedIElementsHierarchy.get(elemid);
//        sEFTable[pos-1]=ndof;
//    }
    
    public void sets_posEFTable(int ndof, int elemid){
        int pos = connectedIElementsHierarchy.get(elemid);
        s_posEFTable[pos-1]=ndof;
    }
    
    public void sets_negEFTable(int ndof, int elemid){
        int pos = connectedIElementsHierarchy.get(elemid);
        s_negEFTable[pos-1]=ndof;
    }

//    public void setsEFTable(int ndof){
//        sEFTable[0]=ndof;
//    }
    
    public void sets_posEFTable(int ndof){
        s_posEFTable[0]=ndof;
    }
    
    public void sets_negEFTable(int ndof){
        s_negEFTable[0]=ndof;
    }
    
//    public void setpEFTable(int ndof, int elemid){
//        int pos = connectedIElementsHierarchy.get(elemid);
//        pEFTable[pos-1]=ndof;
//    }
    
    public void setp_posEFTable(int ndof, int elemid){
        int pos = connectedIElementsHierarchy.get(elemid);
        p_posEFTable[pos-1]=ndof;
    }
    
    public void setp_negEFTable(int ndof, int elemid){
        int pos = connectedIElementsHierarchy.get(elemid);
        p_negEFTable[pos-1]=ndof;
    }

//    public void setpEFTable(int ndof){
//        pEFTable[0]=ndof;
//    }
    
    public void setp_posEFTable(int ndof){
        p_posEFTable[0]=ndof;
    }
    
    public void setp_negEFTable(int ndof){
        p_negEFTable[0]=ndof;
    }

    public int getutEFTable(){
        return this.utEFTable;
    }

    public int getunEFTable(){
        return this.unEFTable;
    }
    
    public int getun_posEFTable(){
        return this.un_posEFTable;
    }
    
    public int getun_negEFTable(){
        return this.un_negEFTable;
    }
    
    public int getut_posEFTable(){
        return this.ut_posEFTable;
    }
    
    public int getut_negEFTable(){
        return this.ut_negEFTable;
    }

    public int getzEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.zEFTable.length==1) pos=1;
        return this.zEFTable[pos-1];
    }

    public int[] getzEFTable(){
        return this.zEFTable;
    }

    public int[] getsEFTable(){
        return this.sEFTable;
    }
    
    public int[] gets_posEFTable(){
        return this.s_posEFTable;
    }
    
    public int[] gets_negEFTable(){
        return this.s_negEFTable;
    }

    public int getsEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.sEFTable.length==1) pos=1;
        return this.sEFTable[pos-1];
    }
    
    public int gets_posEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.s_posEFTable.length==1) pos=1;
        return this.s_posEFTable[pos-1];
    }
    
    public int gets_negEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.s_negEFTable.length==1) pos=1;
        return this.s_negEFTable[pos-1];
    }
    
    public int[] getpEFTable(){
        return this.pEFTable;
    }
    
    public int[] getp_posEFTable(){
        return this.p_posEFTable;
    }
    
    public int[] getp_negEFTable(){
        return this.p_negEFTable;
    }

    public int getpEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.pEFTable.length==1) pos=1;
        return this.pEFTable[pos-1];
    }
    
    public int getp_posEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.p_posEFTable.length==1) pos=1;
        return this.p_posEFTable[pos-1];
    }
    
    public int getp_negEFTable(int elemID){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.p_negEFTable.length==1) pos=1;
        return this.p_negEFTable[pos-1];
    }
    
    public Element getOnElement(){return this.onElement;}

    public void setInEquation(){
        if(this.BaseNode!=null)this.inequation=true;
    }

    public boolean AskIfInEquation(){return this.inequation;}

    public InterfaceElement getSomeElement(){
        InterfaceElement relem = null;
        if(this.getNumOfConnectedIElements()>0) relem = this.connectedIElements.values().iterator().next();
        return relem;
    }

    public void updates_un(int step, double val){
        this.un[step]=val;
        if(un[step]>=0){
            un_neg[step]=0.0; un_pos[step]=un[step];
        }else{
            un_neg[step]=un[step]; un_pos[step]=0.0;
        }
    }
    
    public void updates_un_pos(int step, double val){
        this.un_pos[step]=val;
    }
    
    public void updates_un_neg(int step, double val){
        this.un_neg[step]=val;
    }
    
    public void updates_un_split(int step){
        this.un[step]=this.un_neg[step]+this.un_pos[step];
        if(un[step]>=0){
            un_neg[step]=0.0; un_pos[step]=un[step];
        }else{
            un_neg[step]=un[step]; un_pos[step]=0.0;
        }
    }
    

    public void updates_ut(int step, double val){
        this.ut[step]=val;
    }
    
    public void updates_ut_pos(int step, double val){
        this.ut_pos[step]=val;
    }
    
    public void updates_ut_neg(int step, double val){
        this.ut_neg[step]=val;
    }
    
    public void updates_unh(double val){
        this.unh=val;
    }

    public void updates_uth(double val){
        this.uth=val;
    }

    public void updates_z(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.zEFTable.length==1) pos=1;
        this.z[step][pos-1]=val;
    }

    public void updates_s(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.sEFTable.length==1) pos=1;
        this.s[step][pos-1]=val;
    }
    
    public void updates_s_pos(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.s_posEFTable.length==1) pos=1;
        this.s_pos[step][pos-1]=val;
    }
    
    public void updates_s_neg(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.s_negEFTable.length==1) pos=1;
        this.s_neg[step][pos-1]=val;
    }
    
    public void updates_p(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.pEFTable.length==1) pos=1;
        this.p[step][pos-1]=val;
    }
    
    public void updates_p_pos(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.p_posEFTable.length==1) pos=1;
        this.p_pos[step][pos-1]=val;
    }
    
    public void updates_p_neg(int step, int elemID, double val){
        int pos = this.connectedIElementsHierarchy.get(elemID);
        if(this.p_negEFTable.length==1) pos=1;
        this.p_neg[step][pos-1]=val;
    }

    public int[] getConnectedIElementsIds(){
        int[] Ids = new int[this.connectedIElements.size()];
        int i=0;
        for(Iterator<InterfaceElement> it=this.connectedIElements.values().iterator(); it.hasNext();){
            InterfaceElement theElement = it.next();
            Ids[i]=theElement.getID();
            i++;
        }
        return Ids;
    }

    public int getHierOfElement(int elemID){
        return this.connectedIElementsHierarchy.get(elemID);
    }
    
    public boolean AskIfCollapsed(int step){
        boolean ans=true;
        if(z!=null){
            ans=false;
            if(step>0){
                for(int i=0;i<this.z[step].length;i++){
                    if(z[step][i]<1.e-10)ans=true;
                }
            }else{
                for(int i=0;i<this.z[0].length;i++){
                    if(z[0][i]<1.e-10)ans=true;
                }
            }
        }
        return ans;
    }

    void setZeroResponse(int itime) {
        if(un!=null)un[itime]=0;
        if(un_pos!=null)un_pos[itime]=0;
        if(un_neg!=null)un_neg[itime]=0;
        if(ut!=null)ut[itime]=0;
        if(ut_pos!=null)ut_pos[itime]=0;
        if(ut_neg!=null)ut_neg[itime]=0;
        if(s!=null)for(int i=0;i<s[0].length;i++){s[itime][i]=1.0;}
        if(p!=null)for(int i=0;i<p[0].length;i++){p[itime][i]=1.0;}
    }
    
    public void setPlasticIncrement(boolean what){
        InterfaceNode.PlasticIncrement=what;
    }
    
    public double[] getun_(ViscousMaterial vMaterial, double tau){
        double[] vals;
        double Phi1 = vMaterial.getPhi_1(tau);
        double Phi2 = vMaterial.getPhi_2(tau);
        double Phi3 = vMaterial.getPhi_3(tau);
        double nu=vMaterial.getDispRelaxationTime_1();
        vals=new double[this.un.length];
        for(int i=0;i<un.length;i++){
            vals[i]=0.;
        }
        for(int i=1;i<un.length;i++){
//            if(i>1){
//                vals[i]=getun()[i]/Phi1-vals[i-1]*Phi2/Phi1-vals[i-2]*Phi3/Phi1;
//            }else{
//                vals[i]=getun()[i]/Phi1-vals[i-1]*Phi2/Phi1;
//            }
            if(i>1)vals[i]=(1./(tau+nu))*(tau*this.getun()[i]+nu*vals[i-1]);
        }
        return vals;
    }
    
    public double[] getut_(ViscousMaterial vMaterial, double tau){
        double[] vals;
        double Phi1 = vMaterial.getPhi_1(tau);
        double Phi2 = vMaterial.getPhi_2(tau);
        double Phi3 = vMaterial.getPhi_3(tau);
        double nu=vMaterial.getDispRelaxationTime_1();
        vals=new double[this.un.length];
        for(int i=0;i<ut.length;i++){
            vals[i]=0.;
        }
        for(int i=1;i<ut.length;i++){
//            if(i>1){
//                vals[i]=getut()[i]/Phi1-vals[i-1]*Phi2/Phi1-vals[i-2]*Phi3/Phi1;
//            }else{
//                vals[i]=getut()[i]/Phi1-vals[i-1]*Phi2/Phi1;
//            }
            if(i>1)vals[i]=(1./(tau+nu))*(tau*this.getut()[i]+nu*vals[i-1]);
        }
        return vals;
    }
    
    public double[] getvn_(ViscousMaterial vMaterial, double tau){
        double[] vals;
        double Phi1 = vMaterial.getPhi_1(tau);
        double Phi2 = vMaterial.getPhi_2(tau);
        double Phi3 = vMaterial.getPhi_3(tau);
        double nu=vMaterial.getDispRelaxationTime_1();
        vals=new double[this.un.length];
        for(int i=0;i<un.length;i++){
            vals[i]=0.;
        }
        for(int i=1;i<un.length;i++){
//            if(i>1){
//                vals[i]=Phi1*getun()[i]+Phi2*getun()[i-1]+Phi3*getun()[i-2];
//            }else{
//                vals[i]=Phi1*getun()[i]+Phi2*getun()[i-1];
//            }
            if(i>1)vals[i]=(1./(tau+nu))*(this.getun()[i]-this.getun()[i-1]);
        }
        return vals;
    }
    
    public double[] getvt_(ViscousMaterial vMaterial, double tau){
        double[] vals;
        double Phi1 = vMaterial.getPhi_1(tau);
        double Phi2 = vMaterial.getPhi_2(tau);
        double Phi3 = vMaterial.getPhi_3(tau);
        double nu=vMaterial.getDispRelaxationTime_1();
        vals=new double[this.un.length];
        for(int i=0;i<ut.length;i++){
            vals[i]=0.;
        }
        for(int i=1;i<ut.length;i++){
//            if(i>1){
//                vals[i]=Phi1*getut()[i]+Phi2*getut()[i-1]+Phi3*getut()[i-2];
//            }else{
//                vals[i]=Phi1*getut()[i]+Phi2*getut()[i-1];
//            }
            if(i>1)vals[i]=(1./(tau+nu))*(this.getut()[i]-this.getut()[i-1]);
        }
        return vals;
    }
}
