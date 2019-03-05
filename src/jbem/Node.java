/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author pchr
 */
public class Node extends geom.Point{
    private static int numNodes=0;
    private int numOfNodeElements=0;
    private double[][][] u;
    private double[][] u_eig;
    private double[][][] v;
    private double[][][] p;
    private double[] uh;
    private double[] ph;
    private Map<Integer,Element> connectedElements = new HashMap<Integer,Element>();
    private Map<Integer,Integer> connectedElementsHierarchy = new HashMap<Integer,Integer>();
    private int[] pEFTable;
    private int[] uEFTable;
    private int[] vEFTable;
    
    // auxiliary displacement and tractions for Kelvin-Voigt materials
    private double[][][] u_aux;
    private double[][][] p_aux;
    
    // constructor
    public Node(){}
    
    public Node(int id, double x){
        if(id>0){
            ++numNodes;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            this.coordinates[0]=x;
            this.coordinates[1]=0.0;
            this.coordinates[2]=0.0;
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public Node(int id, double x, double y){
        if(id>0){
            ++numNodes;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            this.coordinates[0]=x;
            this.coordinates[1]=y;
            this.coordinates[2]=0.0;
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public Node(int id, double x, double y, double z){
        if(id>0){
            ++numNodes;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            this.coordinates[0]=x;
            this.coordinates[1]=y;
            this.coordinates[2]=z;
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    public Node(int id, double[] coords){
        if(id>0){
            ++numNodes;
            this.id=id;
            //this.coordinates = new double[coords.length];
            this.coordinates = new double[3];
            for(int i=0; i<3; i++){
                if(i<coords.length){this.coordinates[i]=coords[i];}
                else{this.coordinates[i]=0.;}
            }
        }else{
            System.err.println("invalid id= "+id+" node did not constructed !!");
        }
    }
    
    // methods
    public double[][][] getu(){
        return this.u;
    }
    
    public double[][][] getu_aux(){
        return this.u_aux;
    }

    public double[][] getu_eig(){
        return this.u_eig;
    }
    
    public double[][][] getv(){
        return this.v;
    }
    
    public double[][][] getp(){
        return this.p;
    }
    
    public double[][][] getp_aux(){
        return this.p_aux;
    }

    public double[] getuh(){
        return this.uh;
    }

    public double[] getph(){
        return this.ph;
    }
    
    public double[][] getp(Element theElement, int ndof){
        double[][] elemp;
        //int[] efelem=this.getpEFTable(theElement, ndof);
        int j=this.p[0].length;
        elemp = new double[ndof][j];

        int ElemHierarchy=this.connectedElementsHierarchy.get(theElement.getID());
        for(int jj=0; jj<j; jj++){
            for(int ii=0; ii<ndof; ii++){
                elemp[ii][jj]=this.p[(ElemHierarchy-1)*ndof+ii][jj][0];
            }
        }
        
        return elemp;
    }
    
    public double[][] getp(Element theElement, int ndof, int state){
        double[][] elemp;
        //int[] efelem=this.getpEFTable(theElement, ndof);
        int j=this.p[0].length;
        elemp = new double[ndof][j];

        int ElemHierarchy=this.connectedElementsHierarchy.get(theElement.getID());
        for(int jj=0; jj<j; jj++){
            for(int ii=0; ii<ndof; ii++){
                elemp[ii][jj]=this.p[(ElemHierarchy-1)*ndof+ii][jj][state];
            }
        }
        
        return elemp;
    }
    
    public double[][] getp_aux(Element theElement, int ndof){
        double[][] elemp;
        //int[] efelem=this.getpEFTable(theElement, ndof);
        int j=this.p_aux[0].length;
        elemp = new double[ndof][j];

        int ElemHierarchy=this.connectedElementsHierarchy.get(theElement.getID());
        for(int jj=0; jj<j; jj++){
            for(int ii=0; ii<ndof; ii++){
                elemp[ii][jj]=this.p_aux[(ElemHierarchy-1)*ndof+ii][jj][0];
            }
        }
        
        return elemp;
    }
    
    public double[][] getp_aux(Element theElement, int ndof, int state){
        double[][] elemp;
        //int[] efelem=this.getpEFTable(theElement, ndof);
        int j=this.p_aux[0].length;
        elemp = new double[ndof][j];

        int ElemHierarchy=this.connectedElementsHierarchy.get(theElement.getID());
        for(int jj=0; jj<j; jj++){
            for(int ii=0; ii<ndof; ii++){
                elemp[ii][jj]=this.p_aux[(ElemHierarchy-1)*ndof+ii][jj][state];
            }
        }
        
        return elemp;
    }

    public double[] getph(Element theElement, int ndof){
        double[] elemp;
        int[] efelem=this.getpEFTable(theElement, ndof);
        int j=this.ph.length;
        elemp = new double[ndof];

        int ElemHierarchy=this.connectedElementsHierarchy.get(theElement.getID());
        for(int ii=0; ii<ndof; ii++){
            elemp[ii]=this.ph[(ElemHierarchy-1)*ndof+ii];
        }

        return elemp;
    }
    
    public void putElement(Element theElem){
        ++numOfNodeElements;
        this.connectedElements.put(theElem.getID(), theElem);
        this.connectedElementsHierarchy.put(theElem.getID(),numOfNodeElements);
    }
    
    public void init_u(int n, int steps){
        this.u=new double[n][steps][1];
        for(int i=0;i<n;i++)u[i][0][1]=0.0;
    }
    
    public void init_u(int n, int steps, int nstates){
        this.u=new double[n][steps][nstates];
        for(int i=0;i<n;i++){
            for(int j=0;j<nstates;j++){
                u[i][0][j]=0.0;
            }
        }
    }
    
    public void init_u_aux(int n, int steps){
        this.u_aux=new double[n][steps][1];
        for(int i=0;i<n;i++)u_aux[i][0][1]=0.0;
    }
    
    public void init_u_aux(int n, int steps, int nstates){
        this.u_aux=new double[n][steps][nstates];
        for(int i=0;i<n;i++){
            for(int j=0;j<nstates;j++){
                u_aux[i][0][j]=0.0;
            }
        }
    }

    public void init_uh(int n, int steps){
        this.uh=new double[n];
    }
    
    public void init_v(int n, int steps){
        this.v=new double[n][steps][1];
        for(int i=0;i<n;i++)v[i][0][1]=0.0;
    }
    
    public void init_v(int n, int steps, int nstates){
        this.v=new double[n][steps][nstates];
        for(int i=0;i<n;i++){
            for(int j=0;j<nstates;j++){
                v[i][0][j]=0.0;
            }
        }
    }
    
    public void init_p(int n, int steps){
        int nn=this.connectedElements.size()*n;
        this.p=new double[nn][steps][1];
        for(int i=0;i<nn;i++)p[i][0][1]=0.0;
    }
    
    public void init_p(int n, int steps, int nstates){
        int nn=this.connectedElements.size()*n;
        this.p=new double[nn][steps][nstates];
        for(int i=0;i<nn;i++){
            for(int j=0;j<nstates;j++){
                p[i][0][j]=0.0;
            }
        }
    }
    
    public void init_p_aux(int n, int steps){
        int nn=this.connectedElements.size()*n;
        this.p_aux=new double[nn][steps][1];
        for(int i=0;i<nn;i++)p_aux[i][0][1]=0.0;
    }
    
    public void init_p_aux(int n, int steps, int nstates){
        int nn=this.connectedElements.size()*n;
        this.p_aux=new double[nn][steps][nstates];
        for(int i=0;i<nn;i++){
            for(int j=0;j<nstates;j++){
                p_aux[i][0][j]=0.0;
            }
        }
    }

    public void init_ph(int n, int steps){
        int nn=this.connectedElements.size()*n;
        this.ph=new double[nn];
    }


    public void init_u_eig(int numEIGS){
        this.u_eig=new double[this.u.length][numEIGS];
    }

    
    public void setpEFTable(int n, int size){
        this.pEFTable=new int[size];
        for(int i=0; i<size; i++){
            this.pEFTable[i]=n+i;
        }
    }
    
    public void setuEFTable(int n, int size){
        this.uEFTable=new int[size];
        for(int i=0; i<size; i++){
            this.uEFTable[i]=n+i;
        }
    }
    public void setvEFTable(int n, int size){
        this.vEFTable=new int[size];
        for(int i=0; i<size; i++){
            this.vEFTable[i]=n+i;
        }
    }
    
    public int getNumOfConnectedElements(){
        return this.connectedElements.size();
    }
    
    public int[] getuEFTable(){
        return this.uEFTable;
    }
    
    public int[] getvEFTable(){
        return this.vEFTable;
    }
    
    public int[] getpEFTable(){
        return this.pEFTable;
    }
    
    /**
     * ndof is the num of "traction" dofs per node,
     * this information is given by fundamental solution class
     */
    public int[] getpEFTable(int ElemHierarchy, int ndof){
        
        int[] eft = new int[ndof];
        for(int i=0; i<ndof; i++){
            eft[i]=this.pEFTable[(ElemHierarchy-1)*ndof+i];
        }
        return eft;
    }
    
    public int[] getpEFTable(Element theElem, int ndof){
        int[] eft = new int[ndof];
        int ElemHierarchy=this.connectedElementsHierarchy.get(theElem.getID());
        for(int i=0; i<ndof; i++){
            eft[i]=this.pEFTable[(ElemHierarchy-1)*ndof+i];
        }
        return eft;
    }
    
    public void setp(int dof, int step, double val){
        this.p[dof][step][0]=val;
    }
    
    public void setp(int dof, int step, double val, int state){
        this.p[dof][step][state]=val;
    }
    
    public void setu(int dof, int step, double val){
        this.u[dof][step][0]=val;
    }
    
    public void setu(int dof, int step, double val, int state){
        this.u[dof][step][state]=val;
    }
    
    public void setp_aux(int dof, int step, double val){
        this.p_aux[dof][step][0]=val;
    }
    
    public void setp_aux(int dof, int step, double val, int state){
        this.p_aux[dof][step][state]=val;
    }
    
    public void setu_aux(int dof, int step, double val){
        this.u_aux[dof][step][0]=val;
    }
    
    public void setu_aux(int dof, int step, double val, int state){
        this.u_aux[dof][step][state]=val;
    }

    public void setu_eig(int dof, int eig, double val){
        this.u_eig[dof][eig]=val;
    }

    public void setph(int dof, double val){
        this.ph[dof]=val;
    }

    public void setuh(int dof, double val){
        this.uh[dof]=val;
    }
    
    public void setv(int dof, int step, double val){
        this.v[dof][step][0]=val;
    }
    
    public void setv(int dof, int step, double val, int state){
        this.v[dof][step][state]=val;
    }
    
    public void print(){
        DecimalFormat Places = new DecimalFormat("0.000000000000000000");
        String vv;
        System.out.print(id+" ");
        for(int i=0; i<this.coordinates.length; i++){
            vv= Places.format(this.coordinates[i]);
            System.out.print(vv+" ");
        }
        System.out.println();
    }
    
//    public double getDist(Node anotherNode){
//        double dist=0.;
//        for(int i=0;i<this.coordinates.length;i++){
//            dist+=(this.coordinates[i]-anotherNode.coordinates[i])*(this.coordinates[i]-anotherNode.coordinates[i]);
//        }
//        dist=Math.sqrt(dist);
//        return dist;
//    }
    
//    public double getDist(Point anotherNode){
//        double dist=0.;
//        for(int i=0;i<this.coordinates.length;i++){
//            dist+=(this.coordinates[i]-anotherNode.coordinates[i])*(this.coordinates[i]-anotherNode.coordinates[i]);
//        }
//        dist=Math.sqrt(dist);
//        return dist;
//    }
    
    
    
    public int[] getConnectedElementsIds(){
        int[] Ids = new int[this.connectedElements.size()];
        int i=0;
        for(Iterator<Element> it=this.connectedElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            Ids[i]=theElement.getID();
            i++;
        }
        return Ids;
    }

    public boolean isMain(){
        boolean request=false;
        for(Iterator<Element> it=this.connectedElements.values().iterator(); it.hasNext();){
            if(it.next().isMain())request=true;
        }
        return request;
    }
    
    public void setResponse(double coef, int fromWhichStep, int toWhichStep, int WhichState){
        if(this.u!=null){
            for(int i=0;i<u.length;i++)u[i][toWhichStep][WhichState]=coef*u[i][fromWhichStep][WhichState];
        }
        if(this.p!=null){
            for(int i=0;i<p.length;i++)p[i][toWhichStep][WhichState]=coef*p[i][fromWhichStep][WhichState];
        }
        if(this.v!=null){
            for(int i=0;i<v.length;i++)v[i][toWhichStep][WhichState]=coef*v[i][fromWhichStep][WhichState];
        }
    }
    
    public void setZeroResponse(int toWhichStep, int WhichState){
        if(this.u!=null){
            for(int i=0;i<u.length;i++)u[i][toWhichStep][WhichState]=0.0;
        }
        if(this.p!=null){
            for(int i=0;i<p.length;i++)p[i][toWhichStep][WhichState]=0.0;
        }
        if(this.v!=null){
            for(int i=0;i<v.length;i++)v[i][toWhichStep][WhichState]=0.0;
        }
    }
    
    public void setResponse_aux(double coef, int fromWhichStep, int toWhichStep, int WhichState){
        if(this.u_aux!=null){
            for(int i=0;i<u_aux.length;i++)u_aux[i][toWhichStep][WhichState]=coef*u_aux[i][fromWhichStep][WhichState];
        }
        if(this.p_aux!=null){
            for(int i=0;i<p_aux.length;i++)p_aux[i][toWhichStep][WhichState]=coef*p_aux[i][fromWhichStep][WhichState];
        }
    }
    
    public void setZeroResponse_aux(int toWhichStep, int WhichState){
        if(this.u!=null){
            for(int i=0;i<u_aux.length;i++)u_aux[i][toWhichStep][WhichState]=0.0;
        }
        if(this.p!=null){
            for(int i=0;i<p_aux.length;i++)p_aux[i][toWhichStep][WhichState]=0.0;
        }
    }
    
    
    public Element getConnectedElement(int elemID){
        return this.connectedElements.get(elemID);
    }
    
    public double[] getNormalMean(){
        double[] vec = new double[3];
        vec[0]=0.0;
        vec[1]=0.0;
        vec[2]=0.0;
        ELine elem;
        for(Iterator<Element> it=this.connectedElements.values().iterator(); it.hasNext();){
            elem = (ELine) it.next();
            vec[0]=vec[0]+elem.getNormal(this.id)[0];
            vec[1]=vec[1]+elem.getNormal(this.id)[0];
            vec[2]=vec[2]+elem.getNormal(this.id)[0];
        }
        vec[0]=vec[0]/this.connectedElements.size();
        vec[1]=vec[1]/this.connectedElements.size();
        vec[2]=vec[2]/this.connectedElements.size();
        return vec;
    }

    public void Aux2MainVariables(int itime,double tau, Material theViscousMaterial) {
        if(theViscousMaterial.getClass().getSimpleName().equalsIgnoreCase("ViscousMaterial")){
            double x0=((ViscousMaterial)theViscousMaterial).getDispRelaxationTime_0();
            double x1=((ViscousMaterial)theViscousMaterial).getDispRelaxationTime_1();
            double x2=((ViscousMaterial)theViscousMaterial).getDispRelaxationTime_2();
            double y0=((ViscousMaterial)theViscousMaterial).getStressRelaxationTime_0();
            double y1=((ViscousMaterial)theViscousMaterial).getStressRelaxationTime_1();
            double y2=((ViscousMaterial)theViscousMaterial).getStressRelaxationTime_2();

            for(int i=0;i<u_aux.length;i++){
                if(itime>1){
                    u[i][itime][3]=((y2+tau*y1+y0*tau*tau)*u_aux[i][itime][3]
                        +(2*x2+tau*x1)*u[i][itime-1][3]
                        -(x2)*u[i][itime-2][3])/(x2+x1*tau+x0*tau*tau);
                }else if(itime>0){
                    u[i][itime][3]=((y2+tau*y1+y0*tau*tau)*u_aux[i][itime][3]
                        +(2*x2+tau*x1)*u[i][itime-1][3])/(x2+x1*tau+x0*tau*tau);
                }else{
                    u[i][itime][3]=0.0;
//                    u[i][itime][3]=((y2+tau*y1+y0*tau*tau)*u_aux[i][itime][3])/(x2+x1*tau+x0*tau*tau);
                }
                u[i][itime][2]=u_aux[i][itime][2];
                u[i][itime][1]=u_aux[i][itime][1];
                u[i][itime][0]=u_aux[i][itime][0];
            }
            for(int i=0;i<p_aux.length;i++){
                if(itime>1){
                    p[i][itime][3]=p_aux[i][itime][3]
                            +((2*y2+y1*tau)/(y2+y1*tau+y0*tau*tau))*p[i][itime-1][3]
                            -((y2)/(y2+y1*tau+y0*tau*tau))*p[i][itime-2][3];
                }else if(itime>0){
                    p[i][itime][3]=p_aux[i][itime][3]
                            +((2*y2+y1*tau)/(y2+y1*tau+y0*tau*tau))*p[i][itime-1][3];
                }else{
                    p[i][itime][3]=0.0;
//                    p[i][itime][0]=p_aux[i][itime][0];
                }
                p[i][itime][2]=p_aux[i][itime][2];
                p[i][itime][1]=p_aux[i][itime][1];
                p[i][itime][0]=p_aux[i][itime][0];
            }
        }else{
            for(int i=0;i<u_aux.length;i++){
                u[i][itime][3]=u_aux[i][itime][3];
                u[i][itime][2]=u_aux[i][itime][2];
                u[i][itime][1]=u_aux[i][itime][1];
                u[i][itime][0]=u_aux[i][itime][0];
            }
            for(int i=0;i<p_aux.length;i++){
                p[i][itime][3]=p_aux[i][itime][3];
                p[i][itime][2]=p_aux[i][itime][2];
                p[i][itime][1]=p_aux[i][itime][1];
                p[i][itime][0]=p_aux[i][itime][0];
            }
        }
    }
    
    public void Aux2MainVariables_(int itime, double tau, Material theViscousMaterial) {
        if(theViscousMaterial.getClass().getSimpleName().equalsIgnoreCase("ViscousMaterial")){
            double x0=((ViscousMaterial)theViscousMaterial).getDispRelaxationTime_0();
            double x1=((ViscousMaterial)theViscousMaterial).getDispRelaxationTime_1();
            double x2=((ViscousMaterial)theViscousMaterial).getDispRelaxationTime_2();
            double y0=((ViscousMaterial)theViscousMaterial).getStressRelaxationTime_0();
            double y1=((ViscousMaterial)theViscousMaterial).getStressRelaxationTime_1();
            double y2=((ViscousMaterial)theViscousMaterial).getStressRelaxationTime_2();

            for(int i=0;i<u_aux.length;i++){
                if(itime>1){
                    u[i][itime][0]=((y2+tau*y1+y0*tau*tau)*u_aux[i][itime][0]
                        +(2*x2+tau*x1)*u[i][itime-1][0]
                        -(x2)*u[i][itime-2][0])/(x2+x1*tau+x0*tau*tau);
                }else if(itime>0){
                    u[i][itime][0]=((y2+tau*y1+y0*tau*tau)*u_aux[i][itime][0]
                        +(2*x2+tau*x1)*u[i][itime-1][0])/(x2+x1*tau+x0*tau*tau);
                }else{
//                    u[i][itime][0]=0.0;
//                    u[i][itime][0]=u_aux[i][itime][0];
                    u[i][itime][0]=((y2+tau*y1+y0*tau*tau)*u_aux[i][itime][0])/(x2+x1*tau+x0*tau*tau);
                }
            }
            for(int i=0;i<p_aux.length;i++){
                if(itime>1){
//                    p[i][itime][0]=((y2+tau*y1+y0*tau*tau)*p_aux[i][itime][0]
//                        +(2*x2+tau*x1)*p[i][itime-1][0]
//                        -(x2)*p[i][itime-2][0])/(x2+x1*tau+x0*tau*tau)
//                            +((2*y2+y1*tau)/(y1*tau+y0*tau*tau))*p_aux[i][itime-1][0]
//                            -((y2)/(y1*tau+y0*tau*tau))*p_aux[i][itime-2][0];
                    p[i][itime][0]=p_aux[i][itime][0]
                            +((2*y2+y1*tau)/(y2+y1*tau+y0*tau*tau))*p[i][itime-1][0]
                            -((y2)/(y2+y1*tau+y0*tau*tau))*p[i][itime-2][0];
                }else if(itime>0){
//                    p[i][itime][0]=((y2+tau*y1+y0*tau*tau)*p_aux[i][itime][0]
//                        +(2*x2+tau*x1)*p[i][itime-1][0])/(x2+x1*tau+x0*tau*tau)
//                            +((2*y2+y1*tau)/(y1*tau+y0*tau*tau))*p_aux[i][itime-1][0];
                    p[i][itime][0]=p_aux[i][itime][0]
                            +((2*y2+y1*tau)/(y2+y1*tau+y0*tau*tau))*p[i][itime-1][0];
                }else{
//                    p[i][itime][0]=0.0;
                    p[i][itime][0]=p_aux[i][itime][0];
//                    p[i][itime][0]=((y2+tau*y1+y0*tau*tau)*p_aux[i][itime][0])/(x2+x1*tau+x0*tau*tau);
                }
            }         
        }else{
            if(itime>=0){
                for(int i=0;i<u_aux.length;i++){
                    u[i][itime][0]=u_aux[i][itime][0];
                }
                for(int i=0;i<p_aux.length;i++){
                    p[i][itime][0]=p_aux[i][itime][0];
                }
            }else{
                for(int i=0;i<u_aux.length;i++){
                    u[i][itime][0]=0.;
    //                u[i][itime][0]=u_aux[i][itime][0];
                }
                for(int i=0;i<p_aux.length;i++){
                    p[i][itime][0]=0.0;
    //                p[i][itime][0]=p_aux[i][itime][0];
                }

            }
        }
    }
    
    public void setInitialDisplacement(double val, int dof, int state){
        this.u[dof][0][state]=val;
    }
    
    public void setInitialTraction(double val, int dof, int state){
        this.p[dof][0][state]=val;
    }

}
