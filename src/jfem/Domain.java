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
import java.awt.Color;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;
import java.util.Iterator;
import jmat.AbstractMatrix;

/**
 *
 * @author pchr
 */
public class Domain {
    private int id;
    private static int numberOfDomains = 0;
    private Map<Integer,Node> theNodes = new TreeMap<Integer,Node>();
    private Map<Integer,Element> theElements = new TreeMap<Integer,Element>();
    private Map<Integer,Material> theMaterials = new TreeMap<Integer,Material>();
    private Map<Integer,CrossSection> theCrossSections = new TreeMap<Integer,CrossSection>();
    private Map<Integer,InitialCondition> theInitialConditions = new TreeMap<Integer,InitialCondition>();
    private Map<Integer,LoadCase> theLoadCases = new TreeMap<Integer,LoadCase>();
    private Map<Integer,ConstraintElement> theConstraintElements = new TreeMap<Integer,ConstraintElement>();
    
    private double[] eigenvalues;
    
    private double maximum_X_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_X_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Y_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Y_coordinate=Double.POSITIVE_INFINITY;
    private double maximum_Z_coordinate=Double.NEGATIVE_INFINITY;
    private double minimum_Z_coordinate=Double.POSITIVE_INFINITY;
    
    public Color ColorNodesID=Color.LIGHT_GRAY;
    public Color ColorNodesUND=Color.orange;
    public Color ColorElemsUND=Color.blue;
    public Color ColorNodesDEF=Color.GREEN;
    public Color ColorElemsDEF=Color.RED;
    
    private int NumGridX=0;
    private int NumGridY=0;
    
    // constructors
    public Domain(){
        id = ++numberOfDomains;
    }
    
    public void clsNumberOfDomains(){numberOfDomains=0;}
    
    public static Domain dom(){
        return new Domain();
    }
    
    // Methods of Domain class
    public void putNode(Node aNode) {
        if(aNode.getCoords()[0]>maximum_X_coordinate)maximum_X_coordinate=aNode.getCoords()[0];
        if(aNode.getCoords()[0]<minimum_X_coordinate)minimum_X_coordinate=aNode.getCoords()[0];
        if(aNode.getCoords().length>1){
            if(aNode.getCoords()[1]>maximum_Y_coordinate)maximum_Y_coordinate=aNode.getCoords()[1];
            if(aNode.getCoords()[1]<minimum_Y_coordinate)minimum_Y_coordinate=aNode.getCoords()[1];
        }
        if(aNode.getCoords().length>2){
            if(aNode.getCoords()[2]>maximum_Z_coordinate)maximum_Z_coordinate=aNode.getCoords()[2];
            if(aNode.getCoords()[2]<minimum_Z_coordinate)minimum_Z_coordinate=aNode.getCoords()[2];
        }
        theNodes.put(aNode.getID(), aNode);
    }
    
    public void putMaterial(Material aMaterial) {
        // check for id
        if(theMaterials.containsKey(aMaterial.getID())){
            System.out.println("warning:    Material with id:"+aMaterial.getID());
            System.out.println("            already defined in Domain with id:"+this.id);
        }else{
            theMaterials.put(aMaterial.getID(), aMaterial);
        }
        
    }
    
    public void putCrossSection(CrossSection aCrossSection) {
        theCrossSections.put(aCrossSection.getID(), aCrossSection);
    }
    
    public void putElement(Element aElement) {
        theElements.put(aElement.getID(), aElement);
    }
    
    public void putConstraintElement(ConstraintElement aConstraintElement) {
        theConstraintElements.put(aConstraintElement.getID(), aConstraintElement);
    }
    
    public void putLoadCase(LoadCase aLoadCase) {
        theLoadCases.put(aLoadCase.getID(), aLoadCase);
        this.setNodesForLoadCase(aLoadCase.getID(), aLoadCase.getNumOfIncrements());
    }
    
    public void putInitialCondition(InitialCondition aInitialCondition) {
        theInitialConditions.put(aInitialCondition.getID(), aInitialCondition);
    }
    
    public void setEigenValues(double[] evs) {
        int order=evs.length;
        this.eigenvalues = new double[order];
        for(int i=0; i<order; i++){
            this.eigenvalues[i]=evs[i];
        }
    }
    
    public int getID(){
        return this.id;
    }
    
    public Node getNode(int id_){
        return theNodes.get(id_);
    }
    
    public Map getNodes(){
        return this.theNodes;
    }

    public void setNodes(Map<Integer,Node> theNodes){
        this.theNodes=theNodes;
    }
    
    public Material getMaterial(int id_){
        return theMaterials.get(id_);
    }
    
    public CrossSection getCrossSection(int id_){
        return theCrossSections.get(id_);
    }
    
    public Map getElements(){
        return this.theElements;
    }

    public void setElements(Map<Integer,Element> theElements){
        this.theElements=theElements;
    }
    
    public Element getElement(int id_){
        return this.theElements.get(id_);
    }
    
    
    public Map getConstraintElements(){
        return this.theConstraintElements;
    }
    
    public Element getConstraintElement(int id_){
        return this.theConstraintElements.get(id_);
    }
    
    public LoadCase getLoadCase(int id_){
        return theLoadCases.get(id_);
    }
    
    public Map getLoadCases(){
        return this.theLoadCases;
    }
    
    public InitialCondition getInitialCondition(int id_){
        return theInitialConditions.get(id_);
    }
    
    public Map getInitialConditions(){
        return this.theInitialConditions;
    }
    
    public int getNdofs(){
        return this.theNodes.get(1).getNdofs();
    }
    
    public int getNeigs(){
        int n=0;
        if(this.eigenvalues!=null)n=eigenvalues.length;
        return n;
    }
    
    public int getNumElement(){
        return this.theElements.size();
    }
    
    public int getNumQuadElement(){
        int count=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            if(theElement.getClass()==PlaneStressQuad.class || theElement.getClass()==PlaneStrainQuad.class)count+=1;
        }
        return count;
    }
    
    public int getNumLineElement(){
        int count=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            if(theElement.getClass()==Truss2d.class || theElement.getClass()==EBeam2d.class
                     || theElement.getClass()==EBeam3d.class)count+=1;
        }
        return count;
    }

    public int getNumNodes(){
        return this.theNodes.size();
    }
    
    public void clearDomain(){
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.clear();
        }
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.clear();
        }
    }
    
    public void clearDeltaDisps(){
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            theNode.clearDeltaDisps();
        }
    }
    
    public void DomainInitialConditions(){
        double[] initDisps = new double[this.getNdofs()];
        for(int i=0; i<initDisps.length; i++){initDisps[i]=0.;}
        for(Iterator<InitialCondition> it=this.theInitialConditions.values().iterator(); it.hasNext();){
            InitialCondition theInitialCondition = it.next();
            int dof = theInitialCondition.getNode().getFtable()[theInitialCondition.getDOF()-1];
            initDisps[dof-1]=theInitialCondition.getInitialValue();
            switch(theInitialCondition.getType()){
                case 0: theInitialCondition.getNode().updateDisps(initDisps); break;
                case 1: theInitialCondition.getNode().updateVelcs(initDisps); break;
                default: System.out.println("non properly defined initial condition "+theInitialCondition.getID());
            }
            
        } 
    }
    
    public int getNumOfPlastifiedElems(){
        int theNum=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            if(theElement.ElementPlastified()){theNum++ ;}
        }
        return theNum;
    }
    
    public double getPlasticNRG(){
        double val=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            val+=theElement.PlasticNRG();
        }
        return val;
    }

    public double getuKu(){
        double val=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            val+=theElement.getuKu();
        }
        return val;
    }
    
    public double getvMv(){
        double val=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            val+=theElement.getvMv();
        }
        return val;
    }
    
    public double getuCv(double am, double ak){
        double val=0;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            val+=theElement.getuCv(am, ak);
        }
        return val;
    }

    public double getuKu_constraints(){
        double val=0;
        for(Iterator<ConstraintElement> it=this.theConstraintElements.values().iterator(); it.hasNext();){
            val+=it.next().getuKu();
        }
        return val;
    }

    public void printNodesGeom(){
        System.out.println("Nodes of domain with id= "+this.id);
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            aNode.print();
        }
    }

    public void printElementsConnectivity(){
        System.out.println("Elements of Domain with id:"+this.id);
        for(Iterator<Element>it=this.theElements.values().iterator();it.hasNext();){
            it.next().print();
        }
    }

    public double getmaxX(){
        double val=0.;
        double max,min;
        max=0.; min=0.;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(aNode.getCoords()[0]>max)max=aNode.getCoords()[0];
            if(aNode.getCoords()[0]<min)min=aNode.getCoords()[0];
        }
        val=max-min;
        return val;
    }
    
    public double[] getminmaxX(){
        double[] val=new double[2];
        double max,min;
        max=0.; min=0.;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(aNode.getCoords()[0]>max)max=aNode.getCoords()[0];
            if(aNode.getCoords()[0]<min)min=aNode.getCoords()[0];
        }
        val[0]=min; val[1]=max;
        return val;
    }

    public double getmaxY(){
        double val=0.;
        double max,min;
        max=0.; min=0.;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(aNode.getCoords()[1]>max)max=aNode.getCoords()[1];
            if(aNode.getCoords()[1]<min)min=aNode.getCoords()[1];
        }
        val=max-min;
        return val;
    }
    
    public double[] getminmaxY(){
        double[] val=new double[2];
        double max,min;
        max=0.; min=0.;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            if(aNode.getCoords()[1]>max)max=aNode.getCoords()[1];
            if(aNode.getCoords()[1]<min)min=aNode.getCoords()[1];
        }
        val[0]=min; val[1]=max;
        return val;
    }

    public void printConstraintElements(){
        for(Iterator<ConstraintElement>it=this.theConstraintElements.values().iterator();it.hasNext();){
            ConstraintElement aConstraintElement=it.next();
            aConstraintElement.printC();
        }
    }
    
    public double[] getDomainBVNorm(){
        double[] norm = new double[2];
        norm[0]=0.; norm[1]=0.;
        double[] elemnorm;
        double maxU=-10.e100;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node theNode=it.next();
//            double curU=Math.sqrt(theNode.getDisp(1)*theNode.getDisp(1)+theNode.getDisp(2)*theNode.getDisp(2));
            double curU=Math.sqrt(theNode.getDispsTrial()[0]*theNode.getDispsTrial()[0]+theNode.getDispsTrial()[1]*theNode.getDispsTrial()[1]);
            if(curU>=maxU)maxU=curU;
        }
        for(Iterator<Element>it=this.theElements.values().iterator();it.hasNext();){
            Element theElem=it.next();
            elemnorm = theElem.getBVNorm(maxU);
            norm[0]+=elemnorm[0];
            norm[1]+=elemnorm[1];
        }
        return norm;
    }
    
    public int find(double x){
        double[] array = new double[2];
        array[0]=x; array[1]=0.0;
        return find(array);
    }
    
    public int find(double x, double y){
        double[] array = new double[2];
        array[0]=x; array[1]=y;
        return find(array);
    }
    
    public int find(double[] array){
        int nid=-1;
        double minL=10.e100;
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node theNode=it.next();
            double dist=0.;
            for(int i=0;i<array.length;i++){
                dist+=(theNode.getCoords()[i]-array[i])*(theNode.getCoords()[i]-array[i]);
            }
            dist=Math.sqrt(dist);
            if( dist<minL){minL=dist; nid=theNode.getID();}
        }
        return nid;
    }
    
    public double getDissipation(){
        double val=0.;
        for(Iterator<Element>it=this.theElements.values().iterator();it.hasNext();){
            Element theElem=it.next();
            val+=theElem.DamageDissipation();
        }
        return val;
    }
    
    public double[] getDisplacements(){
        double[] disps =new double[this.getNdofs()];
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node theNode=it.next();
            for(int i=0;i<theNode.getFtable().length;i++){
                if(theNode.getFtable()[i]>0 && theNode.getFtable()[i]<=this.getNdofs()){
//                    System.err.println(theNode.getFtable()[i]);
                    disps[theNode.getFtable()[i]-1]=theNode.getDisp(i+1);
                }
            }
        }
        return disps;
    }
    
    public double[] getVelocities(){
        double[] disps =new double[this.getNdofs()];
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node theNode=it.next();
            for(int i=0;i<theNode.getFtable().length;i++){
                if(theNode.getFtable()[i]>0 && theNode.getFtable()[i]<=this.getNdofs()){
//                    System.err.println(theNode.getFtable()[i]);
                    disps[theNode.getFtable()[i]-1]=theNode.getVelc(i+1);
                }
            }
        }
        return disps;
    }
    
    public double maximum_X_coordinate(){return this.maximum_X_coordinate;}
    
    public double maximum_Y_coordinate(){return this.maximum_Y_coordinate;}
    
    public double maximum_Z_coordinate(){return this.maximum_Z_coordinate;}
    
    public double minimum_X_coordinate(){return this.minimum_X_coordinate;}
    
    public double minimum_Y_coordinate(){return this.minimum_Y_coordinate;}
    
    public double minimum_Z_coordinate(){return this.minimum_Z_coordinate;}
    
    public void setNodesForLoadCase(int LC, int numSteps){
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            it.next().putLoadCase(LC, numSteps);
        }
    }

    public int getNumGridX() {
        return this.NumGridX;
    }

    public int getNumGridY() {
        return this.NumGridY;
    }
    
    public void setNumGridXY(int n){
        setNumGridXY(n,n);
    }
    
    public void setNumGridXY(int nx, int ny){
        NumGridX=nx;
        NumGridY=ny;
    }

    void setEigenLoadCases(int order) {
        LoadCase aLC;
        for(int i=1;i<=order;i++){aLC = new LoadCase(-i); this.putLoadCase(aLC);}
    }
    
    public int getNumOfLoadCases(){return this.theLoadCases.size();}
    
    public double[] getEigenValues(){
        return this.eigenvalues;
    }
    
    public void printEigenValues(){
        System.out.println("Eigenvalues of domain with id: "+this.id);
        System.out.println("number of eigval"+"     "+"eigenvalues     "+"     "+"eigenfrequencies"+"     "+"eigenperiods    ");
        System.out.println("----------------"+"     "+"----------------"+"     "+"----------------"+"     "+"----------------");
        for(int i=0; i<this.eigenvalues.length; i++){
             //System.out.println(Math.sqrt(eig_vals[i]));
             DecimalFormat df = new DecimalFormat("####");
             //System.out.print(df.format(i+1)+"     ");
             System.out.print((i+1)+"     ");
             df = new DecimalFormat("##.00#");
//             System.out.print(df.format(eigenvalues[i])+"     ");
//             System.out.print(df.format(Math.sqrt(eigenvalues[i]))+"     ");
//             System.out.println(df.format(2.0*Math.PI/Math.sqrt(eigenvalues[i]))+"     ");
             System.out.print((eigenvalues[i])+"     ");
             System.out.print((Math.sqrt(eigenvalues[i]))+"     ");
             System.out.println((2.0*Math.PI/Math.sqrt(eigenvalues[i]))+"     ");
         }
    }
    
    public double getMaxKcomponent(){
        double val=Double.NEGATIVE_INFINITY;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            for(int i=0;i<theElement.getElementNdofs();i++){
                val=Math.max(val, theElement.getK().get(i, i));
            }
        }
        return val;
    }
    
    public double[] getMaxKMcomponent(){
        double[] val=new double[2];
        val[0]=Double.NEGATIVE_INFINITY;
        val[1]=Double.NEGATIVE_INFINITY;
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            for(int i=0;i<theElement.getElementNdofs();i++){
                val[0]=Math.max(val[0], theElement.getK().get(i, i));
                val[1]=Math.max(val[1], theElement.getM().get(i, i));
            }
        }
        return val;
    }
    
    public void activateConsMass(){
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.activateConsMass();
        }
    }
    
    public void deactivateConsMass(){
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            theElement.deactivateConsMass();
        }
    }
    
    public double[][] getInitDispVec(){
        double[][] vals = new double[this.getNumDOFS()][1];
        for(Iterator<InitialCondition> it=this.theInitialConditions.values().iterator(); it.hasNext();){
            InitialCondition theIC = it.next();
            if(theIC.type==0){
                double val=theIC.getInitialValue();
                int dof=this.getNode(theIC.getNode().getID()).getFtable()[theIC.getDOF()-1];
                if(dof>0)vals[dof-1][0]=val;
            }
        }
        return vals;
    }
    
    public double[][] getInitVelcVec(){
        double[][] vals = new double[this.getNumDOFS()][1];
        for(Iterator<InitialCondition> it=this.theInitialConditions.values().iterator(); it.hasNext();){
            InitialCondition theIC = it.next();
            if(theIC.type==1){
                double val=theIC.getInitialValue();
                int dof=this.getNode(theIC.getNode().getID()).getFtable()[theIC.getDOF()-1];
                if(dof>0)vals[dof-1][0]=val;
            }
        }
        return vals;
    }
    
    public double[][] getEigenVec(int eig){
        double[][] vals = new double[this.getNumDOFS()][1];
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node theNode = it.next();
            for(int i=0;i<7;i++){
                int dof=theNode.getFtable()[i];
                double val=theNode.getLoadCaseDisps(-eig)[i];
                if(dof>0)vals[dof-1][0]=val;
            }
        }
        return vals;
    }
    
    public int getNumDOFS(){
        int num=0;
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            num+= it.next().getNdofs_ofNode();
        }
        return num;
    }
    
    public AbstractMatrix getK(){
        int N=getNumDOFS();
        AbstractMatrix theMatrix = new AbstractMatrix(N,N);
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getK();
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        
        for(Iterator<ConstraintElement> it=this.theConstraintElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getK();
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        return theMatrix;
    }
    
    public AbstractMatrix getDK(int param){
        return getDK(param, 0);
    }
    
    public AbstractMatrix getDK(int param, int pid){
        int N=getNumDOFS();
        AbstractMatrix theMatrix = new AbstractMatrix(N,N);
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getDK(param,pid);
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        
        for(Iterator<ConstraintElement> it=this.theConstraintElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getDK(param,pid);
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        return theMatrix;
    }
    
    public AbstractMatrix getM(){
        int N=getNumDOFS();
        AbstractMatrix theMatrix = new AbstractMatrix(N,N);
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM();
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        
        for(Iterator<ConstraintElement> it=this.theConstraintElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getM();
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        return theMatrix;
    }
    
    public AbstractMatrix getDM(int param){
        return getDM(param, 0);
    }
    
    public AbstractMatrix getDM(int param, int pid){
        int N=getNumDOFS();
        AbstractMatrix theMatrix = new AbstractMatrix(N,N);
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getDM(param,pid);
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        
        for(Iterator<ConstraintElement> it=this.theConstraintElements.values().iterator(); it.hasNext();){
            Element elem = it.next();
            AbstractMatrix mat ;
            mat=elem.getDM(param,pid);
            int[] theFtable=elem.getFtable();
            
            int s=theFtable.length;
            for(int i=0; i<s; i++){
                for(int j=0; j<s; j++){
                    theMatrix.addVal(theFtable[i]-1, theFtable[j]-1, mat.get(i, j));
                }
            }
        }
        return theMatrix;
    }
}
