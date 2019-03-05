/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class Interface {
    protected int id;
    protected Map<Integer,Node> theNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,Element> theElements = new TreeMap<Integer,Element>();
    protected Material theMaterial;
    protected Spring NormalSpring;
    protected Spring TangentSpring;
    protected EnergyDissipator theDissipator;
    protected Map<Integer,Node> theConnectedNodes = new TreeMap<Integer,Node>();
    protected Map<Integer,InterfaceNode> theInterfaceNodes = new TreeMap<Integer,InterfaceNode>();
    protected Map<Integer,InterfaceElement> theInterfaceElements = new TreeMap<Integer,InterfaceElement>();
    protected Map<Integer,MinimizationConstraintEquation> theMinimizationEquations = new TreeMap<Integer,MinimizationConstraintEquation>();
    protected Domain Domain1,Domain2=null;
    protected AbstractMatrix NMatrix,bMatrix,TMatrix;
    private boolean setMats=true;
    
    protected boolean xProjection=true;
    protected boolean closedCurve=false;

    protected boolean zconstant = false;
    protected boolean sconstant = false;
    protected boolean pconstant = false;

    protected boolean zcontiuous = true;
    protected boolean scontiuous = true;
    protected boolean pcontiuous = true;

    protected boolean tangential = false;
    protected boolean normal = false;

    protected boolean damage = false;
    protected boolean slip = false;
    protected boolean plast = false;
    
    protected boolean friction=false;
    
    protected boolean SplitNormalDisp = false;
    protected boolean SplitTangentialDisp = false;

    protected int sumDOFS;
    protected int utDOFS;
    protected int unDOFS;
    protected int zDOFS;
    
    protected int sDOFS;
    protected int s_posDOFS;
    protected int s_negDOFS;
    
    protected int pDOFS;
    protected int p_posDOFS;
    protected int p_negDOFS;
    
    protected int un_posDOFS;
    protected int un_negDOFS;
    
    protected int ut_posDOFS;
    protected int ut_negDOFS;
    
    protected int z_posDOFS;
    protected int z_negDOFS;
    
    private int beginDOF=0;
    private int beginMinEqId=0;
    
    private Double SlipRefernce=null;
    private Double PlastRefernce=null;
    private Double NormalRefernce=null;
    private Double TangentRefernce=null;
    
    private double minlength;
    private double eps=1.e-2;
    
    private boolean TangentPairs=false;
    
    private int INodeIDstart, INodeIDend;
    
    private boolean rigidnormal=false;
    private boolean zeroEndTangentTraction=false;
    private boolean zeroStartTangentTraction=false;
    
    private boolean ReleaseNormal=false;
    
    private boolean permitNCmesh=false;
    
    private boolean PlasticIncrement=true;
    
    private boolean OppositeNormal=false;
    
    private double tau;

    public Interface(){}

    public Interface(int id){
        this.id=id;
    }

    public Interface(int id, Domain A, Domain B){
        this.id=id;
        this.Domain1=A;
        this.Domain2=B;
        MinimizationConstraintEquation.setMulti();
    }

    public Interface(int id, Domain A){
        this.id=id;
        this.Domain1=A;
    }

    public void setDomains(Domain A, Domain B){
        this.Domain1=A;
        this.Domain2=B;
        MinimizationConstraintEquation.setMulti();
    }

    public void setDomains(Domain A){
        this.Domain1=A;
    }
    
    public void setTimeStep(double dt){this.tau=dt;}

    public int getID(){
        return this.id;
    }

    public void addMaterial(Material aMaterial){
        this.theMaterial=aMaterial;
    }


    public void putNode(Node aNode){
        this.theNodes.put(aNode.getID(), aNode);
    }

    public void putElement(Element anElement){
        this.theElements.put(anElement.getID(), anElement);
        //for(int i=1;i<=anElement.getNumNodes();i++){
        //    if(!this.theNodes.containsKey(anElement.getNodeHier(i).getID()))putNode(anElement.getNodeHier(i));
        //}
    }

    public Element getElement(int id){
        return this.theElements.get(id);
    }

    public Node getNode(int id){
        return this.theNodes.get(id);
    }

    public Map<Integer,Element> getElements(){
        return this.theElements;
    }

    public Map<Integer,InterfaceElement> getInterfaceElements(){
        return this.theInterfaceElements;
    }

    public Map<Integer,Node> getNodes(){
        return this.theNodes;
    }

    public Map<Integer,InterfaceNode> getInterfaceNodes(){
        return this.theInterfaceNodes;
    }

    public void setSpaceIntegrator(SpaceIntegrator SI){
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement theElement = it.next();
            theElement.setSpaceIntegrator(SI);
        }
    }

    public void printElementsConectivity(){
        int nnode;
        System.out.println("Interface id="+this.id+", has "+this.theElements.size()+" elements with id and nodes:");
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            Element theElement = it.next();
            System.out.print(theElement.getID()+" ");
            nnode=theElement.getNumNodes();
            for(int i=1;i<=nnode;i++){System.out.print(theElement.getNodeHier(i).getID()+" ");}
            System.out.println();
        }
    }
    
    public void printElementsConectivity(javax.swing.JTextArea TextArea){
        int nnode;
        TextArea.append("Domain id="+this.id+", has "+this.theElements.size()+" elements with id and nodes:"+'\n');
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            ELine theElement = (ELine) it.next();
            //TextArea.append(theElement.getID()+" "+" of Length="+theElement.getLength()+'\n');
            nnode=theElement.getNumNodes();
            for(int i=1;i<=nnode;i++){
                TextArea.append(theElement.getNodeHier(i).getID()+" ");
            }
            TextArea.append(""+'\n');
        }
    }
    
    public void printNodesGeom(){
        System.out.println("Nodes of interface with id= "+this.id);
        for(Iterator<Node>it=this.theNodes.values().iterator();it.hasNext();){
            Node aNode=it.next();
            aNode.print();
        }
    }

    public void printINodesGeom(){
        System.out.println("Interface Nodes of interface with id= "+this.id);
        System.out.println("id"+" "+"mainNode"+" "+"inequation"+" "+"x"+" "+"y"+" "+"z"+" "+"base"+" "+"twin"+" "+"onElem"+" "+"gap"+" "+"nx"+" "+"ny"+" "+"nz");
        for(Iterator<InterfaceNode>it=this.theInterfaceNodes.values().iterator();it.hasNext();){
            InterfaceNode aNode=it.next();
            aNode.print();
        }
    }
    
    public void printINodesResponse(){
        System.out.println("Interface Nodes of interface with id= "+this.id);
        System.out.println("--------------------------------------------------------------");
        for(Iterator<InterfaceNode>it=this.theInterfaceNodes.values().iterator();it.hasNext();){
            InterfaceNode aNode=it.next();
            aNode.printResponse();
        }
    }

    public void printINodeszDOFS(){
        System.out.println("zDOFS of Nodes of interface with id= "+this.id);
        for(Iterator<InterfaceElement>it=this.theInterfaceElements.values().iterator();it.hasNext();){
            InterfaceElement aElement=it.next();
            aElement.getINodeHier(1).printzDOFS();
            aElement.getINodeHier(2).printzDOFS();
        }
    }

    public void setConnectedNodes(){
        InterfaceNode INode;
        Node Node1;
        double length;
        minlength=10.e20;
        Element theElement1;
        
        for(Iterator<Element> it=this.theElements.values().iterator(); it.hasNext();){
            theElement1 = it.next();
            length = ((ELine) theElement1).getLength();
            if(length<minlength)minlength=length;
        }
        
        double xsi, dist;
        int[] ElemsOfNode;
        double[] adoubl;
        int tid = 0;

        
        boolean onElem,sameDomain;
        double coef;
        
        double u0=1.;
        double u1=1.;
        double u2=-1.;
            
        for(Iterator<Node> it=this.theNodes.values().iterator(); it.hasNext();){
            Node1 = it.next();
            if(this.Domain2==null){
                if(Node1.isMain()==false){
                    System.err.println("Interface with single domain of non-Main nodes");
                    System.err.println("Analysis terminated");
                    System.err.println("Set domain into Main and restart analysi");
                    System.exit(1500001);
                }
            }
            ElemsOfNode = Node1.getConnectedElementsIds();
            if(this.theElements.containsKey(ElemsOfNode[0]) || this.theElements.containsKey(ElemsOfNode[1])){
                INode = new InterfaceNode(Node1.getCoordinates(),Node1);

                if(this.theElements.containsKey(ElemsOfNode[0]) && (!this.theElements.containsKey(ElemsOfNode[1])) ){
                    adoubl = new double[3];
                    for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(Node1.getID())[i];}
                    INode.setNormal(adoubl, u0 );
                }else if((!this.theElements.containsKey(ElemsOfNode[0]))&& this.theElements.containsKey(ElemsOfNode[1])){
                    adoubl = new double[3];
                    for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(Node1.getID())[i];}
                    INode.setNormal(adoubl, u0 );
                }else{
                    //System.err.println("!!!!!!!!!!!!!!!!! "+INode.getBaseID()+" suspicious normal generation (1) !!!!!!!!!!!!!!!!!");
                    adoubl = new double[3];
                    for(int i=0;i<3;i++){
                        adoubl[i]=0.;
                        adoubl[i]=0.5*((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(Node1.getID())[i]+
                           0.5*((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(Node1.getID())[i];
                    }
                    INode.setNormal(adoubl, u0 );
                }

                if(Node1.isMain()){INode.setMain();}
                tid = INode.getID();
                
                this.theInterfaceNodes.put(INode.getID(), INode);
                for(Iterator<Element> et=this.theElements.values().iterator(); et.hasNext();){
                    theElement1 = et.next();
                    onElem=false; sameDomain=false;
                    for(int e=0;e<ElemsOfNode.length;e++){if(ElemsOfNode[e]==theElement1.getID())onElem=true;}
                    sameDomain=this.InSameDomain(INode, theElement1);
                    if(!onElem && !sameDomain){
                        xsi=theElement1.getCoordMinDistOfNode(Node1)[0];
                        if(xsi>(-1.)+eps && xsi<(1.)-eps){
                            dist = theElement1.getMinDistOfNode(Node1);
                            if(this.getINodeWithTwinID(tid)!=null){
                                if(INode.getDist(this.getINodeWithTwinID(tid))>=dist){
                                    INode = new InterfaceNode(theElement1.getGlobalCoordinates(theElement1.getCoordMinDistOfNode(Node1)),theElement1);
                                    if(this.Domain2!=null)INode.setTwinINode(theInterfaceNodes.get(tid));
                                    adoubl = ((ELine) theElement1).getNormal(xsi); if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                                    INode.setNormal(adoubl,coef);
                                    if(theElement1.isMain()){INode.setMain();}
                                    this.theInterfaceNodes.put(INode.getID(), INode);
                                    if(this.Domain2!=null)theInterfaceNodes.get(tid).setTwinINode(INode);
                                    theInterfaceNodes.get(tid).setNormal(adoubl,coef);
                                }
                            }else{
                                INode = new InterfaceNode(theElement1.getGlobalCoordinates(theElement1.getCoordMinDistOfNode(Node1)),theElement1);
                                if(this.Domain2!=null)INode.setTwinINode(theInterfaceNodes.get(tid));
                                adoubl = ((ELine) theElement1).getNormal(xsi);if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                                INode.setNormal(adoubl,coef);
                                if(theElement1.isMain()){INode.setMain();}
                                this.theInterfaceNodes.put(INode.getID(), INode);
                                if(this.Domain2!=null)theInterfaceNodes.get(tid).setTwinINode(INode);
                                theInterfaceNodes.get(tid).setNormal(adoubl,coef);
                            }

                        }else if(Math.abs(xsi-1.)<=eps){
                            if(this.getINodeWithBaseID(theElement1.getNodeHier(2).getID())!=null){
                                if(INode.getTwinID()!=0){
                                    if(INode.getDist(theElement1.getNodeHier(2))<=INode.getDist(INode.getTwin())){
                                        ElemsOfNode = theElement1.getNodeHier(2).getConnectedElementsIds();
                                        if(this.theElements.containsKey(ElemsOfNode[0]) && !this.theElements.containsKey(ElemsOfNode[1])){
                                            adoubl = new double[3];
                                            for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(2).getID())[i];}
                                            INode.setNormal(adoubl, u0 );
                                        }else if(!this.theElements.containsKey(ElemsOfNode[0]) && this.theElements.containsKey(ElemsOfNode[1])){
                                            adoubl = new double[3];
                                            for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(2).getID())[i];}
                                            INode.setNormal(adoubl, u0 );
                                        }else{
                                            //System.err.println("!!!!!!!!!!!!!!!!! "+INode.getBaseID()+" suspicious normal generation (2) !!!!!!!!!!!!!!!!!");
                                            adoubl = new double[3];
                                            for(int i=0;i<3;i++){
                                                adoubl[i]=0.;
                                                adoubl[i]=0.5*((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(2).getID())[i]+
                                                   0.5*((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(2).getID())[i];
                                            }
                                        }
                                        if(this.getINodeWithBaseID(theElement1.getNodeHier(2).getID()).getID()!=INode.getID()){
                                        if(this.Domain2!=null)this.getINodeWithBaseID(theElement1.getNodeHier(2).getID()).setTwinINode(INode);
                                        if(this.Domain2!=null)INode.setTwinINode(this.getINodeWithBaseID(theElement1.getNodeHier(2).getID()));}
                                        if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                                        INode.setNormal(adoubl,coef);
                                        INode.getTwin().setNormal(adoubl,coef);
                                    }
                                }else{
                                    ElemsOfNode = theElement1.getNodeHier(2).getConnectedElementsIds();
                                    if(this.theElements.containsKey(ElemsOfNode[0]) && (!this.theElements.containsKey(ElemsOfNode[1]))){
                                        adoubl = new double[3];
                                        for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(2).getID())[i];}
                                        INode.setNormal(adoubl, u0 );
                                    }else if((!this.theElements.containsKey(ElemsOfNode[0]))&& this.theElements.containsKey(ElemsOfNode[1])){
                                        adoubl = new double[3];
                                        for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(2).getID())[i];}
                                        INode.setNormal(adoubl, u0 );
                                    }else{
                                        //System.err.println("!!!!!!!!!!!!!!!!! "+INode.getBaseID()+" suspicious normal generation (3) !!!!!!!!!!!!!!!!!");
                                        adoubl = new double[3];
                                        for(int i=0;i<3;i++){
                                            adoubl[i]=0.;
                                            adoubl[i]=0.5*((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(2).getID())[i]+
                                               0.5*((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(2).getID())[i];
                                        }
                                    }
                                    if(this.getINodeWithBaseID(theElement1.getNodeHier(2).getID()).getID()!=INode.getID()){
                                    if(this.Domain2!=null)this.getINodeWithBaseID(theElement1.getNodeHier(2).getID()).setTwinINode(INode);
                                    if(this.Domain2!=null)INode.setTwinINode(this.getINodeWithBaseID(theElement1.getNodeHier(2).getID()));}
                                    if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                                    INode.setNormal(adoubl,coef); 
                                    if(this.Domain2!=null)INode.getTwin().setNormal(adoubl,coef);
                                }
                                
                            
                            }
                        }else if(Math.abs(xsi+1.)<=eps){
                            if(this.getINodeWithBaseID(theElement1.getNodeHier(1).getID())!=null){
                                if(INode.getTwinID()!=0){
                                    if(INode.getDist(theElement1.getNodeHier(1))<=INode.getDist(INode.getTwin())){
                                        ElemsOfNode = theElement1.getNodeHier(1).getConnectedElementsIds();
                                        if(this.theElements.containsKey(ElemsOfNode[0]) && !this.theElements.containsKey(ElemsOfNode[1])){
                                            adoubl = new double[3];
                                            for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(1).getID())[i];}
                                            INode.setNormal(adoubl, u0 );
                                        }else if(!this.theElements.containsKey(ElemsOfNode[0]) && this.theElements.containsKey(ElemsOfNode[1])){
                                            adoubl = new double[3];
                                            for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(1).getID())[i];}
                                            INode.setNormal(adoubl, u0 );
                                        }else{
                                            //System.err.println("!!!!!!!!!!!!!!!!! "+INode.getBaseID()+" suspicious normal generation (4) !!!!!!!!!!!!!!!!!");
                                            adoubl = new double[3];
                                            for(int i=0;i<3;i++){
                                                adoubl[i]=0.;
                                                adoubl[i]=0.5*((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(1).getID())[i]+
                                                   0.5*((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(1).getID())[i];
                                            }
                                        }
                                        if(this.getINodeWithBaseID(theElement1.getNodeHier(1).getID()).getID()!=INode.getID()){
                                        if(this.Domain2!=null)this.getINodeWithBaseID(theElement1.getNodeHier(1).getID()).setTwinINode(INode);
                                        if(this.Domain2!=null)INode.setTwinINode(this.getINodeWithBaseID(theElement1.getNodeHier(1).getID()));}
                                        if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                                        INode.setNormal(adoubl,coef);
                                        if(this.Domain2!=null)INode.getTwin().setNormal(adoubl,coef);
                                    }
                                }else{
                                    ElemsOfNode = theElement1.getNodeHier(1).getConnectedElementsIds();
                                    if(this.theElements.containsKey(ElemsOfNode[0]) && !this.theElements.containsKey(ElemsOfNode[1])){
                                        adoubl = new double[3];
                                        for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(1).getID())[i];}
                                        INode.setNormal(adoubl, u0 );
                                    }else if(!this.theElements.containsKey(ElemsOfNode[0]) && this.theElements.containsKey(ElemsOfNode[1])){
                                        adoubl = new double[3];
                                        for(int i=0;i<3;i++){adoubl[i]=0.; adoubl[i]=((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(1).getID())[i];}
                                        INode.setNormal(adoubl, u0 );
                                    }else{
                                        //System.err.println("!!!!!!!!!!!!!!!!! "+INode.getBaseID()+" suspicious normal generation (5) !!!!!!!!!!!!!!!!!");
                                        adoubl = new double[3];
                                        for(int i=0;i<3;i++){
                                            adoubl[i]=0.;
                                            adoubl[i]=0.5*((ELine) this.theElements.get(ElemsOfNode[0])).getNormal(theElement1.getNodeHier(1).getID())[i]+
                                               0.5*((ELine) this.theElements.get(ElemsOfNode[1])).getNormal(theElement1.getNodeHier(1).getID())[i];
                                        }
                                    }
                                    if(this.getINodeWithBaseID(theElement1.getNodeHier(1).getID()).getID()!=INode.getID()){
                                        if(this.Domain2!=null)this.getINodeWithBaseID(theElement1.getNodeHier(1).getID()).setTwinINode(INode);
                                        if(this.Domain2!=null)INode.setTwinINode(this.getINodeWithBaseID(theElement1.getNodeHier(1).getID()));
                                    }
                                    if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                                    INode.setNormal(adoubl,coef);
                                    if(this.Domain2!=null)INode.getTwin().setNormal(adoubl,coef);
                                }
                                
                            
                            }
                        }
                    }
                }
            }
            
        }
        
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            INode = it.next();
            if(!INode.isMain()){
                if(INode.getTwinID()!=0)this.theInterfaceNodes.get(INode.getTwinID()).setTwinINode(INode);
            }
        }

        // nodes exists, now construct interface elements
        InterfaceElement theIElem;
        dist=10e14;
        xsi=-dist;
        boolean themain;

        int nm=0, nnm=0;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            INode =  it.next();
            // also on some variables needed for uni-nodal constraints and CHECKING reasons
            if(INode.getBaseID()!=0){
                if(INode.isMain()){nm+=1;}else{nnm+=1;}
            }
            //////////////////////////////////////////////////////////////////////
        }

        if(nm>nnm){
            //System.out.println("extra equations on main solid");
            themain = true;
        }else{
            //System.out.println("extra equations on non main solid");
            themain = false;
        }


        int proj=1; if(this.xProjection)proj=0;
        tid=0;
        double xc=0., yc=0.;
        int refID=0;

        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            INode = it.next();
            if(this.getNumOfConnectedDomains()==2){
                if(INode.isMain()==themain && INode.getTwinID()!=0){
                    if(refID==0)refID=INode.getID();
                    tid+=1;
                    xc+=INode.getCoordinates()[0];
                    yc+=INode.getCoordinates()[1];
                }
            }
            else if(this.getNumOfConnectedDomains()==1){
                if(refID==0)refID=INode.getID();
                tid+=1;
                xc+=INode.getCoordinates()[0];
                yc+=INode.getCoordinates()[1];
            }
        }
        xc=xc/tid;
        yc=yc/tid;
        if(this.closedCurve)tid+=1;
        ElemsOfNode = new int[tid];

        if(this.closedCurve){
            ElemsOfNode[0]=refID;
            double angle;
            for(int i=1;i<ElemsOfNode.length;i++){
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    INode = it.next();
                    if(INode.getID()!=refID){
                        angle = Math.atan2(-INode.getCoordinates()[1]+yc, -INode.getCoordinates()[0]+xc)
                                -Math.atan2(-theInterfaceNodes.get(refID).getCoordinates()[1]+yc, -theInterfaceNodes.get(refID).getCoordinates()[0]+xc);
                        //angle = Math.acos( (theInterfaceNodes.get(refID).getCoordinates()[0]-xc)*(INode.getCoordinates()[0]-xc)+
                        //        (theInterfaceNodes.get(refID).getCoordinates()[1]-yc)*(INode.getCoordinates()[1]-yc));
                        if(this.getNumOfConnectedDomains()==2){
                            if(angle>xsi && angle<dist && 
                                    INode.getTwinID()!=0 && INode.isMain()==themain){
                                tid=INode.getID(); dist=angle;
                            }
                        }
                        else if(this.getNumOfConnectedDomains()==1){
                            if(angle>xsi && angle<dist && 
                                    INode.isMain()==themain){
                                tid=INode.getID(); dist=angle;
                            }
                        }
                    }
                }
                ElemsOfNode[i]=tid;
                xsi=dist;
                dist=10.e14;
            }
        }else{
            for(int i=0;i<ElemsOfNode.length;i++){
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    INode = it.next();
                    if(this.getNumOfConnectedDomains()==2){
                        if(INode.getCoordinates()[proj]>xsi && INode.getCoordinates()[proj]<dist && 
                                INode.getTwinID()!=0 && INode.isMain()==themain){
                            tid=INode.getID(); dist=INode.getCoordinates()[proj];
                        }
                    }
                    else if(this.getNumOfConnectedDomains()==1){
                        if(INode.getCoordinates()[proj]>xsi && INode.getCoordinates()[proj]<dist && 
                                INode.isMain()==themain){
                            tid=INode.getID(); dist=INode.getCoordinates()[proj];
                        }
                    }
                }
                ElemsOfNode[i]=tid;
                xsi=dist;
                dist=10.e14;
            }
        }
        if(this.closedCurve)ElemsOfNode[ElemsOfNode.length-1]=ElemsOfNode[0];
        
        InterfaceNode theINode1, theINode2;
        for(int i=0;i<ElemsOfNode.length-1;i++){
            theINode1=this.theInterfaceNodes.get(ElemsOfNode[i]);
            theINode2=this.theInterfaceNodes.get(ElemsOfNode[i+1]);
            if(i==0)this.INodeIDstart=theINode1.getID();
            if(i==(ElemsOfNode.length-2))this.INodeIDend=theINode2.getID();
            if(theINode1.getBaseID()!=0)theINode1.setInEquation();
            if(theINode2.getBaseID()!=0)theINode2.setInEquation();

            if(theINode1.getBaseID()==0){
                if(getINodeWithBaseID(theINode1.getOnElement().getNodeHier(1).getID())!=null)this.getINodeWithBaseID(theINode1.getOnElement().getNodeHier(1).getID()).setInEquation();
                if(getINodeWithBaseID(theINode1.getOnElement().getNodeHier(2).getID())!=null)this.getINodeWithBaseID(theINode1.getOnElement().getNodeHier(2).getID()).setInEquation();
            }
            if(theINode2.getBaseID()==0){
                if(getINodeWithBaseID(theINode2.getOnElement().getNodeHier(1).getID())!=null)this.getINodeWithBaseID(theINode2.getOnElement().getNodeHier(1).getID()).setInEquation();
                if(getINodeWithBaseID(theINode2.getOnElement().getNodeHier(2).getID())!=null)this.getINodeWithBaseID(theINode2.getOnElement().getNodeHier(2).getID()).setInEquation();
            }

            if(this.getNumOfConnectedDomains()==2){
                theINode1=this.theInterfaceNodes.get(this.theInterfaceNodes.get(ElemsOfNode[i]).getTwinID());
                theINode2=this.theInterfaceNodes.get(this.theInterfaceNodes.get(ElemsOfNode[i+1]).getTwinID());
                if(theINode1.getBaseID()!=0)theINode1.setInEquation();
                if(theINode2.getBaseID()!=0)theINode2.setInEquation();
            }

            if(theINode1.getBaseID()==0){
                if(getINodeWithBaseID(theINode1.getOnElement().getNodeHier(1).getID())!=null)this.getINodeWithBaseID(theINode1.getOnElement().getNodeHier(1).getID()).setInEquation();
                if(getINodeWithBaseID(theINode1.getOnElement().getNodeHier(2).getID())!=null)this.getINodeWithBaseID(theINode1.getOnElement().getNodeHier(2).getID()).setInEquation();
            }
            if(theINode2.getBaseID()==0){
                if(getINodeWithBaseID(theINode2.getOnElement().getNodeHier(1).getID())!=null)this.getINodeWithBaseID(theINode2.getOnElement().getNodeHier(1).getID()).setInEquation();
                if(getINodeWithBaseID(theINode2.getOnElement().getNodeHier(2).getID())!=null)this.getINodeWithBaseID(theINode2.getOnElement().getNodeHier(2).getID()).setInEquation();
            }

            theIElem = new IELine2(this.theInterfaceNodes.get(ElemsOfNode[i]),this.theInterfaceNodes.get(ElemsOfNode[i+1]));
            theInterfaceElements.put(theIElem.getID(), theIElem);
        }

        // set undefined normals if some node has orthogonal projection inside an interface element
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            INode = it.next();
            onElem = false;
            if(INode.AskIfInEquation() && INode.getNormal()==null){
                for(Iterator<InterfaceElement> et=this.theInterfaceElements.values().iterator(); et.hasNext();){
                    theIElem = et.next();
                    xsi=theIElem.getCoordMinDistOfINode(INode)[0];
                    if(xsi>=(-1.) && xsi<=(1.)){
                        onElem = true;
                        adoubl = ((IELine) theIElem).getAproximatedNormal(xsi); if(INode.isMain()){coef=u1;}else{coef=u2;}
                        INode.setNormal(adoubl, coef);
                    }
                }
                if(onElem==false){
                    Node1 = INode.getBase();
                    for(int i=0;i<Node1.getConnectedElementsIds().length;i++){
                        if(theElements.containsKey(Node1.getConnectedElementsIds()[i])){
                            theElement1 = theElements.get(Node1.getConnectedElementsIds()[i]);
                            tid=theElement1.getHierOfNode(Node1.getID());
                            switch(tid){
                                case 1: xsi=-1.; break;
                                case 2: xsi=1.; break;
                            }
                            adoubl = ((ELine) theElement1).getNormal(xsi);if(theElement1.isMain()){coef=u1;}else{coef=u2;}
                            INode.setNormal(adoubl, coef);
                        }
                    }
                }
            }
        }
        
//        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
//            InterfaceElement theElement = it.next();
//            for(int i=1;i<=2;i++){
//                theINode1 = theElement.getINodeHier(i);
//                theINode2 = theElement.getINodeHier(i).getTwin();
//                adoubl = new double[3]; adoubl[2]=0.;
//                adoubl[0]=theINode1.getCoordinates()[0]-theINode2.getCoordinates()[0];
//                adoubl[1]=theINode1.getCoordinates()[1]-theINode2.getCoordinates()[1];
//                if(theINode2.isMain()){adoubl[0]=-adoubl[0]; adoubl[1]=-adoubl[1];}
//                theINode1.setNormal(adoubl, u1);theINode2.setNormal(adoubl, u2);
//            }
//        }
        
        if(this.OppositeNormal)for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            INode = it.next();
            double[] OppositeNormalVec = new double[INode.getNormal().length];
            for(int i=0;i<INode.getNormal().length;i++){
                OppositeNormalVec[i]=-INode.getNormal()[i];
            }
            INode.setNormal(OppositeNormalVec, 1.);
        }
    }
    
    public void setOppositeNormal(){this.OppositeNormal=true;}

    public Map<Integer,Node> getConnectedNodes(){
        return this.theConnectedNodes;
    }

    public Map<Integer,MinimizationConstraintEquation> getMinimizationEquations(){
        return this.theMinimizationEquations;
    }

    
    public boolean ExistINodeWithBaseID(int theID){
        boolean respond=false;
        InterfaceNode theINode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            theINode = it.next();
            if(theINode.getBaseID()==theID)respond=true;
        }
        return respond;
    }

    public void printPath(){
        System.out.println("Intermediate Path of interface with id = "+this.id);
        //InterfaceNode theINode;
        int nnode;
        InterfaceElement theElement;
        /*for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            theINode = it.next();
            if(theINode.getTwinID()!=0 && theINode.isMain()){
                System.out.println(theINode.getID()+" "+theINode.getIntermediateCoords()[0]+" "+theINode.getIntermediateCoords()[1]);
            }
        }*/
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            theElement = it.next();
            //System.out.print(theElement.getID()+" ");
            nnode=theElement.getNumNodes();
            for(int i=1;i<=nnode;i++){
                //System.out.println(theElement.getINodeHier(i).getCoordinates()[0]+" "+theElement.getINodeHier(i).getCoordinates()[1]);
                System.out.println(theElement.getINodeHier(i).getIntermediateCoords()[0]+" "+theElement.getINodeHier(i).getIntermediateCoords()[1]);
            }
            //System.out.println();
        }
    }

    public InterfaceNode getINodeWithBaseID(int bid){
        InterfaceNode theRNode = null;
        InterfaceNode theINode ;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            theINode = it.next();
            if(theINode.getBaseID()==bid){
                theRNode=theINode;
            }
        }
        return theRNode;
    }

    public InterfaceNode getINodeWithTwinID(int bid){
        InterfaceNode theRNode = null;
        InterfaceNode theINode ;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            theINode = it.next();
            if(theINode.getTwinID()==bid){
                theRNode=theINode;
            }
        }
        return theRNode;
    }

    public void set_xProjection(boolean b){this.xProjection=b;}

    public void set_zConstant(boolean b){this.zconstant=b;}

    public void set_sConstant(boolean b){this.sconstant=b;}
    
    public void set_pConstant(boolean b){this.pconstant=b;}

    public void set_zContiuous(boolean b){this.zcontiuous=b;}

    public void set_sContiuous(boolean b){this.scontiuous=b;}
    
    public void set_pContiuous(boolean b){this.pcontiuous=b;}

    public void set_Damage(boolean b){this.damage=b;}

    public void set_Slip(boolean b){this.slip=b;}
    
    public void set_Plast(boolean b){this.plast=b;}

    public void set_Tangential(boolean b){this.tangential=b;}

    public void set_Normal(boolean b){this.normal=b;}
    
    public void set_Friction(boolean b){
        this.friction=b; 
//        this.SplitNormalDisp=b;
    }

    public void printIElementsConectivity(){
        int nnode;
        InterfaceElement theElement;
        double lenght;
        System.out.println("Interface id="+this.id+", has "+this.theInterfaceElements.size()+" elements with id and nodes:");
        System.out.println("Start INode id = "+this.INodeIDstart+" || End INode id = "+this.INodeIDend);
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            theElement = it.next();
            lenght=( (IELine) theElement).getLength();
            System.out.print(theElement.getID()+" "+lenght+" " );
            nnode=theElement.getNumNodes();
            for(int i=1;i<=nnode;i++){System.out.print(theElement.getINodeHier(i).getID()+" ");}
            if(lenght<this.eps*this.minlength)System.out.print(" !!!!!!!!!!!");
            System.out.println();
        }
    }

    public void setNormalSpring(Spring aSpring){
        this.NormalSpring=aSpring;
    }

    public void setTangentSpring(Spring aSpring){
        this.TangentSpring=aSpring;
    }

    public void setDissipator(EnergyDissipator aDis){
        this.theDissipator=aDis;
    }
    
    public EnergyDissipator getDissipator(){
        return this.theDissipator;
    }

    public Spring getNormalSpring(){
        return this.NormalSpring;
    }

    public Spring getTangentSpring(){
        return this.TangentSpring;
    }

    public EnergyDissipator getEnergyDissipator(){return this.theDissipator;}
    
    public void setBeginDOF(int beg){
        this.beginDOF=beg;
    }

    public int getBeginDOFs(){
        return this.beginDOF;
    }
    
    public void setBeginMinEqId(int beg){
        this.beginMinEqId=beg;
    }

    public int getBeginMinEqIds(){
        return this.beginMinEqId;
    }

    public void setDOFs(int steps){
        // initial check
        if(this.slip && this.tangential==false){
            System.err.println("ERROR ! Inteface of id = "+this.id+" exist slip without tangential dofs");
            System.exit(1);
        }
        if(this.plast && this.normal==false){
            System.err.println("ERROR ! Inteface of id = "+this.id+" exist plast without normal dofs");
            System.exit(1);
        }
        if(this.normal==false){
            System.err.println("WARNING ! Inteface of id = "+this.id+" without normal dofs");
        }
        if(this.zconstant==true && this.zcontiuous==true){
//            System.err.println("WARNING ! Inteface of id = "+this.id+" when constant z supposed, discontinuity between elements by default");
        }
        if(this.sconstant==true && this.scontiuous==true){
            System.err.println("WARNING ! Inteface of id = "+this.id+" when constant s supposed, discontinuity between elements by default");
        }
        if(this.pconstant==true && this.pcontiuous==true){
            System.err.println("WARNING ! Inteface of id = "+this.id+" when constant p supposed, discontinuity between elements by default");
        }
        ////////////////
        this.sumDOFS=0;
        this.utDOFS=0;
        this.unDOFS=0;
        this.zDOFS=0;
        this.sDOFS=0;
        this.s_posDOFS=0;
        this.s_negDOFS=0;
        this.pDOFS=0;
        this.p_posDOFS=0;
        this.p_negDOFS=0;
        this.un_negDOFS=0;
        this.un_posDOFS=0;
        this.ut_negDOFS=0;
        this.ut_posDOFS=0;
        
        int numElems;
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(this.normal){
                anInterfaceNode.init_un(steps);
                //if(this.SplitNormalDisp){anInterfaceNode.init_un_pos(steps); anInterfaceNode.init_un_neg(steps);}
                anInterfaceNode.init_un_pos(steps); anInterfaceNode.init_un_neg(steps);
            }
            
            if(this.tangential){
                anInterfaceNode.init_ut(steps);
                if(this.SplitTangentialDisp){
                    anInterfaceNode.init_ut_pos(steps); 
                    anInterfaceNode.init_ut_neg(steps);
                }
            }
            
            if(this.slip){
                numElems = anInterfaceNode.getNumOfConnectedIElements();
                if(numElems>0){
                    if(this.sconstant){anInterfaceNode.init_s(steps, numElems);}else{if(this.scontiuous){anInterfaceNode.init_s(steps, 1);}else{anInterfaceNode.init_s(steps, numElems);}}
                }
            }
            
            if(this.plast){
                numElems = anInterfaceNode.getNumOfConnectedIElements();
                if(numElems>0){
                    if(this.pconstant){anInterfaceNode.init_p(steps, numElems);}else{if(this.pcontiuous){anInterfaceNode.init_p(steps, 1);}else{anInterfaceNode.init_p(steps, numElems);}}
                }
            }
            
            if(this.damage){
                numElems = anInterfaceNode.getNumOfConnectedIElements();
                if(numElems>0){
                    if(this.zconstant){anInterfaceNode.init_z(steps, numElems);}else{if(this.zcontiuous){anInterfaceNode.init_z(steps, 1);}else{anInterfaceNode.init_z(steps, numElems);}}
                }
            }
        }


        if(this.normal){
            if(!this.SplitNormalDisp){
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode = it.next();
                    if(anInterfaceNode.AskIfInEquation()) {
                        anInterfaceNode.setunEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; unDOFS+=1;
                    }
                }
            }else{
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode = it.next();
                    if(anInterfaceNode.AskIfInEquation()) {
                        anInterfaceNode.setun_posEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; un_posDOFS+=1;
                    }
                }
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode = it.next();
                    if(anInterfaceNode.AskIfInEquation()) {
                        anInterfaceNode.setun_negEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; un_negDOFS+=1;
                    }
                }
            }
        }


        if(this.tangential){
            if(!this.SplitTangentialDisp){
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode = it.next();
                    if(anInterfaceNode.AskIfInEquation()) {
                        anInterfaceNode.setutEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; utDOFS+=1;
                    }
                }
            }else{
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode = it.next();
                    if(anInterfaceNode.AskIfInEquation()) {
                        anInterfaceNode.setut_posEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; ut_posDOFS+=1;
                    }
                }
                
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode = it.next();
                    if(anInterfaceNode.AskIfInEquation()) {
                        anInterfaceNode.setut_negEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; ut_negDOFS+=1;
                    }
                }
            }
        }

        IELine2 theIElem;
        // FROM HERE I TOOK DAMAGE

        if(this.slip){
            if(this.sconstant){
                for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                    theIElem = (IELine2) it.next();
                    theIElem.getINodeHier(1).sets_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    theIElem.getINodeHier(2).sets_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    sumDOFS+=1; s_posDOFS+=1;
                    
                    theIElem.getINodeHier(1).sets_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    theIElem.getINodeHier(2).sets_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    sumDOFS+=1; s_negDOFS+=1;
                }
            }else{
                if(this.scontiuous){
                    for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                        anInterfaceNode = it.next();
                        if(anInterfaceNode.getNumOfConnectedIElements()>0) {
                            anInterfaceNode.sets_posEFTable(sumDOFS+1+beginDOF); 
                            sumDOFS+=1; s_posDOFS+=1;
                            anInterfaceNode.sets_negEFTable(sumDOFS+1+beginDOF); 
                            sumDOFS+=1; s_negDOFS+=1;
                        }
                    }
                }else{
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        theIElem.getINodeHier(1).sets_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; s_posDOFS+=1;
                        theIElem.getINodeHier(2).sets_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; s_posDOFS+=1;
                        
                        theIElem.getINodeHier(1).sets_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; s_negDOFS+=1;
                        theIElem.getINodeHier(2).sets_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; s_negDOFS+=1;
                    }
                }
            }
        }
        
        if(this.plast){
            if(this.pconstant){
                for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                    theIElem = (IELine2) it.next();
                    theIElem.getINodeHier(1).setp_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    theIElem.getINodeHier(2).setp_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    sumDOFS+=1; p_posDOFS+=1;
                    
                    theIElem.getINodeHier(1).setp_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    theIElem.getINodeHier(2).setp_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    sumDOFS+=1; p_negDOFS+=1;
                }
            }else{
                if(this.pcontiuous){
                    for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                        anInterfaceNode = it.next();
                        if(anInterfaceNode.getNumOfConnectedIElements()>0) {
                            anInterfaceNode.setp_posEFTable(sumDOFS+1+beginDOF); 
                            sumDOFS+=1; p_posDOFS+=1;
                            
                            anInterfaceNode.setp_negEFTable(sumDOFS+1+beginDOF); 
                            sumDOFS+=1; p_negDOFS+=1;
                        }
                    }
                }else{
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        theIElem.getINodeHier(1).setp_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; p_posDOFS+=1;
                        theIElem.getINodeHier(2).setp_posEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; p_posDOFS+=1;
                        
                        theIElem.getINodeHier(1).setp_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; p_negDOFS+=1;
                        theIElem.getINodeHier(2).setp_negEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; p_negDOFS+=1;
                    }
                }
            }
        }
        
        if(this.damage){
            if(this.zconstant){
                for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                    theIElem = (IELine2) it.next();
                    theIElem.getINodeHier(1).setzEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    theIElem.getINodeHier(2).setzEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                    sumDOFS+=1; zDOFS+=1;
                }
            }else{
                if(this.zcontiuous){
                    for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                        anInterfaceNode = it.next();
                        if(anInterfaceNode.getNumOfConnectedIElements()>0) {anInterfaceNode.setzEFTable(sumDOFS+1+beginDOF); sumDOFS+=1; zDOFS+=1;}
                    }
                }else{
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        theIElem.getINodeHier(1).setzEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; zDOFS+=1;
                        theIElem.getINodeHier(2).setzEFTable(sumDOFS+1+beginDOF,theIElem.getID());
                        sumDOFS+=1; zDOFS+=1;
                    }
                }
            }
        }

        // set here the minimization equations
        MinimizationConstraintEquation theMinEq = null;
        MinimizationConstraintTerm aMinTerm;
        ELine theBEMElem;
        int countMinEq=0;
        double xsi;
        double coef;
        boolean check=true;
        int[] ExtraEqOnNondes = null;
        double phi1m=1.; double phi1s=1.;
//        if(this.Domain1.getMaterial().getClass().getName().equalsIgnoreCase("jbem.ViscousMaterial")){
//            if(this.Domain1.isMain()){phi1m=((ViscousMaterial)this.Domain1.getMaterial()).getPhi_1(tau);}else{
//                phi1s=((ViscousMaterial)this.Domain1.getMaterial()).getPhi_1(tau);
//            }
//        }
//        if(Domain2!=null)if(this.Domain2.getMaterial().getClass().getName().equalsIgnoreCase("jbem.ViscousMaterial")){
//            if(this.Domain2.isMain()){phi1m=((ViscousMaterial)this.Domain2.getMaterial()).getPhi_1(tau);}else{
//                phi1s=((ViscousMaterial)this.Domain2.getMaterial()).getPhi_1(tau);
//            }
//        }

        if(this.normal){
            // check for uni-nodal constraints
            int nn,nb=0;
            nn=theInterfaceNodes.size();
            int ne=0, nm=0, nnm=0,nee = 0;
            boolean themain;
            int mctype=0; if(this.SplitNormalDisp)mctype=7;
            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode =  it.next();

                // also on some variables needed for uni-nodal constraints and CHECKING reasons
                if(anInterfaceNode.getBaseID()!=0){
                    nb+=1;
                    if(anInterfaceNode.AskIfInEquation())nee+=1;
                    if(anInterfaceNode.isMain()){nm+=1;}else{nnm+=1;}
                }
                if(anInterfaceNode.AskIfInEquation())ne+=1;
                //////////////////////////////////////////////////////////////////////

                if(anInterfaceNode.getNumOfConnectedIElements()>0){
                    if(NormalRefernce!=null){
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.NormalRefernce);
                        theMinEq.setRigidNormal(rigidnormal);
                        theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }else{
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                        theMinEq.setRigidNormal(rigidnormal);
                        theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);

                    if(anInterfaceNode.getBaseID()!=0){
                        if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                        aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,mctype,coef);
                        theMinEq.addNodeConstraintTerm(aMinTerm);
                        theMinEq.setRigidNormal(rigidnormal);
                        theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }else{
                         theBEMElem = (ELine) anInterfaceNode.getOnElement();
                         if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                         xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                         if(getINodeWithBaseID(theBEMElem.getNodeHier(1).getID())!=null){
                             aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),mctype,coef*theBEMElem.getShapeFunction(1, xsi));
                             theMinEq.addNodeConstraintTerm(aMinTerm);
                             theMinEq.setRigidNormal(rigidnormal);
                             theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                         }
                         if(getINodeWithBaseID(theBEMElem.getNodeHier(2).getID())!=null){
                             aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),mctype,coef*theBEMElem.getShapeFunction(2, xsi));
                             theMinEq.addNodeConstraintTerm(aMinTerm);
                             theMinEq.setRigidNormal(rigidnormal);
                             theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                         }

                    }

                    if(this.getNumOfConnectedDomains()==2){
                        anInterfaceNode = anInterfaceNode.getTwin();
                        if(anInterfaceNode.getBaseID()!=0){
                            if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                            aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,mctype,coef);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }else{
                             theBEMElem = (ELine) anInterfaceNode.getOnElement();
                             if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                             xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                             if(getINodeWithBaseID(theBEMElem.getNodeHier(1).getID())!=null){
                                 aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),mctype,coef*theBEMElem.getShapeFunction(1, xsi));
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                                theMinEq.setRigidNormal(rigidnormal);
                                theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                             }
                             if(getINodeWithBaseID(theBEMElem.getNodeHier(2).getID())!=null){
                                 aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),mctype,coef*theBEMElem.getShapeFunction(2, xsi));
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                                theMinEq.setRigidNormal(rigidnormal);
                                theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                             }

                        }
                    }
                }
            }

            System.out.println("Interface with id = "+this.id);
            System.out.println("============================");
            System.out.println("inodes = "+nn);
            System.out.println("inodes coinsiding with bem node= "+nb);
            System.out.println("inodes coinsiding with bem node and participate to equation= "+ne);
            if(this.closedCurve){
                System.out.println("nodes of elements of interface= "+(this.theInterfaceElements.size()));
            }else{
                System.out.println("nodes of elements of interface= "+(this.theInterfaceElements.size()+1));
            }
            if(nb!=nn){
                if(this.permitNCmesh){
                    System.out.println("Non conforming mesh used.");
                }else{
                    System.err.println("inodes that do not coincide with bem nodes exist, analysis canceled (! ATTENTION: by choise !)");
                    System.exit(steps);
                }
            }
            if(nb!=ne){
            check=false;
            System.out.println("ATTENTION! :bem node exist with not being in equations");
            }

            if(nm>nnm){
                System.out.println("extra equations on main solid");
                themain = true;
            }else{
                System.out.println("extra equations on non main solid");
                themain = false;
            }
            int nen;
            if(this.closedCurve){
                nen =  nee - (this.theInterfaceElements.size());
            }else{
                nen =  nee - (this.theInterfaceElements.size()+1);
            }
            System.out.println("Extra equations needed 'nen'= "+nen+", while already exist: "+this.theMinimizationEquations.size());
            if(nen>0){ExtraEqOnNondes = new int[nen];}else if(nen<0){
                check=false;
                System.out.println("ATTENTION! :non conformed number of equations");
            }else if(nen==0){
                System.out.println("MESSAGE :no extra equations needed");
            }
            int count=0;
            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode = it.next();
                if(anInterfaceNode.isMain() == themain){
                    if(anInterfaceNode.getTwin()!=null)if(anInterfaceNode.getBaseID()!=0 && anInterfaceNode.getTwin().getBaseID()!=0 && count<nen && anInterfaceNode.AskIfInEquation()){
                        ExtraEqOnNondes[count]=anInterfaceNode.getID();
                        count= count+1;
                    }
                }
            }
            if(count<nen){
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                     anInterfaceNode = it.next();
                     if(anInterfaceNode.isMain() == themain){
                        if(anInterfaceNode.getTwin()!=null)if(anInterfaceNode.getBaseID()!=0 && anInterfaceNode.getTwin().getBaseID()==0 && count<nen && anInterfaceNode.AskIfInEquation()){
                            ExtraEqOnNondes[count]=anInterfaceNode.getID();
                            count= count+1;
                        }else if(anInterfaceNode.getBaseID()!=0 && count<nen && anInterfaceNode.AskIfInEquation()){
                            ExtraEqOnNondes[count]=anInterfaceNode.getID();
                            count= count+1;
                        }
                    }
                }
            }
            
            // here the uni-nodal constraints will be defined
            if(ExtraEqOnNondes!=null){
                mctype=0; if(this.SplitNormalDisp)mctype=7;
                for(int i=0; i<ExtraEqOnNondes.length;i++){
                    if(NormalRefernce!=null){
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.NormalRefernce);
                        theMinEq.setRigidNormal(rigidnormal);
                        theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }else{
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                        theMinEq.setRigidNormal(rigidnormal);
                        theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                    anInterfaceNode = this.theInterfaceNodes.get(ExtraEqOnNondes[i]);
                    if(anInterfaceNode.getBaseID()!=0){
                        if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                        aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,mctype,coef);
                        theMinEq.addNodeConstraintTerm(aMinTerm);
                        theMinEq.setRigidNormal(rigidnormal);
                        theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }else{
                         theBEMElem = (ELine) anInterfaceNode.getOnElement();
                         if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                         xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                         aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),mctype,coef*theBEMElem.getShapeFunction(1, xsi));
                         theMinEq.addNodeConstraintTerm(aMinTerm);
                         theMinEq.setRigidNormal(rigidnormal);
                         theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                         aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),mctype,coef*theBEMElem.getShapeFunction(2, xsi));
                         theMinEq.addNodeConstraintTerm(aMinTerm);
                         theMinEq.setRigidNormal(rigidnormal);
                         theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                    }
                }
            }
            
            
            
            if(this.SplitNormalDisp){
                mctype=8;
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                    anInterfaceNode =  it.next();

                    if(anInterfaceNode.getNumOfConnectedIElements()>0){
                        if(NormalRefernce!=null){
                            theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.NormalRefernce);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }else{
                            theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }
                        this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);

                        if(anInterfaceNode.getBaseID()!=0){
                            if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                            aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,mctype,coef);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }else{
                             theBEMElem = (ELine) anInterfaceNode.getOnElement();
                             if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                             xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                             if(getINodeWithBaseID(theBEMElem.getNodeHier(1).getID())!=null){
                                 aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),mctype,coef*theBEMElem.getShapeFunction(1, xsi));
                                 theMinEq.addNodeConstraintTerm(aMinTerm);
                                 theMinEq.setRigidNormal(rigidnormal);
                                 theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                             }
                             if(getINodeWithBaseID(theBEMElem.getNodeHier(2).getID())!=null){
                                 aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),mctype,coef*theBEMElem.getShapeFunction(2, xsi));
                                 theMinEq.addNodeConstraintTerm(aMinTerm);
                                 theMinEq.setRigidNormal(rigidnormal);
                                 theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                             }

                        }

                        if(this.getNumOfConnectedDomains()==2){
                            anInterfaceNode = anInterfaceNode.getTwin();
                            if(anInterfaceNode.getBaseID()!=0){
                                if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                                aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,mctype,coef);
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                                theMinEq.setRigidNormal(rigidnormal);
                                theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                            }else{
                                 theBEMElem = (ELine) anInterfaceNode.getOnElement();
                                 if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                                 xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                                 if(getINodeWithBaseID(theBEMElem.getNodeHier(1).getID())!=null){
                                     aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),mctype,coef*theBEMElem.getShapeFunction(1, xsi));
                                    theMinEq.addNodeConstraintTerm(aMinTerm);
                                    theMinEq.setRigidNormal(rigidnormal);
                                    theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                                 }
                                 if(getINodeWithBaseID(theBEMElem.getNodeHier(2).getID())!=null){
                                     aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),mctype,coef*theBEMElem.getShapeFunction(2, xsi));
                                    theMinEq.addNodeConstraintTerm(aMinTerm);
                                    theMinEq.setRigidNormal(rigidnormal);
                                    theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                                 }

                            }
                        }
                    }
                }
                
                // here the uni-nodal constraints will be defined
                if(ExtraEqOnNondes!=null){
                    for(int i=0; i<ExtraEqOnNondes.length;i++){
                        if(NormalRefernce!=null){
                            theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.NormalRefernce);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }else{
                            theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }
                        this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                        anInterfaceNode = this.theInterfaceNodes.get(ExtraEqOnNondes[i]);
                        if(anInterfaceNode.getBaseID()!=0){
                            if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                            aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,mctype,coef);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                            theMinEq.setRigidNormal(rigidnormal);
                            theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }else{
                             theBEMElem = (ELine) anInterfaceNode.getOnElement();
                             if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                             xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                             aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),mctype,coef*theBEMElem.getShapeFunction(1, xsi));
                             theMinEq.addNodeConstraintTerm(aMinTerm);
                             theMinEq.setRigidNormal(rigidnormal);
                             theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                             aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),mctype,coef*theBEMElem.getShapeFunction(2, xsi));
                             theMinEq.addNodeConstraintTerm(aMinTerm);
                             theMinEq.setRigidNormal(rigidnormal);
                             theMinEq.setReleaseNormalBound(this.ReleaseNormal);
                        }
                    }
                }
                
                
            }
        }
        
        if( tangential && TangentPairs ){
            // check for uni-nodal constraints
            int nn,nb=0;
            nn=theInterfaceNodes.size();
            int ne=0, nm=0, nnm=0,nee = 0;
            boolean themain;

            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode =  it.next();

                // also on some variables needed for uni-nodal constraints and CHECKING reasons
                if(anInterfaceNode.getBaseID()!=0){
                    nb+=1;
                    if(anInterfaceNode.AskIfInEquation())nee+=1;
                    if(anInterfaceNode.isMain()){nm+=1;}else{nnm+=1;}
                }
                if(anInterfaceNode.AskIfInEquation())ne+=1;
                //////////////////////////////////////////////////////////////////////

                if(anInterfaceNode.getNumOfConnectedIElements()>0){
                    if(TangentRefernce!=null){
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.TangentRefernce);
                    }else{
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                    }
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);

                    if(anInterfaceNode.getBaseID()!=0){
                        if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                        aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,1,coef);
                        theMinEq.addNodeConstraintTerm(aMinTerm);
                    }else{
                         theBEMElem = (ELine) anInterfaceNode.getOnElement();
                         if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                         xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                         if(getINodeWithBaseID(theBEMElem.getNodeHier(1).getID())!=null){
                             aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),1,coef*theBEMElem.getShapeFunction(1, xsi));
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                         }
                         if(getINodeWithBaseID(theBEMElem.getNodeHier(2).getID())!=null){
                             aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),1,coef*theBEMElem.getShapeFunction(2, xsi));
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                         }

                    }

                    if(this.getNumOfConnectedDomains()==2){
                        anInterfaceNode = anInterfaceNode.getTwin();
                        if(anInterfaceNode.getBaseID()!=0){
                            if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                            aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,1,coef);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                        }else{
                             theBEMElem = (ELine) anInterfaceNode.getOnElement();
                             if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                             xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                             if(getINodeWithBaseID(theBEMElem.getNodeHier(1).getID())!=null){
                                 aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),1,coef*theBEMElem.getShapeFunction(1, xsi));
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                             }
                             if(getINodeWithBaseID(theBEMElem.getNodeHier(2).getID())!=null){
                                 aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),1,coef*theBEMElem.getShapeFunction(2, xsi));
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                             }

                        }
                    }
                }
            }
            if(nb!=ne){
            check=false;
            }

            if(nm>nnm){
                themain = true;
            }else{
                themain = false;
            }
            int nen;
            if(this.closedCurve){
                nen =  nee - (this.theInterfaceElements.size());
            }else{
                nen =  nee - (this.theInterfaceElements.size()+1);
            }
            if(nen>0){ExtraEqOnNondes = new int[nen];}else if(nen<0){
                check=false;
            }
            int count=0;
            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode = it.next();
                if(anInterfaceNode.isMain() == themain){
                    if(anInterfaceNode.getTwin()!=null)if(anInterfaceNode.getBaseID()!=0 && anInterfaceNode.getTwin().getBaseID()!=0 && count<nen && anInterfaceNode.AskIfInEquation()){
                        ExtraEqOnNondes[count]=anInterfaceNode.getID();
                        count= count+1;
                    }
                }
            }
            if(count<nen){
                for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                     anInterfaceNode = it.next();
                     if(anInterfaceNode.isMain() == themain){
                        if(anInterfaceNode.getTwin()!=null)if(anInterfaceNode.getBaseID()!=0 && anInterfaceNode.getTwin().getBaseID()==0 && count<nen && anInterfaceNode.AskIfInEquation()){
                            ExtraEqOnNondes[count]=anInterfaceNode.getID();
                            count= count+1;
                        }else if(anInterfaceNode.getBaseID()!=0 && count<nen && anInterfaceNode.AskIfInEquation()){
                            ExtraEqOnNondes[count]=anInterfaceNode.getID();
                            count= count+1;
                        }
                    }
                }
            }
            
            // here the uni-nodal constraints will be defined
            if(ExtraEqOnNondes!=null){
                for(int i=0; i<ExtraEqOnNondes.length;i++){
                    if(TangentRefernce!=null){
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.TangentRefernce);
                    }else{
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                    }
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                    anInterfaceNode = this.theInterfaceNodes.get(ExtraEqOnNondes[i]);
                    if(anInterfaceNode.getBaseID()!=0){
                        if(anInterfaceNode.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                        aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,1,coef);
                        theMinEq.addNodeConstraintTerm(aMinTerm);
                    }else{
                         theBEMElem = (ELine) anInterfaceNode.getOnElement();
                         if(theBEMElem.isMain()){coef=1./phi1m;}else{coef=-1./phi1s;}
                         xsi=theBEMElem.getCoordMinDistOfNode(anInterfaceNode)[0];
                         aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(1).getID()),1,coef*theBEMElem.getShapeFunction(1, xsi));
                         theMinEq.addNodeConstraintTerm(aMinTerm);
                         aMinTerm = new MinimizationConstraintTerm(this.getINodeWithBaseID(theBEMElem.getNodeHier(2).getID()),1,coef*theBEMElem.getShapeFunction(2, xsi));
                         theMinEq.addNodeConstraintTerm(aMinTerm);
                    }

                }
            }
        }

        if( tangential && !TangentPairs ){
            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode =  it.next();
                if(anInterfaceNode.getutEFTable()>0){
                    if(TangentRefernce!=null){
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.TangentRefernce);
                    }else{
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                    }
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                    //if(anInterfaceNode.isMain()){coef=1.;}else{coef=-1.;}
                    aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,1);
                    theMinEq.addNodeConstraintTerm(aMinTerm);
                }
            }
            
            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode =  it.next();
                if(anInterfaceNode.getut_posEFTable()>0){
                    theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                    //if(anInterfaceNode.isMain()){coef=1.;}else{coef=-1.;}
                    aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,9);
                    theMinEq.addNodeConstraintTerm(aMinTerm);
                }
            }
            
            for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                anInterfaceNode =  it.next();
                if(anInterfaceNode.getut_negEFTable()>0){
                    theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                    this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                    //if(anInterfaceNode.isMain()){coef=1.;}else{coef=-1.;}
                    aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,10);
                    theMinEq.addNodeConstraintTerm(aMinTerm);
                }
            }
        }
        
        if(this.slip){
            if(this.sconstant){
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        for(int i2=1;i2<=2;i2++){
                            if(SlipRefernce!=null){
                            theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.SlipRefernce);
                            }else{
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            }
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            if(i2==1){aMinTerm = new MinimizationConstraintTerm(3);}else{aMinTerm = new MinimizationConstraintTerm(4);}
                            aMinTerm.setOnElement(theIElem);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                        }
                    }
            }else{
                if(this.scontiuous){
                    for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                        anInterfaceNode = it.next();
                        if(anInterfaceNode.getNumOfConnectedIElements()>0){
                            for(int i2=1;i2<=2;i2++){
                                if(SlipRefernce!=null){
                                    theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.SlipRefernce);
                                }else{
                                    theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                                }
                                this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                                if(i2==1){aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,3);}else{aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,4);}
                                aMinTerm.setOnElement(anInterfaceNode.getSomeElement());
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                            }
                        }
                    }
                }else{
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        for(int i2=1;i2<=2;i2++){
                            anInterfaceNode=theIElem.getINodeHier(1);
                            if(SlipRefernce!=null){
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.SlipRefernce);
                            }else{
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            }
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            if(i2==1){aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,3);}else{aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,4);}
                            aMinTerm.setOnElement(theIElem);
                            theMinEq.addNodeConstraintTerm(aMinTerm);

                            anInterfaceNode=theIElem.getINodeHier(2);
                            if(SlipRefernce!=null){
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.SlipRefernce);
                            }else{
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            }
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            if(i2==1){aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,3);}else{aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,4);}
                            aMinTerm.setOnElement(theIElem);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                        }
                    }
                }
            }
        }
        
        if(this.plast){
            if(this.pconstant){
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        for(int i2=1;i2<=2;i2++){
                            if(PlastRefernce!=null){
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.PlastRefernce);
                            }else{
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            }
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            if(i2==1){aMinTerm = new MinimizationConstraintTerm(5);}else{aMinTerm = new MinimizationConstraintTerm(6);}
                            aMinTerm.setOnElement(theIElem);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                        }
                    }
            }else{
                if(this.pcontiuous){
                    for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                        anInterfaceNode = it.next();
                        if(anInterfaceNode.getNumOfConnectedIElements()>0){
                            for(int i2=1;i2<=2;i2++){
                                if(PlastRefernce!=null){
                                    theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.PlastRefernce);
                                }else{
                                    theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                                }
                                this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                                if(i2==1){aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,5);}else{aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,6);}
                                aMinTerm.setOnElement(anInterfaceNode.getSomeElement());
                                theMinEq.addNodeConstraintTerm(aMinTerm);
                            }
                        }
                    }
                }else{
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        for(int i2=1;i2<=2;i2++){
                            anInterfaceNode=theIElem.getINodeHier(1);
                            if(PlastRefernce!=null){
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.PlastRefernce);
                            }else{
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            }
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            if(i2==1){aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,5);}else{aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,6);}
                            aMinTerm.setOnElement(theIElem);
                            theMinEq.addNodeConstraintTerm(aMinTerm);

                            anInterfaceNode=theIElem.getINodeHier(2);
                            if(PlastRefernce!=null){
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId,this.PlastRefernce);
                            }else{
                                theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            }
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            if(i2==1){aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,5);}else{aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,6);}
                            aMinTerm.setOnElement(theIElem);
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                        }
                    }
                }
            }
        }
        
        if(this.damage){
            if(this.zconstant){
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                        this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                        aMinTerm = new MinimizationConstraintTerm(2);
                        aMinTerm.setOnElement(theIElem);
                        theMinEq.addNodeConstraintTerm(aMinTerm);
                    }
            }else{
                if(this.zcontiuous){
                    for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
                        anInterfaceNode = it.next();
                        if(anInterfaceNode.getNumOfConnectedIElements()>0){
                            theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                            this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                            aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,2);
                            aMinTerm.setOnElement(anInterfaceNode.getSomeElement());
                            theMinEq.addNodeConstraintTerm(aMinTerm);
                        }
                    }
                }else{
                    for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
                        theIElem = (IELine2) it.next();
                        anInterfaceNode=theIElem.getINodeHier(1);
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                        this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                        aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,2);
                        aMinTerm.setOnElement(theIElem);
                        theMinEq.addNodeConstraintTerm(aMinTerm);

                        anInterfaceNode=theIElem.getINodeHier(2);
                        theMinEq = new MinimizationConstraintEquation(++countMinEq+this.beginMinEqId);
                        this.theMinimizationEquations.put(theMinEq.getID(), theMinEq);
                        aMinTerm = new MinimizationConstraintTerm(anInterfaceNode,2);
                        aMinTerm.setOnElement(theIElem);
                        theMinEq.addNodeConstraintTerm(aMinTerm);
                    }
                }
            }
        }
    }

    public int getSumDOFs(){return this.sumDOFS;}
    public int getunDOFs(){return this.unDOFS;}
    public int getun_posDOFs(){return this.un_posDOFS;}
    public int getun_negDOFs(){return this.un_negDOFS;}
    public int getutDOFs(){return this.utDOFS+ut_posDOFS+ut_negDOFS;}
    public int getut_posDOFs(){return ut_posDOFS;}
    public int getut_negDOFs(){return ut_negDOFS;}
    public int getzDOFs(){return this.zDOFS;}
    public int getsDOFs(){return this.s_posDOFS+this.s_negDOFS;}
    public int getpDOFs(){return this.p_posDOFS+this.p_negDOFS;}
    public int getKinematicDOFs(){return (unDOFS+utDOFS+un_posDOFS+un_negDOFS+ut_posDOFS+ut_negDOFS);}
    public int getKinematicNormalDOFs(){return (unDOFS+un_posDOFS);}
    public int getKinematicTangentDOFs(){return (utDOFS);}

    public void PrintMinimizationConstraintEquations(){
         System.out.println("---- Minimization Constraint Equations Of Interface with id = "+this.id+" , number of MCE = "+theMinimizationEquations.size()+"  -----");
                 System.out.println();
                 for(Iterator<MinimizationConstraintEquation>
                         it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                     MinimizationConstraintEquation theConstraintEquation = it.next();
                     theConstraintEquation.Print();
                 }

        System.out.println("-------------------------------------------------------------");
    }
    
    public void PrintMinimizationConstraintEquations(int step){
         System.out.println("---- Minimization Constraint Equations Of Interface with id = "+this.id+" , number of MCE = "+theMinimizationEquations.size()+"  -----");
                 System.out.println();
                 for(Iterator<MinimizationConstraintEquation>
                         it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                     MinimizationConstraintEquation theConstraintEquation = it.next();
                     theConstraintEquation.Print(step);
                 }

        System.out.println("-------------------------------------------------------------");
    }

    public boolean IsDamageIncluded(){return this.damage;}

    public boolean IsSlipIncluded(){return this.slip;}
    
    public boolean IsPlastIncluded(){return this.plast;}
    
    public boolean IsFrictionIncluded(){return this.friction;}

    public boolean IsNormalIncluded(){return this.normal;}

    public boolean IsTangentialIncluded(){return this.tangential;}

    public int setConditionsOnInterface(double[] point, int[] FreedoMap, int step, int fromID){
        // first of all we will make the update of interface nodes using the point and the FreedoMap
//        for(int i=0; i<point.length;i++){
//            System.out.println("point["+i+"]="+point[i]+" FreedoMap["+i+"]="+FreedoMap[i]);
//        }
//        System.out.println("+++++++++++++++++++++++++++++++++");
        int countN=0;
        int countNN=0;
        int countT=0;
        if(this.normal){
            if(!this.SplitNormalDisp){
                for(int i=0; i<FreedoMap.length;i++){
                    if(this.theMinimizationEquations.containsKey(FreedoMap[i]))if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==0)countN+=1;
                }
                if(setMats)this.NMatrix = new AbstractMatrix(countN,countN,0.0);
                this.bMatrix = new AbstractMatrix(countN,1,0.0);
                MinimizationConstraintEquation aMCE;
                //countN = 0;
                for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                    aMCE = it.next();
                    if(aMCE.Variable==0){
                        //bMatrix.set(count, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        bMatrix.set(aMCE.getID()-1-beginMinEqId, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        if(setMats)for(Iterator<MinimizationConstraintTerm> ct=aMCE.getNodeConstraintTerm().values().iterator(); ct.hasNext();){
                            MinimizationConstraintTerm aMCT = ct.next();
                            //NMatrix.set(count, aMCT.getDOF()-1, aMCT.getCoef());
                            NMatrix.set(aMCE.getID()-1-beginMinEqId, aMCT.getDOF()-1-beginDOF, aMCT.getCoef());
                        }
                    //countN+=1;
                    }
                }
                if(setMats){
//                    System.out.println("NMATRIX FOR INTERFACE:"+this.id);
//                    NMatrix.print(6, 6);
                    NMatrix = NMatrix.inverse();
//                    NMatrix.print(6, 6);
                }
//                System.out.println("b normal (pro)");
//                bMatrix.print(12, 12);
                bMatrix=NMatrix.times(bMatrix);
//                System.out.println("b normal (meta)");
//                bMatrix.print(12, 12);
                if(!this.tangential)setMats=false;
                this.update_Normal(step);
            }else{
                for(int i=0; i<FreedoMap.length;i++){
                    if(this.theMinimizationEquations.containsKey(FreedoMap[i]))if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==7)countN+=1;
                }
                if(setMats)this.NMatrix = new AbstractMatrix(countN,countN,0.0);
                this.bMatrix = new AbstractMatrix(countN,1,0.0);
                MinimizationConstraintEquation aMCE;
                //countN = 0;
                for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                    aMCE = it.next();
                    if(aMCE.Variable==7){
                        //bMatrix.set(count, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        bMatrix.set(aMCE.getID()-1, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        if(setMats)for(Iterator<MinimizationConstraintTerm> ct=aMCE.getNodeConstraintTerm().values().iterator(); ct.hasNext();){
                            MinimizationConstraintTerm aMCT = ct.next();
                            //NMatrix.set(count, aMCT.getDOF()-1, aMCT.getCoef());
                            NMatrix.set(aMCE.getID()-1, aMCT.getDOF()-1, aMCT.getCoef());
                        }
                    //countN+=1;
                    }
                }
                if(setMats){
                    //System.out.println("NMATRIX FOR INTERFACE:"+this.id);
                    //NMatrix.print(6, 6);
                    NMatrix = NMatrix.inverse();
                    //NMatrix.print(6, 6);
                }
                bMatrix=NMatrix.times(bMatrix);
//                System.out.println("b pos");
//                bMatrix.print(12, 12);
                if(!this.tangential)setMats=false;
                this.update_Normal_pos(step);
                
                this.bMatrix = new AbstractMatrix(countN,1,0.0);
                //countN = 0;
                for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                    aMCE = it.next();
                    if(aMCE.Variable==8){
                        //bMatrix.set(count, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        bMatrix.set(aMCE.getID()-1-countN, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        if(setMats)for(Iterator<MinimizationConstraintTerm> ct=aMCE.getNodeConstraintTerm().values().iterator(); ct.hasNext();){
                            MinimizationConstraintTerm aMCT = ct.next();
                            //NMatrix.set(count, aMCT.getDOF()-1, aMCT.getCoef());
                            NMatrix.set(aMCE.getID()-1-countN, aMCT.getDOF()-1-countN, aMCT.getCoef());
                        }
                    //countN+=1;
                    }
                }
                bMatrix=NMatrix.times(bMatrix);
//                System.out.println("b neg");
//                bMatrix.print(12, 12);
                this.update_Normal_neg(step);
                this.update_Normal_split(step);
            }
        }


        countT=0;
        if(this.tangential){
            if(!this.SplitTangentialDisp){
                for(int i=0; i<FreedoMap.length;i++){
                    if(this.theMinimizationEquations.containsKey(FreedoMap[i]))if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==1)countT+=1;
                }
                if(setMats)this.TMatrix = new AbstractMatrix(countT,countT,0.0);
                this.bMatrix = new AbstractMatrix(countT,1,0.0);
                MinimizationConstraintEquation aMCE;
                //countT=0;
                countNN=countN; if(this.SplitNormalDisp) countNN=2*countN;
                for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                    aMCE = it.next();
                    if(aMCE.Variable==1){
                        //bMatrix.set(count, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        bMatrix.set(aMCE.getID()-1-countNN-beginMinEqId, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        if(setMats)for(Iterator<MinimizationConstraintTerm> ct=aMCE.getNodeConstraintTerm().values().iterator(); ct.hasNext();){
                            MinimizationConstraintTerm aMCT = ct.next();
                            //TMatrix.set(count, aMCT.getDOF()-1-this.unDOFS, aMCT.getCoef());
                            TMatrix.set(aMCE.getID()-1-countNN-beginMinEqId, aMCT.getDOF()-1-countNN-beginDOF, aMCT.getCoef());
                        }
                        //countT+=1;
                    }
                }
                if(setMats){
                    //System.out.println("TMATRIX");
                    //TMatrix.print(6, 6);
                    TMatrix = TMatrix.inverse();
                    //TMatrix.print(6, 6);
                }
                //System.out.println("b tang (before)");
                //bMatrix.print(6, 6);
                bMatrix=TMatrix.times(bMatrix);
                //System.out.println("b tang (after)");
                //bMatrix.print(6, 6);
                setMats=false;
                this.update_Tangent(step);
            }else{
                countT=0;
                for(int i=0; i<FreedoMap.length;i++){
                    if(this.theMinimizationEquations.containsKey(FreedoMap[i]))if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==9)countT+=1;
                }
                if(setMats)this.TMatrix = new AbstractMatrix(countT,countT,0.0);
                this.bMatrix = new AbstractMatrix(countT,1,0.0);
                MinimizationConstraintEquation aMCE;
                //countT=0;
                countNN=countN; if(this.SplitNormalDisp) countNN=2*countN;
                for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                    aMCE = it.next();
                    if(aMCE.Variable==9){
                        //bMatrix.set(count, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        bMatrix.set(aMCE.getID()-1-countNN, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        if(setMats)for(Iterator<MinimizationConstraintTerm> ct=aMCE.getNodeConstraintTerm().values().iterator(); ct.hasNext();){
                            MinimizationConstraintTerm aMCT = ct.next();
                            //TMatrix.set(count, aMCT.getDOF()-1-this.unDOFS, aMCT.getCoef());
                            TMatrix.set(aMCE.getID()-1-countNN, aMCT.getDOF()-1-countNN, aMCT.getCoef());
                        }
                        //countT+=1;
                    }
                }
                if(setMats){
                    //System.out.println("TMATRIX");
                    //TMatrix.print(6, 6);
                    TMatrix = TMatrix.inverse();
                    //TMatrix.print(6, 6);
                }
                //System.out.println("b tang pos (before)");
                //bMatrix.print(6, 6);
                bMatrix=TMatrix.times(bMatrix);
                //System.out.println("b tang pos (after)");
                //bMatrix.print(6, 6);
                setMats=false;
                this.update_Tangent_pos(step);
                
                // negative
                countT=0;
                for(int i=0; i<FreedoMap.length;i++){
                    if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==10)countT+=1;
                }
                if(setMats)this.TMatrix = new AbstractMatrix(countT,countT,0.0);
                this.bMatrix = new AbstractMatrix(countT,1,0.0);
                //countT=0;
                countNN=countN+countT; if(this.SplitNormalDisp) countNN=2*countN+countT;
                for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                    aMCE = it.next();
                    if(aMCE.Variable==10){
                        //bMatrix.set(count, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        bMatrix.set(aMCE.getID()-1-countNN, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                        if(setMats)for(Iterator<MinimizationConstraintTerm> ct=aMCE.getNodeConstraintTerm().values().iterator(); ct.hasNext();){
                            MinimizationConstraintTerm aMCT = ct.next();
                            //TMatrix.set(count, aMCT.getDOF()-1-this.unDOFS, aMCT.getCoef());
                            TMatrix.set(aMCE.getID()-1-countNN, aMCT.getDOF()-1-countNN, aMCT.getCoef());
                        }
                        //countT+=1;
                    }
                }
                if(setMats){
                    //System.out.println("TMATRIX");
                    //TMatrix.print(6, 6);
                    TMatrix = TMatrix.inverse();
                    //TMatrix.print(6, 6);
                }
                //System.out.println("b tang neg (before)");
                //bMatrix.print(6, 6);
                bMatrix=TMatrix.times(bMatrix);
                //System.out.println("b tang neg (before)");
                //bMatrix.print(6, 6);
                setMats=false;
                this.update_Tangent_neg(step);
            }
            
        }
        
        // end of update
        ConstraintEquation aConstraintEquation;
        ConstraintTerm aConstraintTerm = null;
        ConstraintTermElement aConstraintTermElement=null;
        double[] cvals= new double[1];
        InterfaceNode anInterfaceNode;
        Domain aDomain;

        // set up continuity of tractions 
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            int id1,id2,elemID;
            if(anInterfaceNode.getBase()!=null){
                id1 = anInterfaceNode.getBase().getConnectedElementsIds()[0];
                id2 = anInterfaceNode.getBase().getConnectedElementsIds()[1];
                if(this.Domain1.ExistNodeWithID(anInterfaceNode.getBaseID())){aDomain = Domain1;}else{aDomain = Domain2;}
                if((this.theElements.containsKey(id1) && this.theElements.containsKey(id2))
//                        &&(
//                        AllNodesOfElementBelongToInterface(this.theElements.get(id1))
//                        &&
//                        AllNodesOfElementBelongToInterface(this.theElements.get(id2))
//                        )
                        ){
                    cvals[0]=0.0;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    aConstraintEquation.equalp(aDomain.getElement(id1), aDomain.getElement(id2),
                                    anInterfaceNode.getBaseID(), 1);
                    aDomain.putConstraintEquation(aConstraintEquation);
                    
                    cvals[0]=0.0;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    aConstraintEquation.equalp(aDomain.getElement(id1), aDomain.getElement(id2),
                                    anInterfaceNode.getBaseID(), 2);
                    aDomain.putConstraintEquation(aConstraintEquation);
                }
                
                if(normal && (anInterfaceNode.getunEFTable()!=0||anInterfaceNode.getun_posEFTable()!=0)){
                    //System.out.println("chaos1");
                    //Constraint that gives the projection of displacmenet vector u onto 'normal' vector n
                    cvals[0]=anInterfaceNode.getun()[step];
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 1, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 2, anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }else{
                    //System.out.println("chaos2");
                    //Constraint of zero traction on the normal
                    // pbug (might fixed on 21/9/2011)
                    elemID = id1;
                    if(!this.theElements.containsKey(id1))elemID=id2;
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            1, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            2, anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }

                if(tangential && (anInterfaceNode.getutEFTable()!=0||anInterfaceNode.getut_posEFTable()!=0) && checkEdgeConditions(anInterfaceNode)){
                    //System.out.println("chaos3");
                    //Constraint that gives the projection of displacmenet vector u onto 'tangent' direction of n
                    // pbug
                    cvals[0]=anInterfaceNode.getut()[step];
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 1, -anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 2, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }else{
                    //System.out.println("chaos4");
                    //Constraint of zero traction on the tangent
                    // pbug
                    elemID = id1;
                    if(!this.theElements.containsKey(id1))elemID=id2;
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    
                    
                    double n0,n1;
                    ELine element =(ELine) anInterfaceNode.getBase().getConnectedElement(elemID);
                    // On 14 July of 2013 I deactivated the following 2 lines and activated the lines [NEW]
//                    n0=element.getNormal(anInterfaceNode.getBaseID())[0];
//                    n1=element.getNormal(anInterfaceNode.getBaseID())[1];
                    
                    // [NEW] To define better n0,n1
                    n0=0.5*(((ELine) anInterfaceNode.getBase().getConnectedElement(id1)).getNormal(anInterfaceNode.getBaseID())[0]
                            +((ELine) anInterfaceNode.getBase().getConnectedElement(id2)).getNormal(anInterfaceNode.getBaseID())[0]);
                    n1=0.5*(((ELine) anInterfaceNode.getBase().getConnectedElement(id1)).getNormal(anInterfaceNode.getBaseID())[1]
                            +((ELine) anInterfaceNode.getBase().getConnectedElement(id2)).getNormal(anInterfaceNode.getBaseID())[1]);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            1, -n1);
                            //1, -anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            2, n0);
                            //2, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }
            }
        }

        return fromID;
    }
    
    public int setConditionsOnInterface_hmg(double[] point, int[] FreedoMap, int fromID){
        int n=0;
        if(this.getNumOfConnectedDomains()==1){
            n=setConditionsOnInterface_hmg_sinlge(fromID);
        }else{
            n=setConditionsOnInterface_hmg_multi(point, FreedoMap,fromID);
        }
        return n;
    }
    
    public int setConditionsOnInterface_hmg_multi(double[] point, int[] FreedoMap, int fromID){
        // first of all we will make the update of interface nodes using the point and the FreedoMap
        /*for(int i=0; i<point.length;i++){
            System.out.println("point["+i+"]="+point[i]+" FreedoMap["+i+"]="+FreedoMap[i]);
        }
        System.out.println("+++++++++++++++++++++++++++++++++");*/
        int countN=0;
        int countT=0;
        if(this.normal){
            
            for(int i=0; i<FreedoMap.length;i++){
                if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==0||this.theMinimizationEquations.get(FreedoMap[i]).Variable==7)countN+=1;
            }
            this.bMatrix = new AbstractMatrix(countN,1,0.0);
            MinimizationConstraintEquation aMCE;
            //countN = 0;
            for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                aMCE = it.next();
                if(aMCE.Variable==0||aMCE.Variable==7){
                    bMatrix.set(aMCE.getID()-1, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                }
            }
            bMatrix=NMatrix.times(bMatrix);
            if(!this.tangential)setMats=false;
            this.update_Normalh();
        }


        countT=0;
        if(this.tangential){
            for(int i=0; i<FreedoMap.length;i++){
                if(this.theMinimizationEquations.get(FreedoMap[i]).Variable==1)countT+=1;
            }
            this.bMatrix = new AbstractMatrix(countT,1,0.0);
            MinimizationConstraintEquation aMCE;
            for(Iterator<MinimizationConstraintEquation> it=this.theMinimizationEquations.values().iterator(); it.hasNext();){
                aMCE = it.next();
                if(aMCE.Variable==1){
                    bMatrix.set(aMCE.getID()-1-countN, 0, point[this.findpos(FreedoMap, aMCE.getID())]);
                }
            }
            bMatrix=TMatrix.times(bMatrix);
            this.update_Tangenth();
        }
        
        // end of update
        ConstraintEquation aConstraintEquation;
        ConstraintTerm aConstraintTerm = null;
        ConstraintTermElement aConstraintTermElement=null;
        double[] cvals= new double[1];
        InterfaceNode anInterfaceNode;
        Domain aDomain;

        // set up continuity of tractions 
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            int id1,id2,elemID;
            if(anInterfaceNode.getBase()!=null){
                id1 = anInterfaceNode.getBase().getConnectedElementsIds()[0];
                id2 = anInterfaceNode.getBase().getConnectedElementsIds()[1];
                if(this.Domain1.ExistNodeWithID(anInterfaceNode.getBaseID())){aDomain = Domain1;}else{aDomain = Domain2;}
                if(this.theElements.containsKey(id1) && this.theElements.containsKey(id2)){
                    cvals[0]=0.0;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    aConstraintEquation.equalp(aDomain.getElement(id1), aDomain.getElement(id2),
                                    anInterfaceNode.getBaseID(), 1);
                    aDomain.putConstraintEquation(aConstraintEquation);
                    
                    cvals[0]=0.0;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    aConstraintEquation.equalp(aDomain.getElement(id1), aDomain.getElement(id2),
                                    anInterfaceNode.getBaseID(), 2);
                    aDomain.putConstraintEquation(aConstraintEquation);
                }
                
                if(normal && ((anInterfaceNode.getunEFTable()!=0)||(anInterfaceNode.getunEFTable()!=7))){
                    //System.out.println("chaos1");
                    //Constraint that gives the projection of displacmenet vector u onto 'normal' vector n
                    cvals[0]=anInterfaceNode.getunh();
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 1, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 2, anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }else{
                    //System.out.println("chaos2");
                    //Constraint of zero traction on the normal
                    // pbug (might fixed on 21/9/2011)
                    elemID = id1;
                    if(!this.theElements.containsKey(id1))elemID=id2;
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            1, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            2, anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }

                if(tangential && anInterfaceNode.getutEFTable()!=0 && checkEdgeConditions(anInterfaceNode)){
                    //System.out.println("chaos3");
                    //Constraint that gives the projection of displacmenet vector u onto 'tangent' direction of n
                    // pbug
                    cvals[0]=anInterfaceNode.getuth();
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 1, -anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 2, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }else{
                    //System.out.println("chaos4");
                    //Constraint of zero traction on the tangent
                    // pbug
                    elemID = id1;
                    if(!this.theElements.containsKey(id1))elemID=id2;
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            1, -anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            2, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }
            }
        }

        return fromID;
    }
    
    public int setConditionsOnInterface_hmg_sinlge(int fromID){
        ConstraintEquation aConstraintEquation;
        ConstraintTerm aConstraintTerm = null;
        ConstraintTermElement aConstraintTermElement=null;
        double[] cvals= new double[1];
        InterfaceNode anInterfaceNode;
        Domain aDomain;

        // set up continuity of tractions 
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            int id1,id2,elemID;
            if(anInterfaceNode.getBase()!=null){
                id1 = anInterfaceNode.getBase().getConnectedElementsIds()[0];
                id2 = anInterfaceNode.getBase().getConnectedElementsIds()[1];
                if(this.Domain1.ExistNodeWithID(anInterfaceNode.getBaseID())){aDomain = Domain1;}else{aDomain = Domain2;}
                if(this.theElements.containsKey(id1) && this.theElements.containsKey(id2)){
                    cvals[0]=0.0;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    aConstraintEquation.equalp(aDomain.getElement(id1), aDomain.getElement(id2),
                                    anInterfaceNode.getBaseID(), 1);
                    aDomain.putConstraintEquation(aConstraintEquation);
                    
                    cvals[0]=0.0;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);
                    aConstraintEquation.equalp(aDomain.getElement(id1), aDomain.getElement(id2),
                                    anInterfaceNode.getBaseID(), 2);
                    aDomain.putConstraintEquation(aConstraintEquation);
                }
                
                if(normal && ((anInterfaceNode.getunEFTable()!=0)||(anInterfaceNode.getunEFTable()!=7))){
                    //System.out.println("chaos1");
                    //Constraint that gives the projection of displacmenet vector u onto 'normal' vector n
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 1, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 2, anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }else{
                    //System.out.println("chaos2");
                    //Constraint of zero traction on the normal
                    // pbug (might fixed on 21/9/2011)
                    elemID = id1;
                    if(!this.theElements.containsKey(id1))elemID=id2;
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            1, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            2, anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }

                if(tangential && anInterfaceNode.getutEFTable()!=0 && checkEdgeConditions(anInterfaceNode)){
                    //System.out.println("chaos3");
                    //Constraint that gives the projection of displacmenet vector u onto 'tangent' direction of n
                    // pbug
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 1, -anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aConstraintTerm = new ConstraintTerm(anInterfaceNode.getBase(), 1, 2, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addNodeConstraintTerm(aConstraintTerm);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }else{
                    //System.out.println("chaos4");
                    //Constraint of zero traction on the tangent
                    // pbug
                    elemID = id1;
                    if(!this.theElements.containsKey(id1))elemID=id2;
                    cvals[0]=0.;
                    aConstraintEquation = new ConstraintEquation(++fromID,cvals);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            1, -anInterfaceNode.getNormal()[1]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aConstraintTermElement = new
                            ConstraintTermElement(aDomain.getElement(elemID),
                            anInterfaceNode.getBaseID(),
                            2, anInterfaceNode.getNormal()[0]);
                    aConstraintEquation.addElemConstraintTerm(aConstraintTermElement);

                    aDomain.putConstraintEquation(aConstraintEquation);
                }
            }
        }
        return fromID;
    }

    private int findpos(int[] FVEC, int val){
        int pos=-1;
        for(int i=0;i<FVEC.length;i++){
            if(FVEC[i]==val)pos=i;
        }
        return pos;
    }

    private void update_Normal(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getunEFTable()!=0)anInterfaceNode.updates_un(step, bMatrix.get(anInterfaceNode.getunEFTable()-1-beginDOF, 0));
        }
    }
    
    private void update_Normal_pos(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getun_posEFTable()!=0)anInterfaceNode.updates_un_pos(step, bMatrix.get(anInterfaceNode.getun_posEFTable()-1, 0));
        }
    }
    
    private void update_Normal_neg(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getun_negEFTable()!=0)anInterfaceNode.updates_un_neg(step, bMatrix.get(anInterfaceNode.getun_negEFTable()-1-un_posDOFS, 0));
        }
    }
    
    private void update_Normal_split(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getun_posEFTable()!=0)anInterfaceNode.updates_un_split(step);
        }
    }
    
    private void update_Normalh(){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getunEFTable()!=0)anInterfaceNode.updates_unh(bMatrix.get(anInterfaceNode.getunEFTable()-1, 0));
        }
    }

    private void update_Tangent(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getutEFTable()!=0)anInterfaceNode.updates_ut(step, bMatrix.get(anInterfaceNode.getutEFTable()-1-unDOFS-un_posDOFS-un_negDOFS-beginDOF, 0));
        }
    }
    
    private void update_Tangent_pos(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getut_posEFTable()!=0)anInterfaceNode.updates_ut_pos(step, bMatrix.get(anInterfaceNode.getut_posEFTable()-1-unDOFS-un_posDOFS-un_negDOFS-utDOFS, 0));
        }
    }
    
    private void update_Tangent_neg(int step){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getut_negEFTable()!=0)anInterfaceNode.updates_ut_neg(step, bMatrix.get(anInterfaceNode.getut_negEFTable()-1-unDOFS-un_posDOFS-un_negDOFS-utDOFS-ut_posDOFS, 0));
        }
    }
    
    private void update_Tangenth(){
        InterfaceNode anInterfaceNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            anInterfaceNode = it.next();
            if(anInterfaceNode.getutEFTable()!=0)anInterfaceNode.updates_uth(bMatrix.get(anInterfaceNode.getutEFTable()-1-unDOFS-un_posDOFS-un_negDOFS, 0));
        }
    }

    public void update_Damage(int step, double[] point, int[] FreedoMap){
        InterfaceElement anInterfaceElement;
        int dof;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            anInterfaceElement = it.next();

            dof = anInterfaceElement.getINodeHier(1).getzEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(1).updates_z(step, anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);

            dof = anInterfaceElement.getINodeHier(2).getzEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(2).updates_z(step,anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);
        }
    }
    
    public void update_Slip(int step, double[] point, int[] FreedoMap){
        InterfaceElement anInterfaceElement;
        int dof;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            anInterfaceElement = it.next();

            dof = anInterfaceElement.getINodeHier(1).gets_posEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(1).updates_s_pos(step, anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);

            dof = anInterfaceElement.getINodeHier(2).gets_posEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(2).updates_s_pos(step,anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);
            
            dof = anInterfaceElement.getINodeHier(1).gets_negEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(1).updates_s_neg(step, anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);

            dof = anInterfaceElement.getINodeHier(2).gets_negEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(2).updates_s_neg(step,anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);
        }
    }
    
    public void update_Plast(int step, double[] point, int[] FreedoMap){
        InterfaceElement anInterfaceElement;
        int dof;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            anInterfaceElement = it.next();

            dof = anInterfaceElement.getINodeHier(1).getp_posEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(1).updates_p_pos(step, anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);

            dof = anInterfaceElement.getINodeHier(2).getp_posEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(2).updates_p_pos(step,anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);
            
            dof = anInterfaceElement.getINodeHier(1).getp_negEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(1).updates_p_neg(step, anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);

            dof = anInterfaceElement.getINodeHier(2).getp_negEFTable(anInterfaceElement.getID());
            anInterfaceElement.getINodeHier(2).updates_p_neg(step,anInterfaceElement.getID(), point[this.findpos(FreedoMap, dof)]);
        }
    }

    public double getWorkNormal(int wstep){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkNormal(this,wstep);
        }
        return work;
    }
    
    public double getWorkNormalAUX(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkNormalAUX(this,wstep,tau);
        }
        return work;
    }
    
    public double getWorkNormal(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkNormal(this,wstep,tau);
        }
        return work;
    }
    
    public double getWorkNormalQuad(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkNormalQuad(this,wstep,tau);
        }
        return work;
    }
    
    public double getWorkNormalLinear(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkNormalLinear(this,wstep,tau);
        }
        return work;
    }
    
    public double getWorkNormal_split(int wstep){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkNormal_split(this,wstep);
        }
        return work;
    }

    public double getWorkTangent(int wstep){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkTangent(this,wstep);
        }
        return work;
    }
    
    public double getWorkTangentAUX(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkTangentAUX(this,wstep,tau);
        }
        return work;
    }
    
    public double getWorkTangent(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkTangent(this,wstep, tau);
        }
        return work;
    }
    
    public double getWorkTangentLinear(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkTangentLinear(this,wstep, tau);
        }
        return work;
    }
    
    public double getWorkTangentQuad(int wstep, double tau){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getWorkTangentQuad(this,wstep, tau);
        }
        return work;
    }

    public AbstractMatrix getNMatrix(){return this.NMatrix;}
    
    public AbstractMatrix getTMatrix(){return this.TMatrix;}

    public double getGeneralisedNormalForce(int wstep, InterfaceNode InterNode) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalForce(this, wstep, InterNode.getID());
            }
        }
        return work;
    }

    public double getGeneralisedTangentForce(int wstep, InterfaceNode InterNode) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedTangentForce(this, wstep, InterNode.getID());
            }
        }
        return work;
    }
    
    public double getGeneralisedNormalForce(int wstep, InterfaceNode InterNode, double tau) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalForce(this, wstep, InterNode.getID(), tau);
            }
        }
        return work;
    }
    
    public double getGeneralisedNormalForceQuad(int wstep, InterfaceNode InterNode, double tau) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalForceQuad(this, wstep, InterNode.getID(), tau);
            }
        }
        return work;
    }
    
    public double getGeneralisedNormalForceLinear(int wstep, InterfaceNode InterNode, double tau) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalForceLinear(this, wstep, InterNode.getID(), tau);
            }
        }
        return work;
    }

    public double getGeneralisedTangentForce(int wstep, InterfaceNode InterNode, double tau) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedTangentForce(this, wstep, InterNode.getID(), tau);
            }
        }
        return work;
    }
    
    public double getGeneralisedTangentForceQuad(int wstep, InterfaceNode InterNode, double tau) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedTangentForceQuad(this, wstep, InterNode.getID(), tau);
            }
        }
        return work;
    }
    
    public double getGeneralisedTangentForceLinear(int wstep, InterfaceNode InterNode, double tau) {
        double work=0.;
        int[] elemIDS = null;
        elemIDS = InterNode.getConnectedIElementsIds();
        if(elemIDS.length==0)elemIDS=InterNode.getTwin().getConnectedIElementsIds();
        if(elemIDS!=null){
            for(int i=0;i<elemIDS.length;i++){
                work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedTangentForceLinear(this, wstep, InterNode.getID(), tau);
            }
        }
        return work;
    }

    public double getDissipatedDamageIEnergy(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getDissipatedDamageIEnergy(this,wstep);
        }
        return work;
    }

    public double getDissipatedSlipIEnergy(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getDissipatedSlipIEnergy(this,wstep);
        }
        return work;
    }
    
    public double getDissipatedSlipIEnergy_min(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getDissipatedSlipIEnergy_min(this,wstep);
        }
        return work;
    }
    
    public double getDissipatedPlastIEnergy(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getDissipatedPlastIEnergy(this,wstep);
        }
        return work;
    }
    
    public double getDissipatedPlastIEnergy_min(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getDissipatedPlastIEnergy_min(this,wstep);
        }
        return work;
    }
    
    public double getDissipatedFrictionIEnergy(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getDissipatedFrictionIEnergy(this,wstep);
        }
        return work;
    }
    
    public double getIntegralDamageDrivingForce(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getIntegralDamageDrivingForce(this,wstep);
        }
        return work;
    }
    
    public double getIntegralSlipDrivingForce(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getIntegralSlipDrivingForce(this,wstep);
        }
        return work;
    }
    
    public double getIntegralPlastDrivingForce(int wstep) {
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getIntegralPlastDrivingForce(this,wstep);
        }
        return work;
    }
    

    public double getDissipatedDamageIForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedDamageIForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getDissipatedDamageIConstant(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedDamageIConstant(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }

    public double getGeneralisedNormalDrivingForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalDrivingForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getGeneralisedNormalDrivingForce(int wstep, InterfaceNode InterNode, int dof, double tau) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalDrivingForce(this, wstep, InterNode.getID(), dof, tau);
        }
        return work;
    }
    
    public double getGeneralisedNormalDrivingForce_split(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedNormalDrivingForce_split(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }

    public double getGeneralisedTangentDrivingForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedTangentDrivingForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getGeneralisedTangentDrivingForce(int wstep, InterfaceNode InterNode, int dof, double tau) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedTangentDrivingForce(this, wstep, InterNode.getID(), dof, tau);
        }
        return work;
    }

    public double getDissipatedSlipIForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedSlipIForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getDissipatedFrictionIForce(int wstep, InterfaceNode InterNode) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedFrictionIForce(this, wstep, InterNode.getID());
        }
        return work;
    }
    
    public double getDissipatedPlastIForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedPlastIForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getDissipatedPlastIConstant(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedPlastIConstant(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getDissipatedSlipIConstant(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getDissipatedSlipIConstant(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }

    public double getGeneralisedSlipDrivingForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedSlipDrivingForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }
    
    public double getGeneralisedPlastDrivingForce(int wstep, InterfaceNode InterNode, int dof) {
        double work=0.;
        int[] elemIDS;
        elemIDS = InterNode.getConnectedIElementsIds();
        for(int i=0;i<elemIDS.length;i++){
            work+=this.theInterfaceElements.get(elemIDS[i]).getGeneralisedPlastDrivingForce(this, wstep, InterNode.getID(), dof);
        }
        return work;
    }

    public int getNumOfConnectedDomains(){
        int num=0;
        if(this.Domain1!=null)num+=1;
        if(this.Domain2!=null)num+=1;
        return num;
    }

    public int getNumOfDamagedElements(int step){
        int n=0;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            if( (elem.getINodeHier(1).getz(elem.getID())[step]<1.) || (elem.getINodeHier(2).getz(elem.getID())[step]<1.))n+=1;
        }
        return n;
    }
    
    public int getNumOfTotalDamagedElements(int step){
        int n=0;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            if( (elem.getINodeHier(1).getz(elem.getID())[step]<1.e-8) && (elem.getINodeHier(2).getz(elem.getID())[step]<1.e-8))n+=1;
        }
        return n;
    }

    public int getNumElements(){
        return this.theElements.size();
    }

    public int getNumIElements(){
        return this.theInterfaceElements.size();
    }

    public int getDomain1_id(){
        int ret=0;
        if(this.Domain1!=null) ret =this.Domain1.getID();
        return ret;
    }
    
    public Domain getDomain1(){
        return this.Domain1;
    }

    public int getDomain2_id(){
        int ret=0;
        if(this.Domain2!=null) ret =this.Domain2.getID();
        return ret;
    }
    
    public Domain getDomain2(){
        return this.Domain2;
    }
    
    public void setSlipReference(double val){this.SlipRefernce=val;}
    
    public void setPlastReference(double val){this.PlastRefernce=val;}
    
    public void setNormalReference(double val){this.NormalRefernce=val;}
    
    public void setTangentReference(double val){this.TangentRefernce=val;}
    
    public void setTangentPairs(boolean t){this.TangentPairs=t;}
    
    public void setRigidNormal(){this.rigidnormal=true;}
    
    public void setZeroEndTangentTraction(){this.zeroEndTangentTraction=true;}
    
    public void setZeroStartTangentTraction(){this.zeroStartTangentTraction=true;}
    
    public void setInitialDamage_x(double lowX, double upperX, double val){
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement theElement = it.next();
            for(Iterator<InterfaceNode>nt=theElement.getNodes().values().iterator();nt.hasNext();){
                InterfaceNode aNode=nt.next();
                if(aNode.getCoordinates()[0]>=lowX && aNode.getCoordinates()[0]<=upperX){
                    if(this.zconstant){
                        theElement.getINodeHier(1).updates_z(0, theElement.getID(), val);
                        theElement.getINodeHier(2).updates_z(0, theElement.getID(), val);
                    }else{
                        aNode.updates_z(0, theElement.getID(), val);
                    }
                }
            }
        }
    }
    
    
    private boolean InSameDomain(InterfaceNode aNode, Element anElement){
        boolean answer=false;
        if(aNode.getConnectedIElementsIds().length>0)
        if(aNode.getBase()!=null){
            if(this.Domain1.getNodes().containsKey(aNode.getBaseID())){
                if(this.Domain1.getElements().containsKey(anElement.getID()))answer=true;
            }else{
                if(this.Domain2.getElements().containsKey(anElement.getID()))answer=true;
            }
        }else{
            if(this.Domain1.getElements().containsKey(aNode.getConnectedIElementsIds()[0])){
                if(this.Domain1.getElements().containsKey(anElement.getID()))answer=true;
            }else{
                if(this.Domain2.getElements().containsKey(anElement.getID()))answer=true;
            }
        }
        
        return answer;
    }
    
    private boolean checkEdgeConditions(InterfaceNode anInterfaceNode){
        boolean ans=true;
        // for ans=true the tangential displacement conditions will be imposed
        if( this.zeroStartTangentTraction && (anInterfaceNode.getID() == this.INodeIDstart || anInterfaceNode.getTwinID() ==this.INodeIDstart) )ans=false;
        if( this.zeroEndTangentTraction   && (anInterfaceNode.getID() == this.INodeIDend   || anInterfaceNode.getTwinID() ==this.INodeIDend) )ans=false;
        return ans;
    }
    
    public void setClosedCurve(){this.closedCurve=true;}
    
    public boolean getClosedCurve(){return this.closedCurve;}
    
    public void setReleaseNormal(){this.ReleaseNormal=true;}
    
    public void setSplitNormalDisp(){this.SplitNormalDisp=true;}
    
    public void setSplitTangentiallDisp(){this.SplitTangentialDisp=true;}
    
    public boolean getSplitTangentialDisp(){return this.SplitTangentialDisp;}
    
    public void setPermitNonConformingMesh(boolean w){this.permitNCmesh=w;}
    
    public double getNewCrackLength(int wstep){
        double work=0.;
        for(Iterator<InterfaceElement> it=this.theInterfaceElements.values().iterator(); it.hasNext();){
            InterfaceElement elem = it.next();
            work+=elem.getNewCrack(wstep);
        }
        return work;
    }
    
    public double getNormalBEMTraction(int IElid, int INodeid, int wstep){
        double val=0.;
        int theElid,baseID,theDomID;
        int count=0;
        double nx,ny;
        if(this.theInterfaceNodes.get(INodeid).getBase()!=null){
            baseID=this.theInterfaceNodes.get(INodeid).getBase().getID();
            theDomID=this.Domain1.getID();
            if(this.Domain2!=null){if(this.Domain2.ExistNodeWithID(baseID)){theDomID=this.Domain2.getID();}}
            for(int i=0;i<theInterfaceNodes.get(INodeid).getBase().getConnectedElementsIds().length;i++){
                theElid = theInterfaceNodes.get(INodeid).getBase().getConnectedElementsIds()[i];
                if(theElements.containsKey(theElid)){
                    nx=((ELine2) theElements.get(theElid)).getNormal(baseID)[0]; ny=((ELine2) theElements.get(theElid)).getNormal(baseID)[1]; 
                    count+=1;
                    val+=-theElements.get(theElid).getNode(baseID).getp(theElements.get(theElid), Domain1.getFundamentalSolution().get_p_DOFs(),3)[0][wstep]*ny
                            +theElements.get(theElid).getNode(baseID).getp(theElements.get(theElid), Domain1.getFundamentalSolution().get_p_DOFs(),3)[1][wstep]*nx;
                }
            }
            val=val/count;
        }else{
            throw new UnsupportedOperationException("Not supported yet.");
        }
        return val;
    }

    public void setZeroSolution(int itime) {
        InterfaceNode theNode;
        for(Iterator<InterfaceNode> it=this.theInterfaceNodes.values().iterator(); it.hasNext();){
            theNode = it.next();
            theNode.setZeroResponse(itime);
        }
    }
    
    public void setPlasticIncrementForINodes(boolean what){
        InterfaceNode anInterfaceNode = new InterfaceNode();
        anInterfaceNode.setPlasticIncrement(what);
        this.PlasticIncrement=what;
    }
    
    public boolean getPlasticIncrement(){return this.PlasticIncrement;}
    
//    private boolean AllNodesOfElementBelongToInterface(Element theElement){
//        boolean bvar=true;
//        Node theNode;
//        for(Iterator<Node> it=theElement.getNodes().values().iterator(); it.hasNext();){
//            theNode = it.next();
//            if(!this.theNodes.containsKey(theNode.getID()))bvar=false;
//        }
//        return bvar;
//    }
    
}