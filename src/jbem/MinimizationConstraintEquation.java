/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author pchr
 */
public class MinimizationConstraintEquation {
    public static boolean stcInitPermitDamage=true;
    protected int id;
    protected Map<Integer,MinimizationConstraintTerm> TermsOfConstraintEquation = new TreeMap<Integer,MinimizationConstraintTerm>();
    protected int numOfTerms=0;
    protected int Variable;
    // Variable is 
    // 0 for normal displ 
    // 1 for tangent 
    // 2 for damage 
    // 3-4 for slip 
    // 5-6 for plast 
    // 7 normal_pos
    // 8 normal_neg
    // 9 tangent_pos
    // 10 tangent_neg
    // 11 damage_pos
    // 12 damage_neg
    protected Double reference;
    protected static boolean multi=false;
    protected boolean RigidNormal=false;
    protected boolean ReleaseNormalBound=false;
    public boolean PermitDamage;

    public MinimizationConstraintEquation(int id){
        PermitDamage=stcInitPermitDamage;
        this.id=id;
        reference=null;
    }

    public MinimizationConstraintEquation(int id, double reference){
        PermitDamage=stcInitPermitDamage;
        this.id=id;
        this.reference=reference;
    }

    public void addNodeConstraintTerm(MinimizationConstraintTerm cterm){
        if(TermsOfConstraintEquation.size()>0){
            for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                MinimizationConstraintTerm ExistMinTerms = it.next();
                if(ExistMinTerms.getVariable()!=this.Variable){this.Variable = 5;}
            }
        }else{
            Variable = cterm.getVariable();
        }
        ++this.numOfTerms;
        this.TermsOfConstraintEquation.put(numOfTerms, cterm);
    }

    public int getID(){
        return this.id;
    }

    public Map<Integer,MinimizationConstraintTerm> getNodeConstraintTerm(){
        return this.TermsOfConstraintEquation;
    }

    void Print() {
        System.out.println("Minimization Costraint equation id : "+this.id+", num of NTerms = "+this.numOfTerms+" || LowerB = "+this.getLower(0)+" || UpperB = "+this.getUpper(0));
         for(Iterator<MinimizationConstraintTerm> nit=TermsOfConstraintEquation.values().iterator() ; nit.hasNext();){
             MinimizationConstraintTerm NodeConstraintTerm = nit.next();
             NodeConstraintTerm.Print();
         }
         System.out.println();
    }
    
    void Print(int step) {
        System.out.println("Minimization Costraint equation id : "+this.id+", num of NTerms = "+this.numOfTerms+" || LowerB = "+this.getLower(step)+" || UpperB = "+this.getUpper(step));
         for(Iterator<MinimizationConstraintTerm> nit=TermsOfConstraintEquation.values().iterator() ; nit.hasNext();){
             MinimizationConstraintTerm NodeConstraintTerm = nit.next();
             NodeConstraintTerm.Print();
         }
         System.out.println();
    }

    public Double getLower(int step){
        Double val = null;
        switch(Variable){
            case 0: val = getLowerNormal(step); break;
            case 1: val = getLowerTangent(step); break;
            case 2: val = getLowerDamage(step); break;
            case 3: val = getLowerSlip(step); break;
            case 4: val = getLowerSlip(step); break;
            case 5: val = getLowerPlast(step); break;
            case 6: val = getLowerPlast(step); break;
            case 7: val = getLowerNormal_pos(step); break;
            case 8: val = getLowerNormal_neg(step); break;
            case 9: val = getLowerTangent_pos(step); break;
            case 10: val = getLowerTangent_neg(step); break;    
            default:
                System.err.println("mixed Variables in minimization constraint equation, type = "+Variable);
                System.exit(step);
                break;
        }
        return val;
    }

    public Double getUpper(int step){
        Double val = null;
        switch(Variable){
            case 0: val = getUpperNormal(step); break;
            case 1: val = getUpperTangent(step); break;
            case 2: val = getUpperDamage(step); break;
            case 3: val = getUpperSlip(step); break;
            case 4: val = getUpperSlip(step); break;
            case 5: val = getUpperPlast(step); break;
            case 6: val = getUpperPlast(step); break;
            case 7: val = getUpperNormal_pos(step); break;
            case 8: val = getUpperNormal_neg(step); break;
            case 9: val = getUpperTangent_pos(step); break;
            case 10: val = getUpperTangent_neg(step); break; 
            default:
                System.err.println("mixed Variables in minimization constraint equation, type = "+Variable);
                System.exit(step);
                break;
        }
        return val;
    }

    public int getVariable(){return this.Variable;}

    public Double getLowerNormal(int step){
        Double val = null;
        if(reference!=null)val = -this.reference;
        if(this.TermsOfConstraintEquation.size()>1 && this.RigidNormal){
            val=0.0;
        }
        return val;
    }

    public Double getUpperNormal(int step){
        Double val = null;
        if(this.TermsOfConstraintEquation.size()==1){
            if(MinimizationConstraintEquation.multi){
                if(reference!=null)val = this.reference;
            }else{
                //val=0.;
                for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                        MinimizationConstraintTerm ExistMinTerms = it.next();
                        val = ExistMinTerms.getINode().getCoordinates()[1];
                        //val =0.;
//                        double px,py;
//                        px=1.5; py=0.75;
//                        double phi;
//                        double xn,yn;
//                        xn=ExistMinTerms.getINode().getCoordinates()[0];
//                        yn=ExistMinTerms.getINode().getCoordinates()[1];
//                        phi = Math.atan(Math.abs(xn-px)/Math.abs(yn-py));
//                        val = Math.cos(phi)*ExistMinTerms.getINode().getCoordinates()[1];
                }
            }
        }else{
            if(this.RigidNormal){
                val=0.;
            }else{
                for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                        MinimizationConstraintTerm ExistMinTerms = it.next();
                        val = ExistMinTerms.getINode().getGap();
                }
            }
            
        }
        if(this.ReleaseNormalBound)val=null;
        return val;
    }
    
    public Double getBoundVV(int step, ViscousMaterial Mat_main, ViscousMaterial Mat_sec, double tau){
        /**
         * Implemented for the vanishing viscosity procedure algorithm.
         */
        Double val = null;
//        if(this.TermsOfConstraintEquation.size()==1){
//            double phi1=Mat_main.getPhi_1(tau);
//            double phi2=Mat_main.getPhi_2(tau);
//            double phi3=Mat_main.getPhi_3(tau);
//            for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
//                    MinimizationConstraintTerm ExistMinTerms = it.next();
//                    if(step==0){ 
//                        val = ExistMinTerms.getINode().getCoordinates()[1]*phi1;
//                    }else if(step==1){
//                        val = ExistMinTerms.getINode().getCoordinates()[1]*phi1+ExistMinTerms.getINode().getun_(Mat_main,tau)[step-1]*phi2;
//                    }else{
//                        val = ExistMinTerms.getINode().getCoordinates()[1]*phi1+ExistMinTerms.getINode().getun_(Mat_main,tau)[step-1]*phi2+ExistMinTerms.getINode().getun_(Mat_main,tau)[step-2]*phi3;
//                    }
//            }
//        }else{
//            double phi1m=Mat_main.getPhi_1(tau);
//            double phi2m=Mat_main.getPhi_2(tau);
//            double phi3m=Mat_main.getPhi_3(tau);
//            double phi1s=Mat_sec.getPhi_1(tau);
//            double phi2s=Mat_sec.getPhi_2(tau);
//            double phi3s=Mat_sec.getPhi_3(tau);
//            MinimizationConstraintTerm ExistMinTerms =TermsOfConstraintEquation.values().iterator().next();
//            if(ExistMinTerms.getINode().isMain()){
//                if(step==0){ 
//                    val = ExistMinTerms.getINode().getGap();
//                }else if(step==1){
//                    val = ExistMinTerms.getINode().getGap()+ExistMinTerms.getINode().getun_(Mat_main,tau)[step-1]*phi2m
//                            +ExistMinTerms.getINode().getTwin().getun_(Mat_sec,tau)[step-1]*phi2s;
//                }else{
//                    val = ExistMinTerms.getINode().getGap()+ExistMinTerms.getINode().getun_(Mat_main,tau)[step-1]*phi2m+ExistMinTerms.getINode().getun_(Mat_main,tau)[step-2]*phi3m
//                            +ExistMinTerms.getINode().getTwin().getun_(Mat_sec,tau)[step-1]*phi2s+ExistMinTerms.getINode().getTwin().getun_(Mat_sec,tau)[step-2]*phi3s;
//                }
//            }else{
//                if(step==0){ 
//                    val = ExistMinTerms.getINode().getGap();
//                }else if(step==1){
//                    val = ExistMinTerms.getINode().getGap()+ExistMinTerms.getINode().getun_(Mat_sec,tau)[step-1]*phi2s
//                            +ExistMinTerms.getINode().getTwin().getun_(Mat_main,tau)[step-1]*phi2m;
//                }else{
//                    val = ExistMinTerms.getINode().getGap()+ExistMinTerms.getINode().getun_(Mat_sec,tau)[step-1]*phi2s+ExistMinTerms.getINode().getun_(Mat_sec,tau)[step-2]*phi3s
//                            +ExistMinTerms.getINode().getTwin().getun_(Mat_main,tau)[step-1]*phi2m+ExistMinTerms.getINode().getTwin().getun_(Mat_main,tau)[step-2]*phi3m;
//                }
//            }
//        }
        double nu=Mat_main.getDispRelaxationTime_1();
        double val1;
        if(this.TermsOfConstraintEquation.size()==1){
            for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                    MinimizationConstraintTerm ExistMinTerms = it.next();
                    double px,py;
                    px=1.5; py=0.75;
                    double phi;
                    double xn,yn;
                    xn=ExistMinTerms.getINode().getCoordinates()[0];
                    yn=ExistMinTerms.getINode().getCoordinates()[1];
                    phi = Math.atan(Math.abs(xn-px)/Math.abs(yn-py));
                    val1 = Math.cos(phi)*ExistMinTerms.getINode().getCoordinates()[1];
                    if(step==0){
                        val = 0.;
//                        val = val1*(1.+nu/tau)-(nu/tau)*ExistMinTerms.getINode().getun_(Mat_main,tau)[step];
                    }else{
                        val = val1*(1.+nu/tau)-(nu/tau)*ExistMinTerms.getINode().getun_(Mat_main,tau)[step-1];
                    }
            }
        }else{
            throw new UnsupportedOperationException("Not supported yet.");
        }
        
        
        return val;
    }
    
    public Double getLowerNormal_pos(int step){
        Double val = null;
        if(reference!=null)val = -this.reference;
        if(MinimizationConstraintEquation.multi){
            if(this.TermsOfConstraintEquation.size()>1){
                val=0.0;
            }
        }else{
            val=0.0;
        }
        
        return val;
    }
    
    public Double getUpperNormal_pos(int step){
        Double val = null;
        if(this.TermsOfConstraintEquation.size()==1){
            if(MinimizationConstraintEquation.multi){
                if(reference!=null)val = this.reference;
            }else{val=0.;}
        }else{
            if(this.RigidNormal){
                val=0.;
            }else{
                for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                        MinimizationConstraintTerm ExistMinTerms = it.next();
                        val = ExistMinTerms.getINode().getGap();
                }
            }
            
        }
        if(this.ReleaseNormalBound)val=null;
        return val;
    }
    
    public Double getLowerNormal_neg(int step){
        Double val = null;
        return val;
    }
    
    public Double getUpperNormal_neg(int step){
        Double val = null;
        if(reference!=null)val = -this.reference;
        if(MinimizationConstraintEquation.multi){
            if(this.TermsOfConstraintEquation.size()>1){
                val=0.0;
            }
        }else{
            val=0.0;
        }
        return val;
    }

    public Double getLowerTangent(int step){
        Double val = null;
        if(reference!=null)val = -this.reference;
        return val;
    }

    public Double getUpperTangent(int step){
        Double val = null;
        if(reference!=null)val = this.reference;
        return val;
    }
    
    public Double getLowerTangent_pos(int step){
        return 0.;
    }
    
    public Double getLowerTangent_neg(int step){
        return 0.;
    }
    
    public Double getUpperTangent_pos(int step){
        return null;
    }
    
    public Double getUpperTangent_neg(int step){
        return null;
    }

    public Double getLowerDamage(int step){
        Double val = null;
        val = 1.0;
        if(this.PermitDamage)val=0.0;
        return val;
    }

    public Double getUpperDamage(int step){
        Double val = null;
        if(this.TermsOfConstraintEquation.size()>1){
            System.err.println("minimization constraint equation on damage variable with greater than one, terms.");
            System.exit(step);
        }else{
            if(step>0){
                for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                    MinimizationConstraintTerm ExistMinTerms = it.next();
                    if(ExistMinTerms.getINode()!=null){
                        val = ExistMinTerms.getINode().getz(ExistMinTerms.getIElement().getID())[step-1];
                    }else{
                        val = ExistMinTerms.getIElement().getINodeHier(1).getz(ExistMinTerms.getIElement().getID())[step-1];
                    }
                }
            }else{
                //val = 1.0;
                for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                    MinimizationConstraintTerm ExistMinTerms = it.next();
                    if(ExistMinTerms.getINode()!=null){
                        val = ExistMinTerms.getINode().getz(ExistMinTerms.getIElement().getID())[0];
                    }else{
                        val = ExistMinTerms.getIElement().getINodeHier(1).getz(ExistMinTerms.getIElement().getID())[0];
                    }
                }
            }
        }
        return val;
    }

    public Double getLowerSlip(int step){
        Double val = null;
        if(reference!=null)val = -this.reference;
        return val;
    }

    public Double getUpperSlip(int step){
        Double val = null;
        if(reference!=null)val = this.reference;
        return val;
    }
    
    public Double getLowerPlast(int step){
        Double val = null;
        if(reference!=null)val = -this.reference;
        return val;
    }

    public Double getUpperPlast(int step){
        Double val = null;
        if(reference!=null) val = this.reference;
        return val;
    }

    public int getNumOfTerms(){return this.numOfTerms;}

    public MinimizationConstraintTerm getTheUniqueMCTerm(){
        if(this.getNumOfTerms()>1)System.err.println("ATTENTION ! getTheUniqueMCTerm has been utilized to a multi term Minimization constraint equation");
        return this.TermsOfConstraintEquation.values().iterator().next();
    }
    
    public MinimizationConstraintTerm getTheFirstMCTerm(){
        return this.TermsOfConstraintEquation.values().iterator().next();
    }
    
    public boolean ExistTermWithDOF(int dof){
        boolean exist = false;
        for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
            MinimizationConstraintTerm ExistMinTerms = it.next();
            if(ExistMinTerms.getDOF()==dof)exist=true;
        }
        return exist;
    }
    
    public static void setMulti(){
        MinimizationConstraintEquation.multi=true;
    }
    
    public boolean isINodesParticipateDamaged(int tstep){
        boolean answer =true;
        for(Iterator<MinimizationConstraintTerm> it=this.TermsOfConstraintEquation.values().iterator(); it.hasNext();){
                MinimizationConstraintTerm ExistMinTerms = it.next();
                if(!ExistMinTerms.getINode().AskIfCollapsed(tstep))answer = false;
        }
        return answer;
    }
    
    public void setReference(Double ref){
        this.reference=ref;
    }
    
    public void setRigidNormal(boolean val){RigidNormal=val;}
    
    public void setReleaseNormalBound(boolean val){this.ReleaseNormalBound=val;}

}
