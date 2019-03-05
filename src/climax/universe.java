/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package climax;

import edu.uta.futureye.core.Mesh;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import jdem.DEMdomain;

/**
 *
 * @author pchr
 */
public class universe {
    protected static Map<Integer,DEMdomain> theDEMDomains = new TreeMap<Integer,DEMdomain>();
    protected static Map<Integer,jbem.Domain> theBEMDomains = new TreeMap<Integer,jbem.Domain>();
    protected static Map<Integer,jfem.Domain> theFEMDomains = new TreeMap<Integer,jfem.Domain>();
    protected static Map<Integer,gendomain.Domain> theGENDomains = new TreeMap<Integer,gendomain.Domain>();
    protected static Map<Integer,jpde> thePDEDomains = new TreeMap<Integer,jpde>();
    protected static Map<Integer,contraption> theContraptions = new TreeMap<Integer,contraption>();
    private static int numofPDEDomains=0;

    // constructor
    public universe(){}
    
    public static DEMdomain DEMDomain(){
        DEMdomain d = new DEMdomain();
        putDEMDomain(d);
        return d;
    }
    
    public static jbem.Domain BEMDomain(){
        jbem.Domain d = new jbem.Domain();
        putBEMDomain(d);
        return d;
    }
    
    public static jfem.Domain FEMDomain(){
        jfem.Domain d = new jfem.Domain();
        putFEMDomain(d);
        return d;
    }
    
    public static gendomain.Domain GENDomain(){
        gendomain.Domain d = new gendomain.Domain();
        putGENDomain(d);
        return d;
    }
    
    public static void putFEMDomain(jfem.Domain aDomain){
        theFEMDomains.put(aDomain.getID(), aDomain);
    }
    
    public static void putBEMDomain(jbem.Domain aDomain){
        theBEMDomains.put(aDomain.getID(), aDomain);
    }
    
    public static void putDEMDomain(DEMdomain aDomain){
        theDEMDomains.put(aDomain.getID(), aDomain);
    }
    
    public static void putGENDomain(gendomain.Domain aDomain){
        theGENDomains.put(aDomain.getID(), aDomain);
    }
    
    public static void putPDEDomain(jpde aDomain){
        thePDEDomains.put(++numofPDEDomains, aDomain);
    }
    
    public static void putContraption(contraption Contraption){
        theContraptions.put(Contraption.getID(), Contraption);
    }
    
    public jbem.Domain getBEMDomain(int id){
        return universe.theBEMDomains.get(id);
    }
    
    public DEMdomain getDEMDomain(int id){
        return universe.theDEMDomains.get(id);
    }
    
    public jfem.Domain getFEMDomain(int id){
        return universe.theFEMDomains.get(id);
    }
    
    public gendomain.Domain getGENDomain(int id){
        return universe.theGENDomains.get(id);
    }
    
    public jpde getPDEDomain(int id){
        return universe.thePDEDomains.get(id);
    }
    
    public contraption getContraption(int id){
        return universe.theContraptions.get(id);
    }
    
    public int getBEMDomainsNum(){
        return universe.theBEMDomains.size();
    }
    
    public int[] getBEMDomainsIDs(){
        int[] Ids = new int[this.theBEMDomains.size()];
        int i=0;
        for(Iterator<jbem.Domain> it=this.theBEMDomains.values().iterator(); it.hasNext();){
            jbem.Domain theDomain = it.next();
            Ids[i]=theDomain.getID();
            i++;
        }
        return Ids;
    }
    
    public int getDEMDomainsNum(){
        return universe.theDEMDomains.size();
    }
    
    public int[] getDEMDomainsIDs(){
        int[] Ids = new int[this.theDEMDomains.size()];
        int i=0;
        for(Iterator<DEMdomain> it=this.theDEMDomains.values().iterator(); it.hasNext();){
            DEMdomain theDomain = it.next();
            Ids[i]=theDomain.getID();
            i++;
        }
        return Ids;
    }
    public int getFEMDomainsNum(){
        return universe.theFEMDomains.size();
    }
    
    public int[] getFEMDomainsIDs(){
        int[] Ids = new int[this.theFEMDomains.size()];
        int i=0;
        for(Iterator<jfem.Domain> it=this.theFEMDomains.values().iterator(); it.hasNext();){
            jfem.Domain theDomain = it.next();
            Ids[i]=theDomain.getID();
            i++;
        }
        return Ids;
    }
    
    public int getGENDomainsNum(){
        return universe.theGENDomains.size();
    }
    
    public int[] getGENDomainsIDs(){
        int[] Ids = new int[this.theGENDomains.size()];
        int i=0;
        for(Iterator<gendomain.Domain> it=this.theGENDomains.values().iterator(); it.hasNext();){
            gendomain.Domain theDomain = it.next();
            Ids[i]=theDomain.getID();
            i++;
        }
        return Ids;
    }
    
    public int getPDEDomainsNum(){
        return universe.thePDEDomains.size();
    }
    
    public int[] getPDEDomainsIDs(){
        int[] Ids = new int[this.thePDEDomains.size()];
        int i=0;
        for (Map.Entry<Integer,jpde> entry : this.thePDEDomains.entrySet()) {
                int key = entry.getKey();
                jpde theDomain = entry.getValue();
                Ids[i]=key;
                i++;
            }
        return Ids;
    }
    
    public int getContraptionsNum(){
        return universe.theContraptions.size();
    }
    
    public int[] getContraptionsIDs(){
        int[] Ids = new int[this.theContraptions.size()];
        int i=0;
        for(Iterator<contraption> it=this.theContraptions.values().iterator(); it.hasNext();){
            contraption theDomain = it.next();
            Ids[i]=theDomain.getID();
            i++;
        }
        return Ids;
    }
    
    public Map getDEMDomains(){
        return universe.theDEMDomains;
    }
    
    public Map getBEMDomains(){
        return universe.theBEMDomains;
    }
    
    public Map getFEMDomains(){
        return universe.theFEMDomains;
    }
    
    public Map<Integer,jpde> getPDEDomains(){
        return universe.thePDEDomains;
    }
    
    public Map getGENDomains(){
        return universe.theGENDomains;
    }
    
    public Map getContraptions(){
        return universe.theContraptions;
    }
    
    public void clsDEM(){
        theDEMDomains.clear();
        (new DEMdomain()).clsNumberOfDomains();
    }
    
    public void clsBEM(){
        theBEMDomains.clear();
        (new jbem.Domain()).clsNumberOfDomains();
    }
    
    public void clsFEM(){
        theFEMDomains.clear();
        (new jfem.Domain()).clsNumberOfDomains();
    }
    
    public void clsPDE(){
        thePDEDomains.clear();
        numofPDEDomains=0;
    }
    
    public void clsGEN(){
        theGENDomains.clear();
        (new gendomain.Domain()).clsNumberOfDomains();
    }
    
    public void clsContraptions(){
        theContraptions.clear();
    }
    
    public void cls(){
        clsDEM();
        clsBEM();
        clsFEM();
        clsPDE();
        clsGEN();
        clsContraptions();
    }
}
