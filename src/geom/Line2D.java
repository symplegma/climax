/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geom;

import java.util.Iterator;
import jbem.Element;
import jbem.Node;

/**
 *
 * @author pchr
 */
public class Line2D extends geom.Curve{
    
    public Line2D(){}
    
    public Line2D(int id, geom.Point p1, geom.Point p2){
        this.id=id;
        putPoint(p1);
        putPoint(p2);
    }
    
    public Line2D(geom.Point p1, geom.Point p2){
        putPoint(p1);
        putPoint(p2);
    }
    
    public geom.Point getPstart(){
        geom.Point p1 = null;
        geom.Point p2 = null;
        int count=0;
        for(Iterator<geom.Point> it=this.getPoints().values().iterator(); it.hasNext();){
            if(count==0){p1 = it.next();}else{p2= it.next();}
            count++;
        }
        return p1;
    }
    
    public geom.Point getPend(){
        geom.Point p1 = null;
        geom.Point p2 = null;
        int count=0;
        for(Iterator<geom.Point> it=this.getPoints().values().iterator(); it.hasNext();){
            if(count==0){p1 = it.next();}else{p2= it.next();}
            count++;
        }
        return p2;
    }
    
    public boolean contains(Element elem){
        boolean lieson=true;
        for(Iterator<Node> nit=elem.getNodes().values().iterator();nit.hasNext();){
            Node elNode=nit.next();
            if(!this.contains(elNode))lieson=false;
        }
        return lieson;
    }
}
