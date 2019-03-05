/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mathman;

/**
 *
 * @author pchr
 */
public class LocalPowerDoubleFunction  implements DoubleFunction{
    double x_ast=0.0;
    double y_ast=0.0;
    double tol=0.00001;
    double power=3;
    double amp=1.0;
    
    public LocalPowerDoubleFunction(){}
    
    public LocalPowerDoubleFunction(double x_ast){
        this.x_ast=x_ast;
    }
    
    public LocalPowerDoubleFunction(double x_ast, double y_ast){
        this.x_ast=x_ast;
        this.y_ast=y_ast;
    }
    
    public void setTol(double val){this.tol=val;}
    
    public void setPower(double val){this.power=val;}
    
    public void setAmplification(double val){this.amp=val;}
    
    @Override
    public double run(double x) {
        double val=0.0;
        double dist=Math.sqrt((x-x_ast)*(x-x_ast));
        if(dist<=tol){val=1.0-(dist*dist)/(tol*tol);}
        return amp*Math.pow(val, power);
    }

    @Override
    public double run(double x, double y) {
         double val=0.0;
        double dist=Math.sqrt((x-x_ast)*(x-x_ast)+(y-y_ast)*(y-y_ast));
        if(dist<=tol){val=1.0-(dist*dist)/(tol*tol);}
        return amp*Math.pow(val, power);
    }
}
