/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdem;

import gendomain.Node;

/**
 *
 * @author pchr
 */
public class InitialConditionVelocity extends InitialCondition{
    
    public InitialConditionVelocity(double InitialValue, Node wnode, int wdof) {
        super(InitialValue, wnode, wdof);
        this.type=1;
    }
    
}
