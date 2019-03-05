/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;
/**
 *
 * @author pchr
 */
public class main {
    
        public static void main(String[] args) {
            //===================================================
            // Parsing
            //===================================================
            Parser theParser = null;
            if(args.length>=1){
                //theParser = new FileParser();
            }else{
                theParser = new InParser();
            }
            if(theParser!=null) {
                theParser.Parse();
            }
        }
}
