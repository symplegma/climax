/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jbem;

import jmat.AbstractMatrix;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;

/**
 *
 * @author pchr
 */
public class MatrixStoreFile extends MatrixStore{
    private Map<Integer,preMatrixInfo> theMatrices = new HashMap<Integer,preMatrixInfo>();
    String directoryName=null;
    // constructor
    public MatrixStoreFile( ){
    }
    
    public MatrixStoreFile(String thePath){
        directoryName=thePath;
        this.checkForFiles();
    }

    @Override
    public void putMatrix(int m, int n, int pos, AbstractMatrix theMatrix) {
        ++numOfMatrices;
        preMatrixInfo info = new preMatrixInfo(m-n,pos);
        this.theMatrices.put(numOfMatrices, info);
        String fname=info.get_Name();
        if(this.directoryName==null){
            try {
                PrintWriter out = new PrintWriter(new FileWriter(fname));
                theMatrix.print(out, 25,15);
                out.close();
            } catch (IOException e) {
                System.exit(300);
            }
        }else{
            try {
                PrintWriter out = new PrintWriter(new FileWriter(new File(directoryName+"/"+fname)));
                theMatrix.print(out, 25,15);
                out.close();
            } catch (IOException e) {
                System.exit(350);
            }
        }
        
    }

    @Override
    public boolean isStored(int m, int n, int pos) {
        boolean check=false;
        int mmn=m-n;
        for(Iterator<preMatrixInfo>it=this.theMatrices.values().iterator();it.hasNext();){
            preMatrixInfo apreMatrix=it.next();
            int the_m=apreMatrix.get_m_minus_n();
            int the_pos=apreMatrix.get_pos();
            if( (the_m==mmn)&&(the_pos==pos) ){
                check=true;
            }
            //if(check!=true)break;
            if(check==true)break;
        }
        return check;
    }

    @Override
    public AbstractMatrix getMatrix(int m, int n, int pos) {
        boolean check=false;
        int mmn=m-n;
        AbstractMatrix Mat = null;
        for(Iterator<preMatrixInfo>it=this.theMatrices.values().iterator();it.hasNext();){
            preMatrixInfo apreMatrix=it.next();
            int the_m=apreMatrix.get_m_minus_n();
            int the_pos=apreMatrix.get_pos();
            if( (the_m==mmn)&&(the_pos==pos) ){
                //Mat=apreMatrix.get_Matrix();
                String fname=apreMatrix.get_Name();
                if(this.directoryName==null){
                    try {
                        BufferedReader in = new BufferedReader(new FileReader(fname));
                        Mat=AbstractMatrix.read(in);
                        in.close();
                    } catch (IOException e) {
                        System.exit(400);
                    }
                }else{
                    try {
                        BufferedReader in = new BufferedReader(new FileReader(new File(directoryName+"/"+fname)));
                        Mat=AbstractMatrix.read(in);
                        in.close();
                    } catch (IOException e) {
                        System.exit(450);
                    }
                }
                check=true;
            }
            if(check==true)break;
        }
        return Mat;
    }
    
    private void checkForFiles(){
        File directory;        // File object referring to the directory.
        String[] files=null;        // Array of file names in the directory.
        Scanner scanner;       // For reading a line of input from the user.
        
        scanner = new Scanner(System.in);  // scanner reads from standard input.
        
        //System.out.print("Enter a directory name: ");
        //directoryName = scanner.nextLine().trim();
        directory = new File(directoryName);
        
        if (directory.isDirectory() == false) {
            if (directory.exists() == false)
                System.out.println("There is no such directory: "+directoryName);
            else
                System.out.println("That file is not a directory.");
        }
        else {
            files = directory.list();
            int mmn; int pos;
            String[] temp;
            for (int i = 0; i < files.length; i++){
                if(files[i].contains(".spm")){
                    String sub=files[i].substring(3, files[i].length()-4);
                    System.out.println(sub);
                    temp = sub.split("_");
                    mmn=Integer.parseInt(temp[0]); pos=Integer.parseInt(temp[1]);
                    preMatrixInfo info = new preMatrixInfo(mmn,pos);
                    this.theMatrices.put(numOfMatrices, info);
                    ++numOfMatrices;
                }
                //System.out.println("   " + files[i]);
            }   
        }
    }

}
