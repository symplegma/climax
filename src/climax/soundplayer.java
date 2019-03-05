/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package climax;

import java.util.ArrayList;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.SourceDataLine;
import javax.sound.sampled.LineUnavailableException;
/**
 * http://digitalsoundandmusic.com/2-3-13-modeling-sound-in-java/
 * @author pchr
 */
public class soundplayer {
    // constructor
    public soundplayer(){}
    
    public static void playsound(double[] data, int volume) throws LineUnavailableException {
        byte[] buf;
        AudioFormat audioF;
        
        buf =  new byte[1];
        audioF = new AudioFormat(data.length,8,1,true,false);
        
        //sampleRate, sampleSizeInBits,channels,signed,bigEndian
        SourceDataLine sourceDL = AudioSystem.getSourceDataLine(audioF);
        sourceDL = AudioSystem.getSourceDataLine(audioF);
        sourceDL.open(audioF);
        sourceDL.start();
        
        for(int i=0; i<data.length; i++){
            buf[0]=(byte)(data[i]*volume);
            sourceDL.write(buf,0,1);
        }
        
        sourceDL.drain();
        sourceDL.stop();
        sourceDL.close();
    }
    
    public static void playsound(ArrayList data, int volume) throws LineUnavailableException {
        double[] x=new double[data.size()]; 
        for(int i=0;i<x.length;i++)x[i]=(double)data.get(i);
        playsound(x, volume);
    }
    
    public static void playsound(double[] data) throws LineUnavailableException {
        playsound(data, 30);
    }
    
    public static void playsound(ArrayList data) throws LineUnavailableException {
        double[] x=new double[data.size()]; 
        for(int i=0;i<x.length;i++)x[i]=(double)data.get(i);
        playsound(x, 30);
    }
}
