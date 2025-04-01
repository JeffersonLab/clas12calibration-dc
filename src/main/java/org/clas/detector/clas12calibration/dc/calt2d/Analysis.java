/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.ArrayList;
import java.util.List;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.hitBank;
import org.jlab.groot.data.H1F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

/**
 *
 * @author ziegler
 */
public class Analysis {
    
    public static void main(String[] args) {
        List<String> inputList = new ArrayList<>();
        inputList.add("TestCalOutPut.hipo");
        H1F h = new H1F("h", "","", 200, 0, 1.2);
        H1F h2 = new H1F("h2", "","", 200, 0, 1.2);
        h.setLineColor(2);
        int counter=0;
        for(String inputFile :  inputList) {
            HipoDataSource reader = new HipoDataSource();
            reader.open(inputFile);
            while (reader.hasEvent()) { 
                counter++;
                DataEvent event = reader.getNextEvent();
                if(event.hasBank(hitBank)) { 
                    DataBank bnkHits = event.getBank(hitBank);

                    for (int i = 0; i < bnkHits.rows(); i++) {
                        int sly = bnkHits.getByte("superlayer", i);
                        double tdoca = bnkHits.getFloat("trkDoca", i);
                        double time = bnkHits.getFloat("trkDoca", i);
                        double B = bnkHits.getFloat("B", i);
                        double resiTime = bnkHits.getFloat("timeResidual", i);
                        if(sly==3 || sly==4) {
                        if(resiTime<-0.15 && B>1)
                            h.fill(tdoca);
                        if(Math.abs(resiTime)<0.08 && B>1)
                            h2.fill(tdoca);
                        }
                    }
                    if(counter%1000 ==1) System.out.println(counter);
                }
            }
        }
        TCanvas can = new TCanvas("can", 1200, 800);
        can.divide(1, 2);
        can.cd(0);
        can.draw(h, "");
        can.cd(1);
        can.draw(h2, "");
    }
}
