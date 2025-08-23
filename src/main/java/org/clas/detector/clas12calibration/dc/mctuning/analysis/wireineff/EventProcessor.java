/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import org.clas.detector.clas12calibration.dc.mctuning.viewer.WireIneffAnalViewer;

import org.clas.detector.clas12calibration.viewer.Driver;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;
/**
 *
 * @author veronique
 */


public class EventProcessor {
    private int runNumber = -1;
    private int invocationCount = 0;
    
    // Extracts the trajectory bank from the event.
    public DataBank extractTrajectoryBank(DataEvent ev) {
        // Assuming 'getTrajectoryBank' is a method to extract the bank from the event
        DataBank bank = ev.getBank("TrajectoryBank");
        return bank;
    }
    public static double[][][]pars;
    
    // Processes the event to populate total and effective hit arrays.
    public void processEvent(DataEvent event, int[][] totLayA, int[][] effLayA, int nBins, 
            FitManager fm) {
        if (!event.hasBank("RUN::config")) {
            return ;
        }
        
        DataBank bank = event.getBank("RUN::config");
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) {
           return ;
        } else {
           invocationCount++;
           //System.out.println("PROCESSING EVENTS..........................................................."+bank.getInt("event", 0));
            if (invocationCount == 1) { System.out.println("PROCESSING EVENTS...........................................................1");
                Driver.init();
                runNumber = newRun; 
                System.out.println("Looking for PARAMETERS "+(fm!=null));
                pars = readFcnPars(runNumber);
                fm.intializeFits(pars);
            }
        }
     
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return ;
        } 
        
        DataBank ebank = SegmentAnalysis.getLayIneffBank(event, nBins, runNumber);
        if (ebank == null) return;
        
        
        int nrows =  ebank.rows();
        //Bank.show(); System.out.println(" NUMBER OF ENTRIES IN BANK = "+nrows);
        for (int i = 0; i < nrows; i++) {
            
            int bb = this.getTrkDocBin(ebank.getFloat("trkDoca", i), nBins); //normalized docas
            if(bb==-1) continue;
            //System.out.println("....................normalized doca "+Bank.getFloat("trkDoca", i)+" bin "+bb);
            if(bb<nBins) {
                //(int)((Math.floor(Math.abs(Bank.getFloat("trkDoca", i))/(2.*Constants.wpdist[Bank.getByte("superlayer", i)-1])/trkDBinning)))
                totLayA[ebank.getByte("superlayer", i)-1][bb]++;
                if(ebank.getShort("matchedHitID", i)!=-1) {
                        effLayA[ebank.getByte("superlayer", i)-1][bb]++;
                }
            }
        }
    }
    public int getRunNumber() {
        return runNumber;
    }

    // Helper function to determine the track doca bin index
    public int getTrkDocBin(double d, int nBins) {
        int bin = -1;
        double binWidth = 1./(double)nBins;
        int n = nBins;
        
        for (int i = 0; i < n; i++) {
            double blo = i*binWidth;
            double bhi = (i+1)*binWidth;
            if(d>=blo && d<bhi)
                bin = i;
        }
        return bin;
    }
    public static double[][][] readFcnPars(int run) {
        IndexedTable tab = WireIneffAnalViewer.ccdb.getConstants(run, 
                    "/calibration/dc/signal_generation/inefficiency");
        double[][][] pars = new double[6][6][4]; //6 sectors, 6 superlayers, 4 parameters    
        
        for(int s = 0; s<6; s++) {
            for(int l = 0; l<6; l++) {
                for(int p = 0; p<4; p++) {
                    String par = "p"+(p+1);
                    pars[s][l][p]=tab.getDoubleValue(par, s+1, l+1, 0);
                    System.out.println("par["+s+"]["+l+"]["+p+"] = "+pars[s][l][p]);
                }
            }
        }
        return pars;
    }
    
}
