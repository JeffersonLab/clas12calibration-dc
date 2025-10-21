/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import cnuphys.magfield.MagneticFieldInitializationException;
import cnuphys.magfield.MagneticFields;
import java.io.FileNotFoundException;
import java.util.List;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.WireIneffAnalViewer;
import static org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.WireIneffAnal.Bsq;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.system.ClasUtilsFile;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
/**
 *
 * @author veronique
 */


public class EventProcessor {
    public static int runNumber = -1;
    private int invocationCount = 0;
    public static Swim   swim     = null;
       
    String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        
    private void initField(Double torusScale, Double solenoidScale){
        
        String magfieldDir = ClasUtilsFile.getResourceDir("CLAS12DIR", "../etc/data/magfield/");

        String torusFileName = System.getenv("COAT_MAGFIELD_TORUSMAP");
        if (torusFileName==null) torusFileName = "Symm_torus_r2501_phi16_z251_24Apr2018.dat";

        String solenoidFileName = System.getenv("COAT_MAGFIELD_SOLENOIDMAP");
        if (solenoidFileName==null) solenoidFileName = "Symm_solenoid_r601_phi1_z1201_13June2018.dat";
        try {
            MagneticFields.getInstance().initializeMagneticFields(magfieldDir,torusFileName,solenoidFileName);
        }
        catch (MagneticFieldInitializationException | FileNotFoundException e) {
            System.err.println("MAGFIELD ERROR..............");
            e.printStackTrace();
            return ;
        }
        
		
        MagneticFields.getInstance().getTorus().setScaleFactor(torusScale);
        MagneticFields.getInstance().getSolenoid().setScaleFactor(solenoidScale);
        
        Swimmer.setMagneticFieldsScales(solenoidScale, torusScale, 0);
        swim = new Swim();
    }    
    // Extracts the trajectory bank from the event.
    public DataBank extractTrajectoryBank(DataEvent ev) {
        // Assuming 'getTrajectoryBank' is a method to extract the bank from the event
        DataBank bank = ev.getBank("TrajectoryBank");
        return bank;
    }
    public static double[][][]pars;
    // Processes the event to populate total and effective hit arrays.
    
    public void processEvent(DataEvent event, int[][][] totLayA, int[][][] effLayA, int nBins, 
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
                fm.intializeFits();
                System.out.println("INITIALIZING THE FIELD...."+bank.getFloat("torus", 0)+" "+bank.getFloat("solenoid", 0));
                initField((double)bank.getFloat("torus", 0),(double)bank.getFloat("solenoid", 0));
            }
        }
        
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return ;
        } 
        
        List<SegmentTrajectoryRecord> trjrecs = SegmentAnalysis.getLayIneffBank(event, nBins, runNumber);
        if (trjrecs == null) return;
        
        
        int nrows =  trjrecs.size();
        //Bank.show(); System.out.println(" NUMBER OF ENTRIES IN BANK = "+nrows);
        for (int i = 0; i < nrows; i++) {
            
            int tb = this.getTrkDocBin(trjrecs.get(i).getTrkDoca(), nBins); //normalized docas
            int bb = getBBin(trjrecs.get(i).getTrkBfield());
            if(trjrecs.get(i).getSuperlayer()==3 || trjrecs.get(i).getSuperlayer()==4)
                HistogramManager.Bhists.get(new Coordinate(trjrecs.get(i).getSuperlayer()-3,tb, bb)).fill(trjrecs.get(i).getTrkBfield()*trjrecs.get(i).getTrkBfield());
            if(tb==-1) continue;
            if(tb<nBins) {
                totLayA[trjrecs.get(i).getSuperlayer()-1][tb][bb]++;
                if(trjrecs.get(i).getMatchedHitID()!=-1) {
                        effLayA[trjrecs.get(i).getSuperlayer()-1][tb][bb]++;
                }
            }
        }
    }
    public int getRunNumber() {
        return runNumber;
    }
    // Helper function to determine the Bfield bin index
    private static int getBBin(double bFieldVal) {
        int v = Bsq.length-1;
        double BSqrBinHalfWidth = Bsq[0];
        for(int i = 0; i<Bsq.length; i++) {
            if(Math.abs(bFieldVal*bFieldVal-Bsq[i])<BSqrBinHalfWidth) {
                v = i;
            }
        }      
        return v ;
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
    
    public static boolean getWireStatus(int run, int sector, int superlayer, int layer, int wire) {
        
        IndexedTable tab = WireIneffAnalViewer.ccdb.getConstants(run, 
                    "/calibration/dc/tracking/wire_status");
        int slayer = layer+(superlayer-1)*6;
        int status = tab.getIntValue("status", sector, slayer, wire);
        
        boolean pass = true;
        if(status!=0) pass = false;
        
        return pass;
    }
    
}
