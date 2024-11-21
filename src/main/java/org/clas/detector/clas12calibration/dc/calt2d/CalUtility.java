/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.io.File;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.AlphaBinHalfWidth;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.AlphaValues;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.AlphaValuesUpd;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BfieldValues;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BfieldValuesUpd;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.alphaBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.field;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.rec.dc.hit.FittedHit;

/**
 *
 * @author ziegler
 */
public class CalUtility {
     
    public static void reCook(Map<Coordinate, H1F> timeResi, Map<Coordinate, H1F> timeResiNew, 
            Map<Coordinate, H1F> timeResiB, Map<Coordinate, H1F> A, Map<Coordinate, H1F> B, Map<Coordinate, H1F> BAlphaBins,
            Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, H2F> Tvscalcdocas, Map<Coordinate, H2F> Tresvstrkdocas, 
            HipoDataSource calreader, List<FittedHit> hits, List<FittedHit> calhits, 
            Map<Coordinate, H1F> ParsVsIter, Map<Coordinate, FitLine> TvstrkdocasFits, Map<Coordinate, GraphErrors> TvstrkdocasProf,
            boolean useBProf) {
        for (int i = 0; i < 6; i++) {
            timeResi.get(new Coordinate(i)).reset();
            timeResiNew.get(new Coordinate(i)).reset();
            if(i==2 || i==3) {
                for (int k = 0; k < BBins; k++) {
                     timeResiB.get(new Coordinate(i,k)).reset();
                }
            }
        }
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < alphaBins; j++) {
                for (int k = 0; k < BBins+1; k++) {
                    if(A.containsKey(new Coordinate(i,j,k))) A.get(new Coordinate(i,j,k)).reset();
                    if(k<BBins && B.containsKey(new Coordinate(i,j,k))) B.get(new Coordinate(i,j,k)).reset();
                }
            }
        }
        for (int s = 0; s < 7; s++) {
            for (int i0 = 0; i0 < 6; i0++) {
                int i = i0+s*6;
                for (int j = 0; j < alphaBins; j++) {
                    for (int k = 0; k < BBins+1; k++) {
                        if((i0<2 || i0>3) && k<BBins) continue;
                        Tvstrkdocas.get(new Coordinate(i,j,k)).reset();
                        Tvscalcdocas.get(new Coordinate(i,j,k)).reset();
                        Tresvstrkdocas.get(new Coordinate(i,j,k)).reset();

                    }
                }
            }
        }
        System.out.println("***********************************************");       
        System.out.println("****** Reprocessing TestCalOutPut.hipo ********");
        System.out.println("***********************************************");
        calreader = new HipoDataSource();
        calreader.open("TestCalOutPut.hipo");
        System.out.println("Events in hipofile " +  calreader.getSize() );  
        int numberofeventsinfile = calreader.getSize();
        int eventcounter = 0;
        while (calreader.hasEvent()) { 
            hits.clear();
            DataEvent event = calreader.getNextEvent();
            eventcounter++;
            if ((eventcounter%10000 == 0) && (eventcounter < numberofeventsinfile) ) {
             	 System.out.println("Processed " + eventcounter + " events from " + numberofeventsinfile);  
            }
            if(event.hasBank("TimeBasedTrkg::TBHits")) { 
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    FittedHit h = HitUtility.getCalHit(bnkHits, i);
                    if(h!=null)
                        hits.add(h);
                }
                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    HitUtility.updateHit(hit, true); //update the hit with previous cal results
                }
                //refit with new constants
                Refit rf = new Refit();
                rf.reFit(hits, true);    //refit to get the parameters
                
                for(FittedHit hit : hits) {
                    if(hit.get_OutOfTimeFlag()==true) continue;
                    //filling the timeResi for the previously calibrated hits
                    timeResi.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                }
            } 
        }
        
        // the newly calibrated hits
        System.out.println("*************************************************");       
        System.out.println("*** Done Reprocessing with initial parameters ***");
        System.out.println("*************************************************");
        
        
        FitUtility.reLoadFitPars(ParsVsIter);
        //Parameters are now fit values
        System.out.println("************  Fit Parameters Reloaded! ************");
        calreader.gotoEvent(0);
        eventcounter = 0;
        while (calreader.hasEvent()) { 
            calhits.clear();
            DataEvent event = calreader.getNextEvent();
            eventcounter++;
            if ((eventcounter%10000 == 0) && (eventcounter < numberofeventsinfile) ) {
            	 System.out.println("Processed " + eventcounter + " events from " + numberofeventsinfile);  
            }
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    FittedHit h = HitUtility.getCalHit(bnkHits, i);
                    if(h!=null)
                        calhits.add(h);
                }
                for(FittedHit hit : calhits) {
                    hit.set_TimeResidual(-999);
                    HitUtility.updateHit(hit, true);
                }
                //refit with new constants
                Refit rf = new Refit();
                rf.reFit(calhits, true);    //use only not out of time hits to get the calib cst
                for(FittedHit hit : calhits) {
                    //if(hit.get_OutOfTimeFlag()) continue;
                    double theta0 = Math.toDegrees(Math.acos(1-0.02*hit.getB()));
                    double alphaUncor = hit.getAlpha()+(double)T2DCalib.polarity*theta0;
                    int alphaBin = getAlphaBin(alphaUncor);
                    double bFieldVal = (double) hit.getB();
                    int slyrIdx = (hit.get_Superlayer()-1)+(hit.get_Sector()-1)*6;
                    int allSlyrIdx = 35+hit.get_Superlayer();
                    if(alphaBin!=-1 && !hit.get_OutOfTimeFlag()) {
                       
                        double calibTime = (double) (hit.get_TDC() - hit.getTProp()
                                            - hit.getTFlight() - hit.getTStart()
                                            - hit.getT0()); 
                        
                        double yf = TvstrkdocasFits.get(new Coordinate(hit.get_Superlayer()+(hit.get_Sector()-1)*6 - 1,
                                alphaBin, BBins)).evaluate(hit.get_ClusFitDoca());
                        
                            Tvstrkdocas.get(new Coordinate(slyrIdx, alphaBin, BBins))
                                        .fill(hit.get_ClusFitDoca(), calibTime);
                            Tvstrkdocas.get(new Coordinate(allSlyrIdx, alphaBin, BBins))
                                        .fill(hit.get_ClusFitDoca(), calibTime);
                            Tvscalcdocas.get(new Coordinate(slyrIdx, alphaBin, BBins))
                                        .fill(hit.get_Doca(), calibTime);
                            Tvscalcdocas.get(new Coordinate(allSlyrIdx, alphaBin, BBins))
                                        .fill(hit.get_Doca(), calibTime);
                        if(!Double.isNaN(yf)) {
                            Tresvstrkdocas.get(new Coordinate(slyrIdx, alphaBin, BBins)).fill(hit.get_ClusFitDoca(), calibTime-yf);
                            Tresvstrkdocas.get(new Coordinate(allSlyrIdx, alphaBin, BBins)).fill(hit.get_ClusFitDoca(), calibTime-yf);
                        }
                        //Fill region 2 for different b-field values
                        if(hit.get_Superlayer() >2 && hit.get_Superlayer() <5) { 
                            int bBin = getBBin(bFieldVal);
                            double r2yf = TvstrkdocasFits.get(new Coordinate(slyrIdx, alphaBin, bBin)).evaluate(hit.get_ClusFitDoca());
                            
                                Tvstrkdocas.get(new Coordinate(slyrIdx, alphaBin, bBin))
                                            .fill(hit.get_ClusFitDoca(), calibTime);
                                Tvstrkdocas.get(new Coordinate(allSlyrIdx, alphaBin, bBin))
                                            .fill(hit.get_ClusFitDoca(), calibTime);
                                Tvscalcdocas.get(new Coordinate(slyrIdx, alphaBin, bBin))
                                            .fill(hit.get_Doca(),calibTime);
                                Tvscalcdocas.get(new Coordinate(allSlyrIdx, alphaBin, bBin))
                                            .fill(hit.get_Doca(),calibTime);
                            if(!Double.isNaN(r2yf)) {
                                Tresvstrkdocas.get(new Coordinate(slyrIdx, alphaBin, bBin))
                                        .fill(hit.get_ClusFitDoca(), calibTime-r2yf);
                                Tresvstrkdocas.get(new Coordinate(allSlyrIdx, alphaBin, bBin))
                                        .fill(hit.get_ClusFitDoca(), calibTime-r2yf);
                            }
                        }
                        
                        if(hit.get_Superlayer()<3 || hit.get_Superlayer()>4) { 
                            A.get(new Coordinate(hit.get_Superlayer()-1, alphaBin, BBins))
                                        .fill(alphaUncor);
                            //System.out.println("Filling alpha for spl "+superlayer+" alpha bin "+alphaBin+" alpha "+alphaUncor);
                        }

                        // fill B values histograms
                        if(hit.get_Superlayer() ==3 || hit.get_Superlayer() ==4) {
                            int bBin = getBBin(bFieldVal);
                            B.get(new Coordinate(hit.get_Superlayer()-3, alphaBin, bBin))
                                    .fill(bFieldVal);
                            A.get(new Coordinate(hit.get_Superlayer()-1, alphaBin, getBBin(bFieldVal)))
                                    .fill(alphaUncor);
                            //System.out.println("Filling alpha for spl "+superlayer+" alpha bin "+alphaBin+" Bin "+this.getBBin(bFieldVal)+" alpha "+alphaUncor);
                        }
                        
                    }
                }
                
                rf.reFit(calhits, false);    //use all hits
                for(FittedHit hit : calhits) {
                    //filling the timeResi for the newly calibrated hits
                    timeResiNew.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                    if(hit.get_Superlayer() ==3 || hit.get_Superlayer() ==4) {
                        double bFieldVal = (double) hit.getB();
                        int bBin = getBBin(bFieldVal);
                        timeResiB.get(new Coordinate(hit.get_Superlayer()-1, bBin)).fill(hit.get_TimeResidual());
                    }
                }
            }
        }
        UpdateAlphaBinCenters(A);
        UpdateBBinCenters(B, BAlphaBins);
        System.out.println("REMAKING PROFILES");
        for(int s =0; s<7; s++) {
            for (int i0 = 0; i0 < 6; i0++) {
                int i = i0+s*6;
                for (int j = 0; j < T2DCalib.alphaBins; j++) {
                    GraphUtility.filltrkDocavsTGraphs(i,j, Tvstrkdocas, TvstrkdocasProf, useBProf);
                    System.out.println("PROFILE "+i+" "+j+" is OK!");
                }
            }   
        }
        System.out.println("RECOOKING DONE WITH THE NEW CONSTANTS!");
        System.out.println("CHECK THE RESIDUALS!");
        calreader.close();    
    }
    
    public static int getAlphaBin(double alpha) {
        int v = -1;
        for(int i = 0; i<T2DCalib.AlphaValues.length; i++) {
            
            if(Math.abs(alpha-T2DCalib.AlphaValues[i])<T2DCalib.AlphaBinHalfWidth)
                v = i;
        } 
        
        return v;
    }

    public static int getBBin(double bFieldVal) {
        
        int v = BfieldValues.length-1;
        //BfieldValues = new double[]{0.0000, 1.0000, 1.4142, 1.7321, 2.0000, 2.2361, 2.4495, 2.6458};
        //BfieldValues^2 = new double[]{0.0000, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000};
        double BSqrBinHalfWidth = 0.5;
        for(int i = 0; i<BfieldValues.length; i++) {
            if(Math.abs(bFieldVal*bFieldVal-T2DCalib.BfieldValues[i]*T2DCalib.BfieldValues[i])<BSqrBinHalfWidth)
                v = i;
        }      
        
        //return bbinIdx ;
        return v ;
    }
    
    public static void UpdateAlphaBinCenters(Map<Coordinate, H1F> A) {
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < alphaBins; j++) {
                for (int k = 0; k <= BBins; k++) {
                    AlphaValuesUpd[i][j][k] = AlphaValues[j];
                    if(A.get(new Coordinate(i,j,k)).getBinContent(A.get(new Coordinate(i,j,k)).getMaximumBin())>0) {
                        AlphaValuesUpd[i][j][k] = A.get(new Coordinate(i,j,k)).getMean();
                    System.out.println("ijk" +i+" "+j+" "+k+" Alpha Bin UPdated "+AlphaValues[j]+" --> "+AlphaValuesUpd[i][j][k] );
                    }
                }
            }
        }
    }
    static boolean filledBspectra = false;
    public static void UpdateBBinCenters(Map<Coordinate, H1F> B, Map<Coordinate, H1F> BAlphaBins) {
        if(field==0) return;
        if(filledBspectra) return;
        TCanvas can1 = new TCanvas("superlayer3 B", 800, 800);
        TCanvas can2 = new TCanvas("superlayer4 B", 800, 800);
        
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < alphaBins; j++) {
                for (int k = 0; k < BBins; k++) {
                    BfieldValuesUpd[i][j][k] = BfieldValues[k];
                    if(B.get(new Coordinate(i,j,k)).getBinContent(B.get(new Coordinate(i,j,k)).getMaximumBin())>10) {
                        BfieldValuesUpd[i][j][k] = B.get(new Coordinate(i,j,k)).getMean();
                    System.out.println("ijk" +i+" "+j+" "+k+" BBin UPdated "+BfieldValues[k]+" --> "+BfieldValuesUpd[i][j][k] );
                    }
                }
            }
        }
        can1.divide(7, 3);
        can2.divide(7, 3);
        int can1idx=0;
        int can2idx=0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < alphaBins; j++) {
                BAlphaBins.get(new Coordinate(i,j)).setOptStat(0);
                BAlphaBins.get(new Coordinate(i,j)).setTitle("alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")");
                
                for (int k = 0; k < BBins; k++) {
                    if(B.get(new Coordinate(i,j,k)).getBinContent(B.get(new Coordinate(i,j,k)).getMaximumBin())>1)
                        BAlphaBins.get(new Coordinate(i,j)).add(B.get(new Coordinate(i,j,k)));
                }    
                if(i==0) {
                        can1.cd(can1idx);
                        can1.draw(BAlphaBins.get(new Coordinate(i,j)));
                        can1idx++;
                }
                if(i==1) {
                    can2.cd(can2idx);
                    can2.draw(BAlphaBins.get(new Coordinate(i,j)));
                    can2idx++;
                }
            }
        }
        filledBspectra=true;
    }

    public static void checkFile(String testCalOutPuthipo) {
        File file = new File(testCalOutPuthipo);
        if (file.exists()) {
            // Delete the file
            if (file.delete()) {
                System.out.println("File deleted successfully.");
            } else {
                System.out.println("Failed to delete the file.");
            }
        } else {
            System.out.println("File does not exist.");
        }
    }
    
    
}
