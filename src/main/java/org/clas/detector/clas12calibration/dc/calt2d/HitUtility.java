/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.avebeta;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.betacnt;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.polarity;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.timetodistance.TableLoader;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;

/**
 *
 * @author ziegler
 */
public class HitUtility {

    static TimeToDistanceEstimator tdee = new TimeToDistanceEstimator();
    
    public static FittedHit getCalHit(DataBank bnkHits, int i) {
        
        int id = bnkHits.getShort("id", i);
        int sector = bnkHits.getByte("sector", i);
        int superlayer = bnkHits.getByte("superlayer", i);
        int layer = bnkHits.getByte("layer", i);
        int wire = bnkHits.getShort("wire", i);
        int TDC = bnkHits.getInt("TDC", i);
        double doca = bnkHits.getFloat("doca", i);
        double docaError = bnkHits.getFloat("docaError", i);
        double trkDoca = bnkHits.getFloat("trkDoca", i);
        int LR = bnkHits.getByte("LR", i);
        double X = bnkHits.getFloat("X", i);
        double Z = bnkHits.getFloat("Z", i);
        double B = bnkHits.getFloat("B", i);
        double TProp = bnkHits.getFloat("TProp", i);
        double TFlight = bnkHits.getFloat("TFlight", i);
        double T0 = bnkHits.getFloat("T0", i);
        double TStart = bnkHits.getFloat("TStart", i);
        int clusterID = bnkHits.getShort("clusterID", i);
        int trkID = bnkHits.getByte("trkID", i);
        double time = bnkHits.getFloat("time", i);
        double beta = bnkHits.getFloat("beta", i);
        double tBeta = bnkHits.getFloat("tBeta", i);
        double resiTime = bnkHits.getFloat("timeResidual", i);
        double resiFit = bnkHits.getFloat("fitResidual", i);
        int jitter = (int) bnkHits.getByte("jitter", i);
        double dDoca = bnkHits.getFloat("dDoca", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        //double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        //double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double alpha = bnkHits.getFloat("Alpha", i);

        FittedHit hit  = new FittedHit(sector, superlayer, layer, wire, TDC, jitter, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
        hit.setB(B);
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(T2DViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(T2DViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
        hit.setAlpha(alpha);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        hit.set_DeltaDocaBeta(dDoca);
        //this.getSegProperty(bnkHits);
        //if(this.passCalibCuts(bnkHits, i)) {            
            return hit;
        //} else {
        //    return null;
        //}
    }
    public static FittedHit getHit(DataBank bnkHits, int i) {
        
        int id = bnkHits.getShort("id", i);
        int sector = bnkHits.getByte("sector", i);
        int superlayer = bnkHits.getByte("superlayer", i);
        int layer = bnkHits.getByte("layer", i);
        int wire = bnkHits.getShort("wire", i);
        int TDC = bnkHits.getInt("TDC", i);
        double doca = bnkHits.getFloat("doca", i);
        double docaError = bnkHits.getFloat("docaError", i);
        double trkDoca = bnkHits.getFloat("trkDoca", i);
        int LR = bnkHits.getByte("LR", i);
        double X = bnkHits.getFloat("X", i);
        double Z = bnkHits.getFloat("Z", i);
        double B = bnkHits.getFloat("B", i);
        double TProp = bnkHits.getFloat("TProp", i);
        double TFlight = bnkHits.getFloat("TFlight", i);
        double T0 = bnkHits.getFloat("T0", i);
        double TStart = bnkHits.getFloat("TStart", i);
        int clusterID = bnkHits.getShort("clusterID", i);
        int trkID = bnkHits.getByte("trkID", i);
        double time = bnkHits.getFloat("time", i);
        double beta = bnkHits.getFloat("beta", i);
        double tBeta = bnkHits.getFloat("tBeta", i);
        double resiTime = bnkHits.getFloat("timeResidual", i);
        double resiFit = bnkHits.getFloat("fitResidual", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        //double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        //double alphaUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double alpha = bnkHits.getFloat("Alpha", i);
        int jitter = (int) bnkHits.getByte("jitter", i);
        double dDoca = bnkHits.getFloat("dDoca", i);
        FittedHit hit = new FittedHit(sector, superlayer, layer, wire, TDC, jitter, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
        hit.setB(B);
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(T2DViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(T2DViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
        hit.setAlpha(alpha);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        hit.set_DeltaDocaBeta(dDoca);
        //this.getSegProperty(bnkHits);
        //double dmax = 2.*Constants.getInstance().wpdist[superlayer-1]; 
        //double cos30minusalpha = Math.cos(Math.toRadians(30.-util.getReducedAngle(alpha)));
        return hit;
    }
    public static DataBank fillTBHitsBank(DataEvent event, List<FittedHit> hitlist) {
        //if(event.hasBank("TimeBasedTrkg::TBHits")) { // for second pass tracking
        //     event.removeBank("TimeBasedTrkg::TBHits");
        //}
        DataBank bank = event.createBank("TimeBasedTrkg::TBHits", hitlist.size());

        for (int i = 0; i < hitlist.size(); i++) {
            bank.setShort("id", i, (short) hitlist.get(i).get_Id());
            bank.setShort("status", i, (short) hitlist.get(i).get_QualityFac());
            bank.setByte("superlayer", i, (byte) hitlist.get(i).get_Superlayer());
            bank.setByte("layer", i, (byte) hitlist.get(i).get_Layer());
            bank.setByte("sector", i, (byte) hitlist.get(i).get_Sector());
            bank.setShort("wire", i, (short) hitlist.get(i).get_Wire());

            bank.setFloat("X", i, (float) hitlist.get(i).get_X());
            bank.setFloat("Z", i, (float) hitlist.get(i).get_Z());
            bank.setByte("LR", i, (byte) hitlist.get(i).get_LeftRightAmb());
            bank.setFloat("time", i, (float) (hitlist.get(i).get_Time())); //time is the fully corrected time
            bank.setFloat("tBeta", i, (float) hitlist.get(i).get_DeltaTimeBeta());
            bank.setFloat("fitResidual", i, (float) hitlist.get(i).get_Residual());
            bank.setFloat("Alpha", i, (float) hitlist.get(i).getAlpha());
            
            bank.setFloat("doca", i, (float) hitlist.get(i).get_Doca());
            bank.setFloat("docaError", i, (float) hitlist.get(i).get_DocaErr());
            bank.setFloat("trkDoca", i, (float) hitlist.get(i).get_ClusFitDoca());

            bank.setShort("clusterID", i, (short) hitlist.get(i).get_AssociatedClusterID());
            bank.setByte("trkID", i, (byte) hitlist.get(i).get_AssociatedHBTrackID());
            bank.setFloat("timeResidual", i, (float) hitlist.get(i).get_TimeResidual());
            
            bank.setInt("TDC",i,hitlist.get(i).get_TDC());
            bank.setFloat("B", i, (float) hitlist.get(i).getB());
            bank.setFloat("TProp", i, (float) hitlist.get(i).getTProp());
            bank.setFloat("TFlight", i, (float) hitlist.get(i).getTFlight());
            bank.setFloat("T0", i, (float) hitlist.get(i).getT0());
            bank.setFloat("TStart", i, (float) hitlist.get(i).getTStart());
            bank.setFloat("dDoca", i, (float) hitlist.get(i).get_DeltaDocaBeta());
            bank.setFloat("beta", i, (float) hitlist.get(i).get_Beta());
            
        }
        //System.out.println(" Created Bank "); bank.show();
        return bank;

    }
    
    static final double chi2OvNDF=5;
    private static int readPID(DataEvent event, int trkId) {
        int pid = 0;
        //fetch the track associated pid from the REC tracking bank
        if (!event.hasBank("REC::Particle") || !event.hasBank("REC::Track"))
            return pid;
        DataBank bank = event.getBank("REC::Track");
        //match the index and return the pid
        int rows = bank.rows();
        for (int i = 0; i < rows; i++) {
            if (bank.getByte("detector", i) == 6 &&
                    bank.getShort("index", i) == trkId - 1) {
                double chi2 = bank.getFloat("chi2", i);
                int NDF = bank.getShort("NDF", i);
                if(chi2/(double)NDF>chi2OvNDF) return pid;
                DataBank bank2 = event.getBank("REC::Particle");
                if(bank2.getByte("charge", bank.getShort("pindex", i))!=0) {
                    pid = bank2.getInt("pid", bank.getShort("pindex", i));
                    double chi2pid = bank2.getFloat("chi2pid", bank.getShort("pindex", i));
                    if(Math.abs(chi2pid)>10) pid =0;
                }
            }
        }
        
        return pid;
    } 
    
    private static boolean isElectronZeroField(DataEvent event, int trackId) {
       
        if(!event.hasBank("REC::Particle") ||
           !event.hasBank("REC::Calorimeter") ||
           !event.hasBank("REC::Cherenkov") ||
           !event.hasBank("REC::Track"))
            return false;
        
        DataBank particleBank    = event.getBank("REC::Particle");
        DataBank calorimeterBank = event.getBank("REC::Calorimeter");
        DataBank cherenkovBank   = event.getBank("REC::Cherenkov");
        DataBank trackBank       = event.getBank("REC::Track");
         
        int pindex = -1;
        for (int i=0; i<trackBank.rows(); i++) {
            if(trackBank.getShort("index", i)==trackId-1) {
                pindex = trackBank.getShort("pindex", i);
                break;
            }
        }    
        
        if(pindex>=0) {
            double beta    = particleBank.getFloat("beta", pindex);
            //double vtx     = particleBank.getFloat("vz", pindex);

            double nphe = 0;
            for(int i=0; i<cherenkovBank.rows(); i++) {
                if(cherenkovBank.getShort("pindex", i)==pindex) {
                    nphe = cherenkovBank.getFloat("nphe", i);
                    break;
                }
            }

            double energy = 0;
            for (int i=0; i<calorimeterBank.rows(); i++) {
                if(calorimeterBank.getShort("pindex", i)==pindex) {
                    energy+=calorimeterBank.getFloat("energy", i);
                }
            }

            if(beta>0 && nphe>2 && energy>0.5) { 
                return true;
            }
        }    
        return false;
    }    

    public static DataBank FixDoubleHitsTFlight(DataEvent event, DataBank HitsBank, T2DCalib t2DCalib, 
            Map<Integer, ArrayList<FittedHit>> segMap) {
        List<FittedHit> bnkHits = new ArrayList<>();
        for (int i = 0; i < HitsBank.rows(); i++) {
            bnkHits.add(getHit(HitsBank, i));
        }
        segMap.clear();
        for (int j = 0; j < bnkHits.size(); j++) {
            Integer cID = bnkHits.get(j).get_AssociatedClusterID();
            if (cID < 0) {
                continue;
            }
            if (segMap.containsKey(cID) == false) {
                segMap.put(cID, new ArrayList<>());
            }
            segMap.get(cID).add(bnkHits.get(j));
        }
        Iterator<Map.Entry<Integer, ArrayList<FittedHit>>> itr = segMap.entrySet().iterator();
        bnkHits.clear();
        while (itr.hasNext()) {
            Map.Entry<Integer, ArrayList<FittedHit>> entry = itr.next();
            Collections.sort(entry.getValue());
            for (int j = 0; j < entry.getValue().size() - 1; j++) {
                if (entry.getValue().get(j).get_Layer() == entry.getValue().get(j + 1).get_Layer()) {
                    if (entry.getValue().get(j).getTFlight() == 0) {
                        entry.getValue().get(j).setTFlight(entry.getValue().get(j + 1).getTFlight());
                    }
                    if (entry.getValue().get(j + 1).getTFlight() == 0) {
                        entry.getValue().get(j + 1).setTFlight(entry.getValue().get(j).getTFlight());
                    }
                }
                bnkHits.add(entry.getValue().get(j));
            }
            bnkHits.add(entry.getValue().get(entry.getValue().size() - 1));
        }
        DataBank newHitsBank = fillTBHitsBank(event, bnkHits);
        return newHitsBank;
    }

    private static double timeToDistance(int sector, int superlayer, double alpha, double beta, double B, double time) {
        return calcTimeToDistance( sector,  superlayer,  alpha,  beta,  B,  time, tdee) ;
    }
    static double calcTimeToDistance(int sector, int superlayer, double alpha, double beta, double B, double time, TimeToDistanceEstimator tdee) {
        double distance = 0;
        //reduce the corrected angle
        double ralpha = (double) FcnUtility.getReducedAngle(alpha);
        if (beta > 1.0) {
            beta = 1.0;
        }
        double correctedTime = time;
        if (correctedTime <= 0) {
            return 0; // fixes edge effects ... to be improved
        }
        distance = tdee.interpolateOnGrid(B, ralpha, beta, correctedTime, sector - 1, superlayer - 1);
        return distance;
    }

    public static void updateHit(FittedHit hit, boolean flagOT) {
        double distbeta = T2DCalib.TvstrkdocasFitPars.get(new Coordinate(hit.get_Sector()-1,hit.get_Superlayer()-1)).value(4);
        double v0 = T2DCalib.TvstrkdocasFitPars.get(new Coordinate(hit.get_Sector()-1,hit.get_Superlayer()-1)).value(0);
        double d = hit.get_ClusFitDoca();
        double beta = hit.get_Beta();
        if (beta > 1.0) {
            beta = 1.0;
        }
        //beta = T2DCalib.getBetaAve();
        double calibTime = (double) (hit.get_TDC() - hit.getTProp() - hit.getTFlight() - hit.getTStart() - hit.getT0());
        double deltatime_beta = FcnUtility.getDeltaTimeBeta(d, beta, distbeta, v0);
        hit.set_DeltaTimeBeta(deltatime_beta);
        
        hit.set_Doca(timeToDistance(hit.get_Sector(), hit.get_Superlayer(), hit.getAlpha(), hit.get_Beta(), hit.getB(), calibTime));
        double ralpha = (double) FcnUtility.getReducedAngle(hit.getAlpha());
        if (flagOT) {
            double calctime = TableLoader.calc_Time(hit.get_Doca(), ralpha, hit.getB(), hit.get_Sector(), hit.get_Superlayer());
            double deltatimebeta = FcnUtility.getDeltaTimeBeta(hit.get_Doca(), beta, distbeta, v0);
            calctime += deltatimebeta;
            if (calibTime - calctime > T2DCalib.DeltaTimeCut) {
                //System.out.println(hit.printInfo()+" alpha "+hit.getAlpha()+" ralpha "+ralpha+ " tDOCA "+d+ " DOCA "+hit.get_Doca()+ " calibTime "+calibTime+" calctime "+calctime);
                hit.set_OutOfTimeFlag(true);
            }
        }
        double x = hit.get_XWire();
        double theta0 = Math.toDegrees(Math.acos(1 - 0.02 * hit.getB()));
        //fix alpha to get the local angle
        double alphaRadUncor = Math.toRadians(hit.getAlpha() + (double) T2DCalib.polarity * theta0);
        double trkAngle = Math.tan(alphaRadUncor);
        double cosTkAng = 1.0 / Math.sqrt(trkAngle * trkAngle + 1.0);
        hit.set_X(x + hit.get_LeftRightAmb() * (hit.get_Doca() / cosTkAng));
    }
    
    public static Set<Integer> getSegNumLayer(List<FittedHit> bnkHits, Map<Integer, ArrayList<CalHit>> segMap2TBHits) {
        Set<Integer> hitBankRows = new HashSet<>();
        segMap2TBHits.clear();
               
        for (int j = 0; j < bnkHits.size(); j++){
            Integer cID = bnkHits.get(j).get_AssociatedClusterID();
            if(cID<0) continue;
            if(segMap2TBHits.containsKey(cID)==false) {
                segMap2TBHits.put(cID, new ArrayList<>()); 
            }
            CalHit ch = new CalHit();
            int wire = bnkHits.get(j).get_Wire();
            ch.idx=bnkHits.get(j).get_Id();
            ch.wire = wire;
            ch.layer = bnkHits.get(j).get_Layer();
            segMap2TBHits.get(cID).add(ch);
            
        }
        
        Iterator<Map.Entry<Integer, ArrayList<CalHit>>> itr = segMap2TBHits.entrySet().iterator(); 
          
        while(itr.hasNext()) { 
            Map.Entry<Integer, ArrayList<CalHit>> entry = itr.next(); 
            List<Integer> lys = new ArrayList<>();
            for(CalHit ch : entry.getValue()) {
                lys.add(ch.layer);
            }
            Collections.sort(lys);
            int size=0;
            int l =0;
            for(int i = 0; i<lys.size(); i++) {
                if(lys.get(i)!=l)
                    size++;
                l = lys.get(i); 
            }
            if(size>=Integer.parseInt(T2DViewer.numLayers.getText())) {
                for(CalHit ch : entry.getValue()) {
                    hitBankRows.add(ch.idx);
                }
            } 
        }
        return hitBankRows;
    }   
    
    public static boolean passCalibCuts(DataEvent event, DataBank bnkHits, int i, Map<Integer, SegmentProperty> segPropMap) {
        boolean pass = false;
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;
        double beta = (double) bnkHits.getFloat("beta", i);
        //if(beta>1) beta=1;
        if (bnkHits.getByte("trkID", i) >0 
                    && beta>= Double.parseDouble(T2DViewer.betaCut.getText()) 
                    && beta<= Double.parseDouble(T2DViewer.betaCut2.getText()) 
                    && selectOnAlpha(superlayer, alphaRadUncor)==true
                    && bnkHits.getFloat("TFlight", i)>0 
                    //&& segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()!=Math.abs(segPropMap.get(bnkHits.getInt("clusterID", i)).getSumLR())
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getNumWireWithinDW()<=Integer.parseInt(T2DViewer.npassWires.getText())
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()>Integer.parseInt(T2DViewer.nWires.getText())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && passPID(event, bnkHits, i)==true
                )
            {   
                avebeta+=beta;
                betacnt++;
                pass = true;
            }
        return pass;
    } 

     
    public static boolean passCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        
        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("TFlight", i)>0 
                    //&& segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()!=Math.abs(segPropMap.get(bnkHits.getInt("clusterID", i)).getSumLR())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && passPID(event, bnkHits, i)==true)
            {
                pass = true;
            }
        return pass;
    }
    
    
    private static boolean passPID(DataEvent event, DataBank bnkHits, int rowIdxinTrkBank) {
        boolean pass = false;
        int trkID = bnkHits.getByte("trkID", rowIdxinTrkBank);
        if(polarity==0) return isElectronZeroField(event, trkID);
        int pid = readPID(event, trkID);
        //pass if the track is identified as an electron or as a hadron
        //if(pid==11 || Math.abs(pid)==211 || Math.abs(pid)==2212 || Math.abs(pid)==321) {
        if(Integer.parseInt(T2DViewer.pid.getText())==-1) {
            pass=true;
        }
        if(pid==Integer.parseInt(T2DViewer.pid.getText())) {    
            pass = true;
        }
        
        return pass;
    }
    
    
    private static boolean selectOnAlpha(int superlayer, double alphaRadUncor) {
        boolean pass = false;
        int i = superlayer - 1;
        if(alphaRadUncor>Double.parseDouble(T2DViewer.alphaCuts1[i].getText()) &&
                alphaRadUncor<Double.parseDouble(T2DViewer.alphaCuts2[i].getText())) {
            pass = true;
        }
        return pass;        
    }
    public static void getSegProperty(DataBank bnkHits, Map<Integer, SegmentProperty> segPropMap, 
            Map<Integer, ArrayList<Integer>> segMapTBHits) {
        
        segMapTBHits.clear();
        segPropMap.clear();
               
        for (int j = 0; j < bnkHits.rows(); j++){
            Integer cID = bnkHits.getInt("clusterID", j);
            if(segMapTBHits.containsKey(cID)==false) {
                segMapTBHits.put(cID, new ArrayList<>());
            }
            
            segMapTBHits.get(cID).add((int)(bnkHits.getShort("wire", j)*bnkHits.getByte("LR", j)));
            
        }
        
        Iterator<Map.Entry<Integer, ArrayList<Integer>>> itr = segMapTBHits.entrySet().iterator(); 
          
        while(itr.hasNext()) { 
            Map.Entry<Integer, ArrayList<Integer>> entry = itr.next(); 
            segPropMap.put(entry.getKey() , 
                    new SegmentProperty(entry.getKey(),entry.getValue(),Integer.parseInt(T2DViewer.deltaWire.getText())));
        } 

    }
    public static class CalHit {

        public int wire;
        public int layer;
        public int idx;

        public CalHit() {
        }
    }
    
}
