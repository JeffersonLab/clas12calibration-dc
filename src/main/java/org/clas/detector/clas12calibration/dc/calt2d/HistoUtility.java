/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.calt2d.FitUtility.MinuitPar;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.AlphaBinHalfWidth;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.AlphaValues;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.alphaBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.docaBinWidth;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.maxx;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.maxy;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.nbinx;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.nbiny;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.timeBinWidth;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.utils.groups.IndexedList;

/**
 *
 * @author ziegler
 */
public class HistoUtility {
   
    public static void fill(DataEvent event, DataBank bnkHits, Map<Integer, FittedHit> hitmap, Map<Integer, FittedHit> calhitmap, 
            Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, H2F> Tvscalcdocas, Map<Coordinate, FitLine> TvstrkdocasFits, 
            Map<Coordinate, H2F> Tresvstrkdocas, Map<Coordinate, H1F> timeResiFromFile, 
            Map<Coordinate, H1F> A, Map<Coordinate, H1F> B,
            Map<Integer, SegmentProperty> segPropMap) {
        
        int allSecIdx = 6;

        for (int i = 0; i < bnkHits.rows(); i++) {
            // Pre-compute commonly used values
            double bFieldVal = bnkHits.getFloat("B", i);
            int superlayer = bnkHits.getInt("superlayer", i);
            int sector = bnkHits.getInt("sector", i);
            int slyrIdx = superlayer - 1;
            int secIdx = sector - 1;
            double theta0 = Math.toDegrees(Math.acos(1 - 0.02 * bFieldVal));
            double alphaUncor = bnkHits.getFloat("Alpha", i) + (double) T2DCalib.polarity * theta0;
            int alphaBin = CalUtility.getAlphaBin(alphaUncor);

            if (alphaBin == -1) continue;  // Skip if alphaBin is invalid

            // Create the current FittedHit
            FittedHit theHit = HitUtility.getHit(bnkHits, i);

            // First pass: Calibration hits (if not already processed)
            if (HitUtility.passCalibCuts(event, bnkHits, i, segPropMap)) {
                if (hitmap.get(theHit.get_Id()) == null) {
                    hitmap.put(theHit.get_Id(), theHit);
                }
            }

            // Second pass: Process calibrated hits
            if (HitUtility.passCuts(event, bnkHits, i) && calhitmap.get(theHit.get_Id()) == null) {
                calhitmap.put(theHit.get_Id(), theHit);

                if (!T2DCalib.refitSegs) {
                    // Common coordinate for all sector bins
                    Coordinate coordBase = new Coordinate(secIdx, slyrIdx, alphaBin, BBins);
                    Coordinate coordAllSec = new Coordinate(allSecIdx, slyrIdx, alphaBin, BBins);

                    // Calibration time calculation
                    double calibTime = (bnkHits.getInt("TDC", i) - bnkHits.getFloat("TProp", i)
                        - bnkHits.getFloat("TFlight", i) - bnkHits.getFloat("TStart", i)
                        - bnkHits.getFloat("T0", i));

                    // Fill the data in the maps
                    Tvstrkdocas.get(coordBase).fill(bnkHits.getFloat("trkDoca", i), calibTime);
                    Tvstrkdocas.get(coordAllSec).fill(bnkHits.getFloat("trkDoca", i), calibTime);
                    Tvscalcdocas.get(coordBase).fill(bnkHits.getFloat("doca", i), calibTime);
                    Tvscalcdocas.get(coordAllSec).fill(bnkHits.getFloat("doca", i), calibTime);

                    double yf = TvstrkdocasFits.get(coordBase).evaluate(bnkHits.getFloat("trkDoca", i));
                    if (!Double.isNaN(yf)) {
                        Tresvstrkdocas.get(coordBase).fill(bnkHits.getFloat("trkDoca", i), calibTime - yf);
                        Tresvstrkdocas.get(coordAllSec).fill(bnkHits.getFloat("trkDoca", i), calibTime - yf);
                    }

                    // Region 2 (superlayer between 3 and 4)
                    if (superlayer > 2 && superlayer < 5) {
                        int bBin = CalUtility.getBBin(bFieldVal);
                        Coordinate coordRegion2 = new Coordinate(secIdx, slyrIdx, alphaBin, bBin);
                        Coordinate coordAllSecRegion2 = new Coordinate(allSecIdx, slyrIdx, alphaBin, bBin);

                        Tvstrkdocas.get(coordRegion2).fill(bnkHits.getFloat("trkDoca", i), calibTime);
                        Tvstrkdocas.get(coordAllSecRegion2).fill(bnkHits.getFloat("trkDoca", i), calibTime);
                        Tvscalcdocas.get(coordRegion2).fill(bnkHits.getFloat("doca", i), calibTime);
                        Tvscalcdocas.get(coordAllSecRegion2).fill(bnkHits.getFloat("doca", i), calibTime);

                        double r2yf = TvstrkdocasFits.get(coordRegion2).evaluate(bnkHits.getFloat("trkDoca", i));
                        if (!Double.isNaN(r2yf)) {
                            Tresvstrkdocas.get(coordRegion2).fill(bnkHits.getFloat("trkDoca", i), calibTime - r2yf);
                            Tresvstrkdocas.get(coordAllSecRegion2).fill(bnkHits.getFloat("trkDoca", i), calibTime - r2yf);
                        }
                    }
                }

                // Fill uncalibrated plot
                Coordinate timeResiCoord = new Coordinate(superlayer - 1);
                timeResiFromFile.get(timeResiCoord).fill(bnkHits.getFloat("timeResidual", i));

                // Region 2 filling for alpha and b-field bins
                if (superlayer < 3 || superlayer > 4) {
                    A.get(new Coordinate(superlayer - 1, alphaBin, BBins)).fill(alphaUncor);
                }

                // Fill B values for superlayer 3 and 4
                if (superlayer == 3 || superlayer == 4) {
                    int bBin = CalUtility.getBBin(bFieldVal);
                    B.get(new Coordinate(superlayer - 3, alphaBin, bBin)).fill(bFieldVal);
                    A.get(new Coordinate(superlayer - 1, alphaBin, bBin)).fill(alphaUncor);
                }
            }
        }
    }



    public static void createHistos(Map<Coordinate, H1F> timeResi, Map<Coordinate, H1F> timeResiFromFile, 
            Map<Coordinate, H1F> timeResiNew, Map<Coordinate, H1F> fitResi, 
            Map<Coordinate, MinuitPar> TvstrkdocasFitPars, Map<Coordinate, GraphErrors> TvstrkdocasProf, 
            Map<Coordinate, GraphErrors> TvstrkdocasInit, Map<Coordinate, H2F> Tvstrkdocas, 
            Map<Coordinate, H2F> Tresvstrkdocas, Map<Coordinate, H2F> Tvscalcdocas, 
            Map<Coordinate, FitLine> TvstrkdocasFits, 
            Map<Coordinate, H1F> A, Map<Coordinate, H1F> B, Map<Coordinate, H1F> BAlphaBins, Map<Coordinate, H1F> timeResiB,
            IndexedList<DataGroup> dataGroup) {
        DataGroup td = new DataGroup(7,2);
        DataGroup tdp = new DataGroup(14,8);
        DataGroup cd = new DataGroup(7,2);
        DataGroup tr = new DataGroup(6,1);
        DataGroup trb = new DataGroup(6,1);
        DataGroup fr = new DataGroup(6,1);
        
        for(int l=0; l<6; l++) {
            nbinx[l] = (int) Math.ceil(maxx[l]/docaBinWidth);
            nbiny[l] = (int) Math.ceil(maxy[l]/timeBinWidth);
        }
        
        int ijk = 0;
        int ij = 0;
        for (int i = 0; i < 6; i++) {
            timeResi.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            timeResiFromFile.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            timeResiNew.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
            fitResi.put(new Coordinate(i), new H1F("fit residual for sly " + (i+1), 100, -0.5, 0.5));
            
            tr.addDataSet(timeResi.get(new Coordinate(i)), i);
            tr.addDataSet(timeResiFromFile.get(new Coordinate(i)), i);
            tr.addDataSet(timeResiNew.get(new Coordinate(i)), i);
            fr.addDataSet(fitResi.get(new Coordinate(i)), i);
        }
        for(int s = 0; s<7; s++) {
            for (int i = 0; i < 6; i++) { 
                DataGroup prfdvst = new DataGroup(1,1);
                TvstrkdocasFitPars.put(new Coordinate(s, i), new MinuitPar());
                for (int j = 0; j < alphaBins; j++) {
                    for (int k = 0; k < BBins+1; k++) {
                        if((i<2 || i>3) && k<BBins ) continue; 
                        
                        String stg = "";
                        stg+= "sl" + (i + 1);
                        stg+=  ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")"
                                +", B "+k;
                        System.out.println(stg);
                        TvstrkdocasProf.put(new Coordinate(s,i,j,k), new GraphErrors());
                        TvstrkdocasInit.put(new Coordinate(s,i,j,k), new GraphErrors());
                        TvstrkdocasProf.get(new Coordinate(s,i,j,k)).setMarkerColor(k+1);
                        TvstrkdocasInit.get(new Coordinate(s,i,j,k)).setMarkerColor(k+1);
                        TvstrkdocasInit.get(new Coordinate(s,i,j,k)).setMarkerStyle(2);
                        TvstrkdocasInit.get(new Coordinate(s,i,j,k)).setTitle( "superlayer" + (i + 1)
                                + ", alpha ("+(AlphaValues[j]-AlphaBinHalfWidth)+", "+(AlphaValues[j]+AlphaBinHalfWidth)+")");
                        
                        Tvstrkdocas.put(new Coordinate(s,i,j,k), new H2F("trkDocavsT" + stg, nbinx[i], 0, maxx[i], nbiny[i], 0, maxy[i]));
                        Tresvstrkdocas.put(new Coordinate(s,i,j,k), new H2F("trkDocavsTres" +stg, nbinx[i], 0, maxx[i], nbiny[i], -50, 50));
                        Tvscalcdocas.put(new Coordinate(s,i,j,k), new H2F("calcDocavsT" + stg, nbinx[i], 0, maxx[i], nbiny[i], 0, maxy[i]));
                        
                        tdp.addDataSet(TvstrkdocasProf.get(new Coordinate(s,i,j,k)), ijk);
                        tdp.addDataSet(Tresvstrkdocas.get(new Coordinate(s,i,j,k)), 0);
                        //trkdvst.addDataSet(Tvstrkdocas.get(new Coordinate(i,j,k)), 0);
                        prfdvst.addDataSet(TvstrkdocasProf.get(new Coordinate(s,i,j,k)), 0);
                        prfdvst.addDataSet(TvstrkdocasInit.get(new Coordinate(s,i,j,k)), 0);

                        TvstrkdocasFits.put(new Coordinate(s,i,j,k), new FitLine());
                        prfdvst.addDataSet(TvstrkdocasFits.get(new Coordinate(s,i,j,k)), 0);
                        dataGroup.add(prfdvst, s+1, i+1, j+1);
                        dataGroup.add(tdp, s+1, i+1, j+1);
                        ijk++; 
                    }
                }
            }
        }
        
        //Alpha centroids
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < alphaBins; j++) {
                for (int k = 0; k <= BBins+1; k++) {
                    A.put(new Coordinate(i,j,k), new H1F("A centroid " +(i + 1)*1000+(j+1)+26, 100, -36, 36));
                }
            }
        }
        //B centroids
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < alphaBins; j++) {
                BAlphaBins.put(new Coordinate(i,j), new H1F("B  " , 100, 0.0, 3.0));
                for (int k = 0; k < BBins; k++) {
                    B.put(new Coordinate(i,j,k), new H1F("B centroid " +(i + 1)*1000+(j+1)+26, 100, 0.0, 3.0));
                }
            }
        }
        int ik=0;
        for (int i = 2; i < 4; i++) {
            for (int k = 0; k < BBins; k++) {
                timeResiB.put(new Coordinate(i,k), new H1F("time residual for sly " + (i+1), 100, -0.3, 0.3)); 
                trb.addDataSet(timeResiB.get(new Coordinate(i,k)), ik++);
                timeResiB.get(new Coordinate(i,k)).setLineColor(k+1);
            }
        }
        dataGroup.add(td, 0,0,0);
        dataGroup.add(tdp,1,0,0);
        dataGroup.add(cd, 2,0,0);
        dataGroup.add(tr, 3,0,0);
        dataGroup.add(trb, 4,0,0);
        dataGroup.add(fr, 5,0,0);
        System.gc();
    }
    
    
}
