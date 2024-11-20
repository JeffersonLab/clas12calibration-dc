/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.io.FileNotFoundException;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.TvstrkdocasFitPars;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.alphaBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.errs;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.getSectorIdx;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.getSuperlayerIdx;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.rec.dc.Constants;

/**
 *
 * @author ziegler
 */
public class FitUtility {
    int maxNfits = 10;
    public void loadFitPars(Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, MnUserParameters> TvstrkdocasFitPars, 
                            Map<Coordinate, H1F> ParsVsIter, 
                            double[][] resetPars, String[] parNames,
                            CalibrationConstants calib) throws FileNotFoundException {
        System.out.println("UPDATING TABLE.....");
        for (int s = 0; s<6; s++) {
            for (int i0 = 0; i0 < 6; i0++) {
                double[] pars = new double[11];
                int i = i0+s*6;
                //T2DFunctions.polyFcnMac(x, alpha, bfield, v0[s][r], vmid[s][r], FracDmaxAtMinVel[s][r], 
                //tmax, dmax, delBf, Bb1, Bb2, Bb3, Bb4, superlayer) ;
                pars[0] = org.jlab.rec.dc.timetodistance.TableLoader.v0[s][i0];
                pars[1] = org.jlab.rec.dc.timetodistance.TableLoader.vmid[s][i0];
                pars[2] = org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[s][i0];
                pars[3] = org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i0];
                pars[4] = org.jlab.rec.dc.timetodistance.TableLoader.distbeta[s][i0];
                pars[5] = org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[s][i0];
                pars[6] = org.jlab.rec.dc.timetodistance.TableLoader.b1[s][i0];
                pars[7] = org.jlab.rec.dc.timetodistance.TableLoader.b2[s][i0];
                pars[8] = org.jlab.rec.dc.timetodistance.TableLoader.b3[s][i0];
                pars[9] = org.jlab.rec.dc.timetodistance.TableLoader.b4[s][i0];
                pars[10] = 2.*Constants.getInstance().wpdist[i0];//fix dmax

                resetPars[i] = pars;
                TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
                for(int p = 0; p < 10; p++) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                    //create graphs of parameters for various iterations
                    ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p, "superlayer "+(i+1)+" par "+p,maxNfits+1, 0.5,maxNfits+1.5));
                }
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10], errs[10]);
            //
                for (int j = 0; j < alphaBins; j++) {
                    if(i0<2 || i0>3) {
                        if(Tvstrkdocas.get(new Coordinate(i, j, BBins))!=null && 
                                Tvstrkdocas.get(new Coordinate(i, j, BBins)).getProfileX().getMax()>0) {
                            updateTable(i, j, calib, TvstrkdocasFitPars);
                        }
                    } else {
                        int c=0;
                        for(int k = 0; k < BBins; k++) { 
                            if(Tvstrkdocas.get(new Coordinate(i, j, k))!=null && 
                                Tvstrkdocas.get(new Coordinate(i, j, k)).getProfileX().getMax()>0) {
                                c++; 
                            }
                        }
                        if(c>1) updateTable(i, j, calib, TvstrkdocasFitPars);
                    }
                }
            }
        }
        
    }
        
    public static void updateTable(int i, int j, CalibrationConstants calib, Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) {
        int slidx = getSuperlayerIdx(i);
        int secidx = getSectorIdx(i);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(0), "v0", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(1), "vmid", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(2), "R", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(3), "tmax", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(4), "distbeta", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(5), "delBf", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(6), "b1", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(7), "b2", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(8), "b3", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(9), "b4", secidx+1, slidx+1, j+1);
   
    
    }  

    public static void reLoadFitPars(Map<Coordinate, H1F> ParsVsIter) {
        for (int s =0; s < 6; s++) {
            for (int i0 = 0; i0 < 6; i0++) {
                int i = i0+s*6;
                System.out.println("s="+s+" i0="+i0);
                org.jlab.rec.dc.timetodistance.TableLoader.v0[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(0);
                org.jlab.rec.dc.timetodistance.TableLoader.vmid[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(1);
                org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(2);
                System.out.println("TMAX "+org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i0]+" ==> ");
                org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(3);
                System.out.println(".................... "+org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i0]+"");
                org.jlab.rec.dc.timetodistance.TableLoader.distbeta[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(4);
                org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(5);
                org.jlab.rec.dc.timetodistance.TableLoader.b1[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(6);
                org.jlab.rec.dc.timetodistance.TableLoader.b2[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(7);
                org.jlab.rec.dc.timetodistance.TableLoader.b3[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(8);
                org.jlab.rec.dc.timetodistance.TableLoader.b4[s][i0] = TvstrkdocasFitPars.get(new Coordinate(i)).value(9);
                for(int p = 0; p<6; p++) {
                    ParsVsIter.get(new Coordinate(i,p)).setBinContent(0, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
                    ParsVsIter.get(new Coordinate(i,p)).setBinError(0, TvstrkdocasFitPars.get(new Coordinate(i)).error(p));
                }
            }
        }
        
        TableLoader.ReFill();
    }

    void resetPars(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars, Map<Coordinate, H1F> ParsVsIter, String[] parNames, 
                   double[][]resetPars) {
        TvstrkdocasFitPars.clear();
        for(int s = 0; s<7; s++) {
            for (int i0 = 0; i0 < 6; i0++) {
                int i = i0+s*6;
                double[] pars = resetPars[i];
                TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
                for(int p = 0; p < 10; p++) {
                    //TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                    TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p]);
                    //create graphs of parameters for various iterations
                   ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p+": " +parNames[p]+", superlayer "+(i+1),this.maxNfits+1, 0.5,this.maxNfits+1.5));
                }
                //TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10], errs[10]);
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[10], pars[10]);
            }
        } }

    public void runParamScan(boolean fixFit[][][], 
            Map<Coordinate, MnUserParameters> TvstrkdocasFitPars, 
            Map<Coordinate, FitFunction> TvstrkdocasFit, Map<Coordinate, GraphErrors> TvstrkdocasProf) {
        for(int s = T2DCalib.minSec; s<T2DCalib.maxSec; s++) {
            this.runParamScan(fixFit, s, TvstrkdocasFitPars, TvstrkdocasFit, TvstrkdocasProf);
        }
    }
    private void runParamScan(boolean fixFit[][][], int s, 
            Map<Coordinate, MnUserParameters> TvstrkdocasFitPars, 
            Map<Coordinate, FitFunction> TvstrkdocasFit, Map<Coordinate, GraphErrors> TvstrkdocasProf) {
        int fitFitS =s;
        if(s==6) fitFitS=0;
        MnMigrad scanner[] = new MnMigrad[6];
        MnMigrad fitter[] = new MnMigrad[6];
        
        for(int i0 =0; i0<6; i0++) {
            int i =i0+s*6;
            TvstrkdocasFitPars.get(new Coordinate(i)).fix(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i0][fitFitS]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).fix(p);
                }
            }
            TvstrkdocasFit.put(new Coordinate(i), 
                                         new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
            scanner[i0] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                TvstrkdocasFitPars.get(new Coordinate(i)),0);
            fitter[i0] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                TvstrkdocasFitPars.get(new Coordinate(i)),1);
        }
        
        
        //double[] pars = this.estimateFixedPars(fixFit, scanner, fitter);
        //System.out.println("PARAMETERS R "+pars[0]+" distbeta "+pars[1]);
        double[][] pars2 = T2DFitter.estimateFixedParsPerRegion(fixFit, scanner, fitter, s);
        for(int i = 0; i<3; i++) {
            double[] pars = new double[2];
            pars[0] = pars2[0][i];
            pars[1] = pars2[1][i];
            System.out.println(i+"] PARAMETERS R "+pars2[0][i]+" distbeta "+pars2[1][i]);
            T2DFitter.fitWithFixedParsPerRegion(s, i,fixFit, pars, scanner, fitter);
        }
        for(int i0 =0; i0<6; i0++) {
            int i =i0+s*6;
            TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
            for (int p = 0; p < 10; p++) {
                if(fixFit[p][i0][fitFitS]==true) {
                    TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
                }
            }
        }
    }

    void runFit(boolean[][][] fixFit, Map<Coordinate, MnUserParameters> TvstrkdocasFitPars, 
            Map<Coordinate, FitFunction> TvstrkdocasFit, 
            Map<Coordinate, GraphErrors> TvstrkdocasProf) {
        MnMigrad scanner[] = new MnMigrad[36];
        MnMigrad fitter[] = new MnMigrad[36];
        double[] pars = new double[2];
        pars[0] = TvstrkdocasFitPars.get(new Coordinate(0)).value(2);
        pars[1] = TvstrkdocasFitPars.get(new Coordinate(0)).value(4);
        
        int fixFitS = 0;
        
        for(int s =T2DCalib.minSec; s<T2DCalib.maxSec; s++) {
            for(int i0 =0; i0<6; i0++) {
                int i =i0+s*6;
                if(T2DCalib.maxSec<6) fixFitS=s;
                TvstrkdocasFitPars.get(new Coordinate(i)).fix(10);
                for (int p = 0; p < 10; p++) {
                    if(fixFit[p][i0][fixFitS]==true) {
                        TvstrkdocasFitPars.get(new Coordinate(i)).fix(p);
                    }
                }
                TvstrkdocasFit.put(new Coordinate(i), 
                                             new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
                scanner[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                    TvstrkdocasFitPars.get(new Coordinate(i)),0);
                fitter[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                                                    TvstrkdocasFitPars.get(new Coordinate(i)),1);

            }
        }
        
        for(int s =T2DCalib.minSec; s<T2DCalib.maxSec; s++) {
            T2DFitter.fitWithFixedPars(fixFit, pars, scanner, fitter, s);
            for(int i0 =0; i0<6; i0++) {
                int i =i0+s*6;
                TvstrkdocasFitPars.get(new Coordinate(i)).release(10);
                for (int p = 0; p < 10; p++) {
                    if(fixFit[p][i0][fixFitS]==true) {
                        TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
                    }
                }
            }
        }
    }
}
