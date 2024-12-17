/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;


import java.io.FileNotFoundException;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.alphaBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.errs;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.parNames;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
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

    public static void fixPar(int i, MinuitPar mp) {
        if(mp.fixed[i]==true) {
        } else {
            mp.fix(i);
            mp.fixed[i]=true;
        }
    }
    public static void releasePar(int i, MinuitPar mp) {
        if(mp.fixed[i]==true) {
            mp.release(i);
            mp.fixed[i]=false;
        } else {
        }
    }

    public static void updatePar(MinuitPar mp, MnUserParameters userParameters) {
        for(int p =0; p<11; p++) {
            mp.setValue(p, userParameters.value(p));
            mp.setError(p, userParameters.error(p));
        }
    }
    
    public static void initParsForFit(Map<Coordinate, MinuitPar> TvstrkdocasFitPars, 
            boolean[][][] fixFit) {
        for(int s =0; s<7; s++) {
            int fitFitS =s;
            if(s==6) fitFitS=0;
            for(int i =0; i<6; i++) {
                MinuitPar TvstrkdocasFitParsClon = TvstrkdocasFitPars.get(new Coordinate(s,i));
                System.out.println(s+"] FOR FIT ...."+TvstrkdocasFitParsClon.toString());
                for (int p = 0; p < 10; p++) { 
                    TvstrkdocasFitParsClon.fixed[p]=false; 
                    if(fixFit[p][i][fitFitS]==true) {
                        FitUtility.fixPar(p, TvstrkdocasFitParsClon);
                    }
                }
                TvstrkdocasFitParsClon.fixed[10]=false; 
                FitUtility.fixPar(10, TvstrkdocasFitParsClon);
                //Fix r, distbeta//
                FitUtility.fixPar(2, TvstrkdocasFitParsClon);
                FitUtility.fixPar(4, TvstrkdocasFitParsClon);

                 System.out.println("FOR FIT ...."+TvstrkdocasFitParsClon.toString());
            }
        }    
    }
    public static void releaseParsAfterFit(Map<Coordinate, MinuitPar> TvstrkdocasFitPars) {
        for(int s =0; s<7; s++) {
            for(int i =0; i<6; i++) {
                MinuitPar TvstrkdocasFitParsClon = TvstrkdocasFitPars.get(new Coordinate(s,i));
                 System.out.println(s+"] A FIT ...."+TvstrkdocasFitParsClon.toString());
                for (int p = 0; p < 11; p++) {
                    FitUtility.releasePar(p, TvstrkdocasFitParsClon);
                }
            }
        }    
    }
    
    int maxNfits = 10;
    public void loadFitPars(Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, MinuitPar> TvstrkdocasFitPars, 
                            Map<Coordinate, FitLine> TvstrkdocasFits,
                            Map<Coordinate, H1F> ParsVsIter, 
                            double[][] resetPars, String[] parNames,
                            CalibrationConstants calib) throws FileNotFoundException {
        System.out.println("UPDATING TABLE.....");
        for (int s = 0; s<7; s++) {
            for (int i = 0; i < 6; i++) {
                double[] pars = new double[11];
                int si = s;
                if(s==6) si = 0; //for all sectors use sector 1
                //T2DFunctions.polyFcnMac(x, alpha, bfield, v0[s][r], vmid[s][r], FracDmaxAtMinVel[s][r], 
                //tmax, dmax, delBf, Bb1, Bb2, Bb3, Bb4, superlayer) ;
                pars[0] = org.jlab.rec.dc.timetodistance.TableLoader.v0[si][i];
                pars[1] = org.jlab.rec.dc.timetodistance.TableLoader.vmid[si][i];
                pars[2] = org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[si][i];
                pars[3] = org.jlab.rec.dc.timetodistance.TableLoader.Tmax[si][i];
                pars[4] = org.jlab.rec.dc.timetodistance.TableLoader.distbeta[si][i];
                pars[5] = org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[si][i];
                pars[6] = org.jlab.rec.dc.timetodistance.TableLoader.b1[si][i];
                pars[7] = org.jlab.rec.dc.timetodistance.TableLoader.b2[si][i];
                pars[8] = org.jlab.rec.dc.timetodistance.TableLoader.b3[si][i];
                pars[9] = org.jlab.rec.dc.timetodistance.TableLoader.b4[si][i];
                pars[10] = 2.*Constants.getInstance().wpdist[i];//fix dmax

                resetPars[i] = pars;
                TvstrkdocasFitPars.put(new Coordinate(s,i), new MinuitPar());
                for(int p = 0; p < 10; p++) {
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).add(parNames[p], pars[p], errs[p]);
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).fixed[p]=false;
                    //create graphs of parameters for various iterations
                    if(s==0) ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p, "superlayer "+(i+1)+" par "+p,maxNfits+1, 0.5,maxNfits+1.5));
                }
                TvstrkdocasFitPars.get(new Coordinate(s,i)).add(parNames[10], pars[10], errs[10]);
                
                for (int j = 0; j < alphaBins; j++) {
                    if(i<2 || i>3) {
                        if(Tvstrkdocas.get(new Coordinate(s, i, j, BBins))!=null && 
                                Tvstrkdocas.get(new Coordinate(s, i, j, BBins)).getProfileX().getMax()>0) {
                            updateTable(s, i, j, calib, TvstrkdocasFitPars);
                        }
                    } else {
                        int c=0;
                        for(int k = 0; k < BBins; k++) { 
                            if(Tvstrkdocas.get(new Coordinate(s, i, j, k))!=null && 
                                Tvstrkdocas.get(new Coordinate(s, i, j, k)).getProfileX().getMax()>0) {
                                c++; 
                            }
                        }
                        if(c>1) updateTable(s, i, j, calib, TvstrkdocasFitPars);
                    }
                }
            }
        }
        for(int s=0; s<6;s++) {
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < alphaBins; j++) {
                    if(i<2 || i>3) {
                        TvstrkdocasFits.put(new Coordinate(s, i,j,BBins), new FitLine("f"+""+s+""+i+""+j+"0", i, j, BBins, 
                        TvstrkdocasFitPars.get(new Coordinate(s,i))));

                    } else {
                        for(int k = 0; k < BBins; k++) { 
                            TvstrkdocasFits.put(new Coordinate(s,i,j,k), new FitLine("f"+""+s+""+i+""+j+""+k, i, j, k, 
                            TvstrkdocasFitPars.get(new Coordinate(s,i))));

                        }
                    }
                }
            }
        }
        System.out.println("TABLE UPDATED");
    }
        
    public static void updateTable(int secidx, int slidx, int j, CalibrationConstants calib, Map<Coordinate, MinuitPar> TvstrkdocasFitPars) {
       
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(0), "v0", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(1), "vmid", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(2), "R", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(3), "tmax", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(4), "distbeta", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(5), "delBf", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(6), "b1", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(7), "b2", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(8), "b3", secidx+1, slidx+1, j+1);
        calib.setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(secidx, slidx)).value(9), "b4", secidx+1, slidx+1, j+1);
   
    
    }  

    public static void reLoadFitPars(Map<Coordinate, H1F> ParsVsIter, Map<Coordinate, MinuitPar> TvstrkdocasFitPars) {
        for (int s =0; s < 6; s++) {
            for (int i = 0; i < 6; i++) {
                //System.out.println("s="+s+" i0="+i0);
                org.jlab.rec.dc.timetodistance.TableLoader.v0[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(0);
                org.jlab.rec.dc.timetodistance.TableLoader.vmid[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(1);
                org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(2);
                //System.out.println("TMAX "+org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i]+" ==> ");
                org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(3);
                //System.out.println(".................... "+org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i]+"");
                org.jlab.rec.dc.timetodistance.TableLoader.distbeta[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(4);
                org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(5);
                org.jlab.rec.dc.timetodistance.TableLoader.b1[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(6);
                org.jlab.rec.dc.timetodistance.TableLoader.b2[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(7);
                org.jlab.rec.dc.timetodistance.TableLoader.b3[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(8);
                org.jlab.rec.dc.timetodistance.TableLoader.b4[s][i] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(9);
                for(int p = 0; p<6; p++) {
                    ParsVsIter.get(new Coordinate(i,p)).setBinContent(0, TvstrkdocasFitPars.get(new Coordinate(s,i)).value(p));
                    ParsVsIter.get(new Coordinate(i,p)).setBinError(0, TvstrkdocasFitPars.get(new Coordinate(s,i)).error(p));
                }
            }
        }
        
        TableLoader.ReFill();
    }

    void resetPars(Map<Coordinate, MinuitPar> TvstrkdocasFitPars, Map<Coordinate, H1F> ParsVsIter, String[] parNames, 
                   double[][]resetPars) {
        TvstrkdocasFitPars.clear();
        for(int s = 0; s<7; s++) {
            for (int i = 0; i < 6; i++) {
                double[] pars = resetPars[i];
                TvstrkdocasFitPars.put(new Coordinate(s,i), new MinuitPar());
                for(int p = 0; p < 10; p++) {
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).add(parNames[p], pars[p]);
                    //create graphs of parameters for various iterations
                   if(s==0) ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p+": " +parNames[p]+", superlayer "+(i+1),this.maxNfits+1, 0.5,this.maxNfits+1.5));
                }
                TvstrkdocasFitPars.get(new Coordinate(s,i)).add(parNames[10], pars[10]);
            }
        } 
    }

    public void runParamScan(boolean fixFit[][][], 
            Map<Coordinate, MinuitPar> TvstrkdocasFitPars, 
            Map<Coordinate, FitFunction> TvstrkdocasFit, Map<Coordinate, GraphErrors> TvstrkdocasProf) {
        for(int s = T2DCalib.minSec; s<T2DCalib.maxSec; s++) {
            this.runParamScan(fixFit, s, TvstrkdocasFitPars, TvstrkdocasFit, TvstrkdocasProf);
            T2DViewer.voice.speak("Parameter Scan done for Sector "+(s+1));
        }
    }
    private void runParamScan(boolean fixFit[][][], int s, 
            Map<Coordinate, MinuitPar> TvstrkdocasFitPars, 
            Map<Coordinate, FitFunction> TvstrkdocasFit, Map<Coordinate, GraphErrors> TvstrkdocasProf) {
       
        MnMigrad scanner[] = new MnMigrad[6];
        MnMigrad fitter[] = new MnMigrad[6];
        
        for(int i =0; i<6; i++) {
            System.out.println("******************************************FIXING PAR10 FOR i = "+i+ " s = "+s+" "
            +TvstrkdocasFitPars.get(new Coordinate(s,i)).toString());
            MinuitPar TvstrkdocasFitParsClon = TvstrkdocasFitPars.get(new Coordinate(s,i)).clone();
            
            TvstrkdocasFit.put(new Coordinate(s,i), 
                                         new FitFunction(s, i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
            scanner[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(s, i)), 
                                                TvstrkdocasFitParsClon,0);
            fitter[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(s, i)), 
                                                TvstrkdocasFitParsClon,1);
            for(int p =0; p<11; p++) {
                if(TvstrkdocasFitPars.get(new Coordinate(s, i)).fixed[p]) {
                    scanner[i].fix(p);
                    //fitter[i0].fix(p);
                }
            }
        }
        
        double[][] pars2 = T2DFitter.estimateFixedParsPerRegion(fixFit, scanner, fitter, s);
        for(int i = 0; i<3; i++) {
            double[] pars = new double[2];
            pars[0] = pars2[0][i];
            pars[1] = pars2[1][i];
           // System.out.println(i+"] PARAMETERS R "+pars2[0][i]+" distbeta "+pars2[1][i]);
            T2DFitter.fitWithFixedParsPerRegion(s, i,fixFit, pars, scanner, fitter, TvstrkdocasFitPars);
        }
        
    }

    void runFit(boolean[][][] fixFit, Map<Coordinate, MinuitPar> TvstrkdocasFitPars, 
            Map<Coordinate, FitFunction> TvstrkdocasFit, 
            Map<Coordinate, GraphErrors> TvstrkdocasProf) {
        MnMigrad scanner[] = new MnMigrad[6*7];
        MnMigrad fitter[] = new MnMigrad[6*7];
        double[][][] pars = new double[7][6][2];
        
        
        for(int s =0; s<7; s++) {
            for(int i =0; i<6; i++) {
                pars[s][i][0] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(2);
                pars[s][i][1] = TvstrkdocasFitPars.get(new Coordinate(s,i)).value(4);
                TvstrkdocasFit.put(new Coordinate(s,i), 
                                             new FitFunction(s,i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
                scanner[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(s,i)), 
                                                    TvstrkdocasFitPars.get(new Coordinate(s,i)),0);
                fitter[i] = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(s,i)), 
                                                    TvstrkdocasFitPars.get(new Coordinate(s,i)),1);
            }
        }
        
        for(int s =T2DCalib.minSec; s<T2DCalib.maxSec; s++) {
            T2DFitter.fitWithFixedPars(pars, scanner, fitter, s);
        }
    }
    

    
    public static class MinuitPar extends MnUserParameters implements Cloneable {

        public boolean[] fixed ;

        public MinuitPar() {
            fixed= new boolean[11];
        }
        
        protected MinuitPar clone() {
            MinuitPar cloned = new MinuitPar();
            for(int p = 0; p<11; p++) {
                cloned.add(parNames[p], this.value(p), this.error(p));
                cloned.fixed=this.fixed;
            }
            return cloned;
        }

    }

    
}
