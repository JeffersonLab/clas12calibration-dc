/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
/**
 *
 * @author ziegler
 */
public class T2DFitter  {
    
    
    static double maxIter = 100;
    //get R, distbeta, delBf, for best combine ch2 for all superlayers
    
    public static double[][] estimateFixedParsPerRegion(boolean fixFit[][][], MnMigrad scanner[], MnMigrad fitter[], 
            int sec) { 
        double errs2 = T2DCalib.errs[2];
        double errs4 = T2DCalib.errs[4];
        int nR = (int) Math.ceil((T2DCalib.Rlimits[0][1]-T2DCalib.Rlimits[0][0])/(double)errs2);
        int ndbeta = (int) Math.ceil((T2DCalib.distbetalimits[0][1]-T2DCalib.distbetalimits[0][0])/(double)errs4);
       
        double[] bestR = new double[3];
        double[] bestDistbeta = new double[3];

        double bestchi2[] = new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};

        double R ;
        double distbeta ;
        System.out.println("Estimating PARS "+nR+" "+ndbeta);
        int cnt =0;
        
        for(int ri =1; ri<nR-1; ri++) {
            for(int di =1; di<ndbeta-1; di++) { 
                R=T2DCalib.Rlimits[0][0]+(double)ri*errs2;
                distbeta=T2DCalib.distbetalimits[0][0]+(double)di*errs4;
                double[] c2=new double[3];
                cnt++;
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta);
                for(int i =0; i<2; i++) {
                    String s="";
                    s+=(" ***********************************************************≈***");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1))+" SECTOR "+(sec+1);
                    s+=(" ***************************************************************");
                    fMin fm = getfMinFixedRDPars(sec, i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2[0]+= fm.getChi2();
                }
                for(int i =2; i<4; i++) {
                    String s="";
                    s+=(" ***********************************************************≈***");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1))+" SECTOR "+(sec+1);
                    s+=(" ***************************************************************");
                    fMin fm = getfMinFixedRDPars(sec, i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2[1]+= fm.getChi2();
                }
                for(int i =4; i<6; i++) {
                    String s="";
                    s+=(" ***********************************************************≈***");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1))+" SECTOR "+(sec+1);
                    s+=(" ***************************************************************");
                    fMin fm = getfMinFixedRDPars(sec, i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2[2]+= fm.getChi2();
                }
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta +" c2 "+c2[0]+" "+c2[1]+" "+c2[2]);
                for(int j = 0; j < 3; j++) {
                    if(c2[j]<bestchi2[j]) {
                        bestR[j] = R;
                        bestDistbeta[j] = distbeta;
                        bestchi2[j] = c2[j];
                    }
                }
            }
        }
        double[][] result = new double[2][3];
        for(int k =0; k<3; k++) {
            result[0][k] = bestR[k];
            result[1][k] = bestDistbeta[k];
            System.out.println("Best results for R "+bestR[k]+" istbeta "+bestDistbeta[k]);
        }
        return result;
    }
    
    public static double[] estimateFixedPars(boolean fixFit[][][], MnMigrad scanner[], MnMigrad fitter[], int sec) { 
        double errs2 = T2DCalib.errs[2];
        double errs4 = T2DCalib.errs[4];
        int nR = (int) Math.ceil((T2DCalib.Rlimits[0][1]-T2DCalib.Rlimits[0][0])/(double)errs2);
        int ndbeta = (int) Math.ceil((T2DCalib.distbetalimits[0][1]-T2DCalib.distbetalimits[0][0])/(double)errs4);
       
        double bestR = T2DCalib.Rlimits[0][0];
        double bestDistbeta = T2DCalib.distbetalimits[0][0];

        double bestchi2 = Double.POSITIVE_INFINITY;

        double R = T2DCalib.Rlimits[0][0];
        double distbeta = T2DCalib.distbetalimits[0][0];
        System.out.println("Estimating PARS "+nR+" "+ndbeta);
        int cnt =0;
        
        for(int ri =1; ri<nR-1; ri++) {
            for(int di =1; di<ndbeta-1; di++) { 
                R=T2DCalib.Rlimits[0][0]+(double)ri*errs2;
                distbeta=T2DCalib.distbetalimits[0][0]+(double)di*errs4;
                double c2=0;
                cnt++;
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta);
                for(int i =0; i<6; i++) {
                    String s="";
                    s+=(" ******************************************");
                    s+=("   RUNNING THE PARAMETER SCAN FOR SUPERLAYER "+(i+1));
                    s+=(" ******************************************");
                    fMin fm = getfMinFixedRDPars(sec, i, fixFit, scanner[i], fitter[i], R, distbeta, true, s);
                    c2+= fm.getChi2();
                }
                System.out.println(cnt+"] R "+R+" disbeta "+distbeta +" c2 "+c2);
                if(c2<bestchi2) {
                    bestR = R;
                    bestDistbeta = distbeta;
                    bestchi2 = c2;
                    System.out.println(cnt+"] best R "+R+" disbeta "+distbeta +" c2 "+c2);
                }
            }
        }
        
        return new double[] {bestR, bestDistbeta};
    }
    
    private static fMin getfMinFixedRDPars(int sec, int i, boolean fixFit[][][], MnMigrad scanner, MnMigrad fitter, 
            double R, double distbeta, boolean reset, String s) {
        
        double edm = Double.POSITIVE_INFINITY;
        double edm2 = Double.POSITIVE_INFINITY;
        double bestchi2 = Double.POSITIVE_INFINITY;
        double bestMchi2 = Double.POSITIVE_INFINITY;
        
        int fitFitS =sec;
        if(sec==6) fitFitS=0;
        
        System.out.println(s); 
        FunctionMinimum min = null ;
        FunctionMinimum bestmin = null ;
        FunctionMinimum min2 = null ;
        FunctionMinimum bestmin2 = null ;
        
        
        for(int pi = 0; pi<3; pi++) {
            scanner.setLimits(pi, T2DCalib.limits[pi][i][0], T2DCalib.limits[pi][i][1]);
            fitter.setLimits(pi, T2DCalib.limits[pi][i][0], T2DCalib.limits[pi][i][1]);
        }

        for(int pi = 4; pi<6; pi++) {
            scanner.setLimits(pi, T2DCalib.limits[pi-1][i][0], T2DCalib.limits[pi-1][i][1]);
            fitter.setLimits(pi, T2DCalib.limits[pi-1][i][0], T2DCalib.limits[pi-1][i][1]);
        }
        System.out.println("tmax"+i+" = "+org.jlab.rec.dc.timetodistance.TableLoader.Tmax[0][i]);
        scanner.setLimits(3, org.jlab.rec.dc.timetodistance.TableLoader.Tmax[0][i]-50, org.jlab.rec.dc.timetodistance.TableLoader.Tmax[0][i]+50);
        fitter.setLimits(3, org.jlab.rec.dc.timetodistance.TableLoader.Tmax[0][i]-50, org.jlab.rec.dc.timetodistance.TableLoader.Tmax[0][i]+50);
        
        
        try {
            min = scanner.minimize();
        } catch (Exception e) {
            // Handle the exception appropriately
            System.err.println("An error occurred during minimization: " + e.getMessage());
            // You may want to log the exception or take other actions depending on your application
        }
        if(fixFit[2][i][fitFitS]==false) {
            min.userParameters().setValue(2,R);
            scanner.fix(2);
            fitter.fix(2);
        }
        
        if(fixFit[4][i][fitFitS]==false) {
            min.userParameters().setValue(4,distbeta);
            scanner.fix(4);
            fitter.fix(4);
        }
        int itercnt=0;
        for(int it = 0; it<maxIter; it++) {
                try {
                        min = scanner.minimize();
                    } catch (Exception e) {
                        // Handle the exception appropriately
                        System.err.println("An error occurred during minimization: " + e.getMessage());
                        // You may want to log the exception or take other actions depending on your application
                    }
                itercnt++;
                if(FitFunction.chi2<bestchi2) {
                    bestchi2 = FitFunction.chi2;
                    bestmin = min;
                    
                }
                if(edm-FitFunction.chi2<0.1 || FitFunction.chi2+10>edm) break;
                edm = FitFunction.chi2;
        } 
        for (int p = 0; p < 10; p++) {
            fitter.setValue(p, bestmin.userParameters().value(p));
        }
        
        int itercnt2=0;
        for(int it = 0; it<maxIter; it++) {
            
            try {
                    min2 = fitter.minimize();
                } catch (Exception e) {
                    // Handle the exception appropriately
                    System.err.println("An error occurred during minimization: " + e.getMessage());
                    // You may want to log the exception or take other actions depending on your application
                }
            itercnt2++;
            
            if(FitFunction.chi2<bestMchi2) {
                bestMchi2 = FitFunction.chi2;
                bestmin2 = min2;
                if(edm2-FitFunction.chi2<0.01) break;
                edm2 = FitFunction.chi2;
            }
        }
        if(bestmin2==null || bestMchi2>bestchi2) {
            bestMchi2 = bestchi2;
            bestmin2 = bestmin;
        }
        
        if(fixFit[2][i][fitFitS]==false) {
            scanner.release(2);
            fitter.release(2);
        }
        if(fixFit[4][i][fitFitS]==false) {
            scanner.release(4);
            fitter.release(4);
        }
       
        if(reset) {
            for (int p = 0; p < 10; p++) {
                scanner.setValue(p, T2DCalib.TvstrkdocasFitPars.get(new Coordinate(i+6*sec)).value(p));
                fitter.setValue(p, T2DCalib.TvstrkdocasFitPars.get(new Coordinate(i+6*sec)).value(p));
            }
        }
        System.out.println(+itercnt+"] SCAN CHI2 "+bestchi2);
        System.out.println(itercnt2+"] MIGRAD CHI2 "+bestMchi2);
        
        System.gc();
        
        return new fMin(min2, bestMchi2);
    }
    
    
    
    public static void fitWithFixedParsPerRegion(int sec, int ridx,boolean fixFit[][][], double pars[], MnMigrad scanner[], MnMigrad fitter[]) { 
        fMin results [] = new fMin[6];
       
        
        for(int i0 =2*ridx; i0<2*ridx+2; i0++) {
            String s2="";
            s2+=(" ******************************************");
            s2+=("   RUNNING THE PARAMETER FIT FOR SUPERLAYER "+(i0+1));
            s2+=(" ******************************************");
            fMin fm2 = getfMinFixedRDPars(sec, i0, fixFit, scanner[i0], fitter[i0], pars[0], pars[1], false, s2);
            FunctionMinimum fmin=null;
            if(fm2.getFcnMin().isValid()) {
                results[i0] = fm2;
                fmin = fm2.getFcnMin();
                System.out.println("UPDATED "+fmin.toString());
                T2DCalib.TvstrkdocasFitPars.put(new Coordinate(i0+6*sec),fmin.userParameters()); 
                if(sec==6) {
                    for(int s = 0; s<6; s++) {
                        T2DCalib.TvstrkdocasFitPars.put(new Coordinate(i0+6*s),fmin.userParameters()); 
                    }
                }
            } 
        }
    
    }
    public static void fitWithFixedPars(boolean fixFit[][][], double pars[], MnMigrad scanner[], MnMigrad fitter[], int sec) { 
        fMin results [] = new fMin[6];
        
        for(int i =0; i<6; i++) {
            
            String s2="";
            s2+=(" ******************************************");
            s2+=("   RUNNING THE PARAMETER FIT FOR SUPERLAYER "+(i+1));
            s2+=(" ******************************************");
            fMin fm2 = getfMinFixedRDPars(sec, i, fixFit, scanner[i], fitter[i], pars[0], pars[1], false, s2);
            FunctionMinimum fmin=null;
            if(fm2.getFcnMin().isValid()) {
                results[i] = fm2;
                fmin = fm2.getFcnMin();
                System.out.println("UPDATED "+fmin.toString());
                T2DCalib.TvstrkdocasFitPars.put(new Coordinate(i+6*sec),fmin.userParameters()); 
                if(sec==6) {
                    for(int s = 0; s<6; s++) {
                        T2DCalib.TvstrkdocasFitPars.put(new Coordinate(i+6*s),fmin.userParameters()); 
                    }
                }
            } 
        }
    }
    
    
    
}
