/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.ArrayList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.rec.dc.Constants;

/**
 *
 * @author ziegler
 */
public class GraphUtility {
    
     public static void filltrkDocavsTGraphs(int s, int i, int j,  Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, GraphErrors> TvstrkdocasProf, boolean useBProf) {
        if(i<2 || i>3) { //region 1 and 3
            try {
                filltrkDocavsTGraphs(s, i, j, BBins, Tvstrkdocas, TvstrkdocasProf,useBProf);
            } catch(Exception e){
	       e.printStackTrace();
            }
        } else {
            for(int k = 0; k < BBins; k++) {    
                try{
                    filltrkDocavsTGraphs(s, i, j, k, Tvstrkdocas, TvstrkdocasProf,useBProf);
                } catch(Exception e){
	       e.printStackTrace();
            }
            }
        }       
    }
     
    static int MINENTRIES = 100;
    static double CHI2CUT = 5;
    static F1D f3 = new F1D("f3","[amp1]*gaus(x,[mean],[sigma1])+[amp2]*landau(x,[mean],[sigma2])", 0, 1.8);
    //F1D f2 = new F1D("f2","[amp1]*gaus(x,[mean1],[sigma1])+[amp2]*gaus(x,[mean2],[sigma2])+[p02]", 0, 1.8);
    static F1D f1 = new F1D("f1","[amp1]*gaus(x,[mean1],[sigma1])", 0, 2);
    
    public static void filltrkDocavsTGraphs(int s, int i, int j, int k,  Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, GraphErrors> TvstrkdocasProf, boolean useBProf) {
        Logger.getLogger("org.freehep.math.minuit").setLevel(Level.WARNING);
        if(TvstrkdocasProf.get(new Coordinate(s, i, j, k))!=null) {
            TvstrkdocasProf.get(new Coordinate(s, i, j, k)).reset();
        } 
        if(Tvstrkdocas.get(new Coordinate(s, i, j, k))!=null) {
            H2F h2 = Tvstrkdocas.get(new Coordinate(s, i, j, k));
            double integ = h2.getEntries(); 
            ArrayList<H1F> hslice = h2.getSlicesX();
            int n=0;
            boolean saveG =false;
            for(int si=0; si<hslice.size(); si++) {
                double x = h2.getXAxis().getBinCenter(si);
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                double dmax = 2.*Constants.getInstance().wpdist[i];
                if(hslice.get(si).getMean()==0 || integ<200 || x>dmax) 
                    continue;
                H1F hn = hslice.get(si).histClone("h");
                int nbrebin=0;
                while( nbrebin<4 && amp<MINENTRIES) {
                    hn= rebin(hn.histClone("hnc"));
                    amp = hn.getBinContent(hn.getMaximumBin());
                    nbrebin++;
                }
                n=0;
                if (amp > MINENTRIES) 
                    n = fitLandau(x, dmax, hn, TvstrkdocasProf.get(new Coordinate(s, i, j, k)));
                    
                if(n>0 && x/dmax >0.8) 
                    saveG = true;  //ensure the large doca entries are filled to avoid biases at small docas
            }
            if(!useBProf && (i==2 || i==3) && k>0)
                saveG = false;
            if(!saveG) {
                TvstrkdocasProf.get(new Coordinate(s, i, j, k)).reset();
            }
        }
    }

    static H1F rebin(H1F h) {
        int nbins = h.getData().length;
        double binWidth = h.getDataX(1)-h.getDataX(0);
        double lowEdge = h.getDataX(0)-binWidth/2.0;
        double highEdge = h.getDataX(nbins-1)+binWidth/2.0;
        int j =0;
        for(int i = 0; i<nbins-1; i++) {
            if(i%2==0) j++;
        } 
       
        H1F hn = new H1F("hn", "","", (int) j, lowEdge, highEdge);
        j =0;
        for(int i = 0; i<nbins-1; i++) {
            if(i%2==0) {
                double y1 = h.getBinContent(i);
                double y2 = h.getBinContent(i+1);
                double e1 = h.getBinError(i);
                double e2 = h.getBinError(i+1);
                double y = y1+y2;
                double e = Math.sqrt(e1*e1+e2*e2);
                hn.setBinContent(j, y);
                hn.setBinError(j, e);
                j++;
            }
            
        }
        
        return hn;
    }
//    private static int fitLandau(double x, double dmax, H1F hi, GraphErrors ge) {
//        
//        H1F h = hi.histClone("hc");
//        double meanh = h.getDataX(h.getMaximumBin());
//        double amp   = h.getBinContent(h.getMaximumBin());
//        double sigma = h.getRMS(); 
//        double mean = h.getMean();
//        double binSize = h.getDataX(1)-h.getDataX(0);
//        double min = meanh-2*sigma;
//        double max1 = meanh+2*sigma;
//        double max2 = meanh+3*sigma;
//        
//        if(min<0) min=0;
//        f1.reset(); 
//        f1.setRange(min, max1); 
//        f1.setParameter(0, amp);
//        f1.setParameter(1, mean);
//        f1.setParameter(2, sigma);
//        //f1.setParLimits(1, mean-2*binSize, mean-2*binSize);
//        boolean f1b=DataFitter.fit(f1, h,"Q"); //No options uses error for sigma 
//        if(!f1b) return 0;
//        double f1mu = f1.getParameter(1);
//        double ch2f1 = f1.getChiSquare();
//        f1.setParameter(1, meanh);
//       // f1.setParLimits(1, meanh-2*binSize, meanh+2*binSize);
//        DataFitter.fit(f1, h,"Q"); //No options uses error for sigma 
//        if(!f1b) return 0;
//        if(f1.getChiSquare()<ch2f1)
//           f1mu = f1.getParameter(1);
//        
//        //f3.reset(); 
//        f3.setRange(min, max2);
//        f3.setParameter(0, amp/3.0);
//        f3.setParameter(1, f1mu);
//        f3.setParameter(2, sigma);
//        f3.setParameter(3, amp);
//        f3.setParameter(4, sigma);
//        f3.setParLimits(0, 0, amp+10000);
//        f3.setParLimits(2, 0, 100);
//        f3.setParLimits(3, 0, amp+10000);
//        f3.setParLimits(4, 0, 100);
//        //f3.setParLimits(1, f1mu-sigma, f1mu+sigma);
//        
//       boolean f3b= DataFitter.fit(f3, h,"Q"); //No options uses error for sigma 
//        
//        int cnt =0;
//        if(f3b && f3.getParameter(0)>0 && (f3.getChiSquare()/(double)f3.getNDF())<CHI2CUT && Math.abs(f3.getParameter(1)-f1mu)<4) {
//            double mu = f3.getParameter(1);
//            double emu ;
//            if(f3.parameter(1).error()<0 || Double.isNaN(f3.parameter(1).error())) {
//                emu = sigma;
//            } else {
//                emu = f3.parameter(1).error();
//            }
//            if(x<0.05) emu = Math.max(mu, emu);
//            if(x/dmax>0.99) emu = sigma;
//            ge.addPoint(x,  mu, 0, emu);
//            cnt++;
//        } else {
//            
//            //f3.setParameter(0, amp/3.0);
//            f3.setParameter(1, f1mu);
//            f3.setParameter(2, sigma);
//            //f3.setParameter(3, amp);
//            f3.setParameter(4, sigma);
//            f3.setParLimits(0, 0, amp+10000);
//            f3.setParLimits(2, 0, 100);
//            f3.setParLimits(3, 0, amp+10000);
//            f3.setParLimits(4, 0, 100);
//            
//            //DataFitter.fit(f1, h,"Q"); //No options uses error for sigma 
//            //f1mu = f1.getParameter(1);
//            
//            f3.setRange(min, max1);
//            
//            f3b = DataFitter.fit(f3, h,"Q"); //No options uses error for sigma 
//            
//            if(f3b && f3.getParameter(0)>0 && (f3.getChiSquare()/(double)f3.getNDF())<CHI2CUT && Math.abs(f3.getParameter(1)-f1mu)<4) {
//                double mu = f3.getParameter(1);
//                double emu ;
//                if(f3.parameter(1).error()<0 || Double.isNaN(f3.parameter(1).error())) {
//                     emu = sigma;
//                } else {
//                    emu = f3.parameter(1).error();
//                }
//                if(x<0.05) emu = Math.max(mu, emu);
//                if(x/dmax>0.99) emu = sigma;
//                ge.addPoint(x,  mu, 0, emu); 
//                cnt++;
//            } else {
//                ge.addPoint(x,  mean, 0, sigma);
//                cnt++;
//            }
//        }
//        
//        return cnt;    
//    } 
    private static final double CHI2_CUT = 10.0; 
    private static final double MAX_PARAM = 10000.0;
    private static final double MAX_SIGMA = 200.0;
    private static final double X_EDGE_LOW = 0.05;
    private static final double X_EDGE_HIGH = 0.99;
    private static final double MAX_MU_DIFF = 4.0;

    private static int fitLandau(double x, double dmax, H1F hi, GraphErrors ge) {
        H1F h = hi.histClone("hc");

        double meanh   = h.getDataX(h.getMaximumBin());
        double amp     = h.getBinContent(h.getMaximumBin());
        double sigma   = h.getRMS();
        double mean    = h.getMean();
        //double binSize = h.getDataX(1) - h.getDataX(0);

        double min  = Math.max(0, meanh - 2 * sigma);
        double max1 = meanh + 2 * sigma;
        double max2 = meanh + 3 * sigma;

        f1.setRange(min, max1);
        f1.setParameter(0, amp);
        f1.setParLimits(0, 0, amp + MAX_PARAM);
        f1.setParameter(1, mean);
        double m0 = mean - sigma;
        if(m0<0) m0 = 0;
        f1.setParLimits(1, 0, mean + sigma);
        f1.setParameter(2, sigma);
        f1.setParLimits(2, 0, MAX_SIGMA);
        double E1[] = new double[]{0.1,0.1,1};
        boolean fit1Success = DataFitter.fit(f1, h, "Q", E1);
        if (!fit1Success) return 0;

        double f1mu   = f1.getParameter(1);
        double chi2f1 = f1.getChiSquare();

        f1.setParameter(1, meanh);
        fit1Success = DataFitter.fit(f1, h, "Q", E1);
        if (fit1Success && f1.getChiSquare() < chi2f1) {
            f1mu = f1.getParameter(1);
        }

        boolean added = tryF3Fit(h, x, dmax, ge, amp, sigma, f1mu, min, max2);
        if (!added) {
            added = tryF3Fit(h, x, dmax, ge, amp, sigma, f1mu, min, max1);
        }
        if (!added) {
            ge.addPoint(x, mean, 0, sigma);
            return 1;
        }

        return 1;
    }

    private static boolean tryF3Fit(H1F h, double x, double dmax, GraphErrors ge,
                                    double amp, double sigma, double f1mu, double min, double max) {

        f3.setRange(min, max);
        f3.setParameter(0, amp / 3.0);
        f3.setParameter(1, f1mu);
        f3.setParameter(2, sigma);
        f3.setParameter(3, amp);
        f3.setParameter(4, sigma);

        f3.setParLimits(0, 0, amp + MAX_PARAM);
        f3.setParLimits(1, f1mu-sigma, f1mu+sigma);
        f3.setParLimits(2, 0, MAX_SIGMA);
        f3.setParLimits(3, 0, amp + MAX_PARAM);
        f3.setParLimits(4, 0, MAX_SIGMA);
        double E3[] = new double[]{0.1,0.1,1, 1, 0.1};
        boolean fitSuccess = DataFitter.fit(f3, h, "Q", E3);

        if (fitSuccess &&
            f3.getParameter(0) > 0 &&
            (f3.getChiSquare() / (double) f3.getNDF()) < CHI2_CUT &&
            Math.abs(f3.getParameter(1) - f1mu) < MAX_MU_DIFF) {

            double mu = f3.getParameter(1);
            double emu = getParameterErrorSafe(f3, 1, sigma);

            if (x < X_EDGE_LOW) emu = Math.max(mu, emu);
            if (x / dmax > X_EDGE_HIGH) emu = sigma;

            ge.addPoint(x, mu, 0, emu);
            return true;
        }

        return false;
    }

    private static double getParameterErrorSafe(F1D func, int index, double fallback) {
        double err = func.parameter(index).error();
        if (err < 0 || Double.isNaN(err)) {
            return fallback;
        }
        return err;
    }

   

}
