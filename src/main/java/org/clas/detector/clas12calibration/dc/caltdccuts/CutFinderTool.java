/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.caltdccuts;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;


/**
 *
 * @author ziegler
 */

public class CutFinderTool {

    private static final int MIN_PLATEAU_WIDTH = 3; // minimum number of bins to qualify as plateau
    public static double findR1UpperEdge(H1F h, int NSigmas) {
        int minBin=0;
        for(int ii =0; ii<h.getMaximumBin(); ii++) {
            if(h.getDataEY(ii)!=0) {
                minBin=ii;
                break;
            }
        }
        
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+[p0]", 
                h.getDataX(minBin), h.getDataX(h.getDataSize(0)-1)); 

        gausFunc.setParameter(0, h.getMax());
        gausFunc.setParameter(1, h.getDataX(h.getMaximumBin()));
        gausFunc.setParameter(2, h.getRMS());
        gausFunc.setParameter(3, 0);
                
        DataFitter.fit(gausFunc, h, "Q");
        
        double c = gausFunc.getParameter(1) + gausFunc.getParameter(2)*(double)NSigmas;
        
        if (c % 5 == 0) {
            return c;
        }

        // Round up to the next multiple of 5
        return Math.ceil(c / 5) * 5;
    }
    
    public static double findLowerEdge(H1F h, int NSigmas) {
        int minBin=0;
        for(int ii =0; ii<h.getMaximumBin(); ii++) {
            if(h.getDataEY(ii)!=0) {
                minBin=ii;
                break;
            }
        }
        
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])", 
                h.getDataX(minBin), h.getDataX(h.getMaximumBin()+1)); 

        gausFunc.setParameter(0, h.getMax());
        gausFunc.setParameter(1, h.getDataX(h.getMaximumBin()));
        gausFunc.setParameter(2, h.getDataX(h.getMaximumBin()+1)-h.getDataX(h.getMaximumBin()-1));
        DataFitter.fit(gausFunc, h, "Q");
        
        double c = gausFunc.getParameter(1) - Math.abs(gausFunc.getParameter(2))*(double)NSigmas;
        
        if(c<0) {
            return 0;
        }
        
        if (c % 5 == 0) {
            return c;
        }

        // Round up to the next multiple of 5
        return Math.ceil(c / 5) * 5;
    }
    
    public static int[] findPlateau(H1F hist) {
        int minBin=0;
        for(int ii =0; ii<hist.getMaximumBin(); ii++) {
            if(hist.getDataEY(ii)!=0) {
                minBin=ii;
                break;
            }
        }
        F1D f = new F1D("f", "[p0]", 
            hist.getDataX(minBin), hist.getDataX(hist.getDataSize(0)-1)); 
        f.setParameter(0, hist.getDataY(minBin));
        f.setLineColor(0);
        f.setLineWidth(0);
        DataFitter.fit(f, hist, "Q");
        
        int maxBin=0;
        for(int ii =hist.getMaximumBin(); ii<hist.getDataSize(0)-1; ii++) {
            if(hist.getDataY(ii)<f.getParameter(0)) {
                maxBin=ii-1;
            }
        }
        for (int start = hist.getMaximumBin(); start < maxBin - MIN_PLATEAU_WIDTH; start++) {
            for (int end = start + MIN_PLATEAU_WIDTH; end < maxBin; end++) {
                double avg = average(hist, start, end);
                boolean isFlat = true;
                for (int i = start; i <= end; i++) { 
                    if (Math.abs(hist.getDataY(i) - avg) / avg > 10) {
                        isFlat = false;
                        break;
                    }
                }
                if (isFlat) {
                    return new int[]{start, end};
                }
            }
        }
        return null;
    }

    private static double average(H1F h, int start, int end) {
        double sum = 0;
        for (int i = start; i <= end; i++) {
            sum += h.getDataY(i);
        }
        return sum / (end - start + 1);
    }
}
