/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.FCNBase;
import org.jlab.groot.data.H2F;

/**
 *
 * @author ziegler
 */ 

public class FitFunction2D implements FCNBase {

    private int i;
    private Map<Coordinate, H2F> _ineffsvstrkdocas;
    private Map<Coordinate, H2F> _ineffsvstrkdocasErrs;

    // Default constructor
    public FitFunction2D() {
    }

    // Constructor for Map case
    public FitFunction2D(int i, Map<Coordinate, H2F> ineffsvstrkdocas,
             Map<Coordinate, H2F> ineffsvstrkdocasErrs) {
        _ineffsvstrkdocas = ineffsvstrkdocas;
        _ineffsvstrkdocasErrs = ineffsvstrkdocasErrs;
        this.i = i; 
    }


    public double eval(double x, double y, double[] par) {
        double term0 = (par[0] + par[1]*y) / Math.pow(x * x + (par[2] + par[3]*y) * x + par[4], (par[5] + par[6]*y + par[7]*y*y));
        double term1 = (par[8] + par[9]*y + par[10]*y*y)/ Math.pow((1 - x) * (1 - x) + par[11], 2);
        return term0 + term1;
    }

    @Override
    public double valueOf(double[] par) {
        double chisq = 0;
        if (_ineffsvstrkdocas != null) {
            H2F g = _ineffsvstrkdocas.get(new Coordinate(this.i));
            H2F gr = _ineffsvstrkdocasErrs.get(new Coordinate(this.i));
            if (g != null && gr.getEntries() > 0) {
                chisq = computeChiSq(g, gr, par);
            }
        }
        
        return chisq;
    }

    // Helper function to compute chi-square from a graph
    private double computeChiSq(H2F g, H2F gr, double[] par) {
        double chisq = 0;
        for (int ix = 0; ix < g.getDataSize(0); ix++) {
            for (int iy = 0; iy < g.getDataSize(1); iy++) {
                double x = g.getDataX(ix);
                double y = HistogramManager.Bhists.get(new Coordinate(i,ix,iy)).getMean();
                double z = g.getData(ix, iy);
                double err = gr.getData(ix, iy);
                if (err > 0 && x>0 && x<1 && z>0) {
                    double f = this.eval(x, y, par);
                    
                    double delta = (z - f) / err;
                    chisq += delta * delta;
                }
            }
        }
        return chisq;
    }
}


