/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.math.Func1D;
/**
 *
 * @author ziegler
 */
public class FitLine2D extends Func1D{
    public int i;
    public double y;
    private FitFunction2D fc ;
    public FitLine2D() {
        super("fcn", 0.0, 1.0);
        fc = new FitFunction2D();
    }
    public static final int nPars = 12;
    private double[] par = new double[nPars];
    public FitLine2D(String name, int i, int j, MnUserParameters pars) {
        super(name, 0.0, 1.0);
        this.i = i;
        //this.y = org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.WireIneffAnal.Bsq[j];
        
        this.y =0;
        for(int k =0; k<WireIneffAnal.nBins; k++) 
            this.y += HistogramManager.Bhists.get(new Coordinate(i,k, j)).getMean();
        this.y /=(double) WireIneffAnal.nBins;
        fc = new FitFunction2D();
        this.initParameters(pars);
    }

    private void initParameters(MnUserParameters pars) {
        for(int p = 0; p< nPars; p++) {
            par[p] = pars.value(p);
        }
    }
    @Override
    public double evaluate(double x) { 
        double term0 = (par[0] + par[1]*y) / Math.pow(x * x + (par[2] + par[3]*y) * x + par[4], (par[5] + par[6]*y + par[7]*y*y));
        double term1 = (par[8] + par[9]*y + par[10]*y*y)/ Math.pow((1 - x) * (1 - x) + par[11], 2);
        return term0 + term1;
    }
    
}
