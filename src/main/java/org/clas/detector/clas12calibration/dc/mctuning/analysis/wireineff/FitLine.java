/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.math.Func1D;
/**
 *
 * @author ziegler
 */
public class FitLine extends Func1D{
    public int i;
    private FitFunction fc ;
    public FitLine() {
        super("fcn", 0.0, 1.0);
        fc = new FitFunction();
    }
    public static final int nPars = 5+2;
    private double[] par = new double[nPars];
    public FitLine(String name, int i, MnUserParameters pars) {
        super(name, 0.0, 1.0);
        this.i = i;
        fc = new FitFunction();
        this.initParameters(pars);
    }

    private void initParameters(MnUserParameters pars) {
        for(int p = 0; p< nPars; p++) {
            par[p] = pars.value(p);
        }
    }
    @Override
    public double evaluate(double x) { 
       
        //return par[0]*(par[1]/Math.pow(x*x + par[2], 2) + par[3]/Math.pow( (1-x) + par[4], 2));

        // par[0] = global scale
        // par[1], par[2] = amplitude, offset near 0
        // par[3], par[4] = amplitude, offset near 1
        // par[5] = exponent near 0
        // par[6] = exponent near 1

        double term0 = par[1] / Math.pow(x*x + par[2], par[5]);
        double term1 = par[3] / Math.pow((1-x) + par[4], par[6]);
        return par[0] * (term0 + term1);
    }
    
}
