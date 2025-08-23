/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.docasmear;

import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.math.Func1D;
/**
 *
 * @author ziegler
 */
public class FitLine extends Func1D{
    public int i;
    public int j;
    public String fcn;
    private FitFunction fc ;
    public FitLine() {
        super("fcn", 0.0, 0.85);
        fc = new FitFunction();
    }
    public static final int nPars = 5;
    private double[] par = new double[nPars];
    public FitLine(String name, int i, int j,  String fcnStg, MnUserParameters pars) {
        super(name, 0.0, 0.95);
        this.i = i;
        this.j = j;
        this.fcn = fcnStg;
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
        
        double beta = DocaSmearAnal.Beta.get(new Coordinate(this.i,this.j)).getMean();
        if(this.fcn.equalsIgnoreCase("fc1")) {
            return fc.smearFcn(x, par, beta);
        } else {
            if(this.fcn.equalsIgnoreCase("fc3")) {
                return fc.smearFcn3(x, par);
            } else {
                return fc.smearFcn2(x, par, beta);
            }
        }
        
    }

    
}
