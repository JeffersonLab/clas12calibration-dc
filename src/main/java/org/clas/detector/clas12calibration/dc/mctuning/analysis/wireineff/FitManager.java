/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import java.util.HashMap;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.calt2d.FitUtility.MinuitPar;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.data.GraphErrors;

/**
 *
 * @author veronique
 */

public class FitManager {
    private final Map<Coordinate, MnUserParameters> fitParams = new HashMap<>();
    private final Map<Coordinate, FitLine> fitLines = new HashMap<>();
    private final Map<Coordinate, FitFunction> fitFunctions = new HashMap<>();
    
    private final String[] parNames = {"p0", "p1", "p2", "p3", "p4", "p5", "p6"};
    private final double[] errs = {0.01,0.001,0.001,0.001,0.001, 0.001, 0.001};
    private final WireIneffAnal _wia;
    public void intializeFits(double[][][]pars) {
        for(int s = 0; s<1;s++) {
            for(int l = 0; l<6; l++) {
                fitParams.put(new Coordinate(l), new MinuitPar());
                fitParams.get(new Coordinate(l)).add(parNames[0], 1.0, errs[0]);
                _wia.getCalib().addEntry(0,l+1,0);
                System.out.println("INIT FOR "+l);
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(0), "p0", 0, l+1, 0);
                for(int p = 0; p<4; p++) {
                    fitParams.get(new Coordinate(l)).add(parNames[p+1], pars[s][l][p], errs[p+1]);
                    String st = "p"+(p+1);
                    _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(p+1), st, 0, l+1, 0);
                }
                fitParams.get(new Coordinate(l)).add(parNames[5], 2, errs[5]);
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(5), "p5", 0, l+1, 0);
                fitParams.get(new Coordinate(l)).add(parNames[6], 4, errs[6]);
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(6), "p6", 0, l+1, 0);
            }
        }
        _wia.getCalib().fireTableDataChanged();
    }

    public void updateTable() {
        for(int s = 0; s<1;s++) {
            for(int l = 0; l<6; l++) {
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(0), "p0", 0, l+1, 0);
                for(int p = 0; p<6; p++) {
                    String st = "p"+(p+1);
                    _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(p+1), st, 0, l+1, 0);
                }
            }
        }
        _wia.getCalib().fireTableDataChanged();
    }
    
    public FitManager(WireIneffAnal wia) {
       _wia = wia;
    }
    
    // Runs the fitting procedure for a given superlayer.
    public void runFit(int sl, Map<Coordinate, GraphErrors>  prof) {
        Coordinate coord = new Coordinate(sl);
        FitFunction ff = new FitFunction(sl, prof);
        MnUserParameters params = fitParams.getOrDefault(coord, new MnUserParameters());
        
        MnScan scan = new MnScan(ff, params, 2);
        scan.fix(0);
        for(int i = 5; i<7; i++) {
            scan.fix(i);
        }
        //scan.setLimits(0, 0.99, 1.01);
        for (int iter = 0; iter < 100; iter++) {
            params = scan.minimize().userParameters();
        }
        MnMigrad migrad = new MnMigrad(ff, params, 1);
        migrad.setCheckAnalyticalDerivatives(true);

        for (int iter = 0; iter < 4; iter++) {
            FunctionMinimum min = migrad.minimize();
            System.err.printf("FIT RESULTS SL %d, iteration %d: %s%n", sl + 1, iter + 1, min);
            params = min.userParameters();
        }

        fitParams.put(coord, params);
        this.updateTable();
        
        FitLine fl = new FitLine("fitSL" + sl, sl, params);
        fl.setLineStyle(4);
        fl.setLineWidth(5);
        fl.setLineColor(8);
        fitLines.put(coord, fl);
    }

    // Gets the fit line for a given superlayer.
    public FitLine getFitLine(int sl) {
        return fitLines.get(new Coordinate(sl));
    }

    // Gets the fit parameters for a given superlayer.
    public MnUserParameters getFitParams(int sl) {
        return fitParams.get(new Coordinate(sl));
    }
}
