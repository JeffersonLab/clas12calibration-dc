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
import org.jlab.groot.data.H2F;

/**
 *
 * @author veronique
 */

public class FitManager {
    private final Map<Coordinate, MnUserParameters> fitParams = new HashMap<>();
    private final Map<Coordinate, MnUserParameters> fitParams2d = new HashMap<>();
    private final Map<Coordinate, FitLine> fitLines = new HashMap<>();
    private final Map<Coordinate, FitFunction> fitFunctions = new HashMap<>();
    
    
    private final String[] parNames = {"p0", "p2", "p4", "p5", "p8", "p11"};
    private final double[][] pars1d = { 
                                        {0.3,1.1,5.6,-2.0,-45,1.8},
                                        {0.3,1.1,5.6,-3.0,-45,1.8}, 
                                        {0.2,-0.5,1.5,-6.9,0.83,0.1},
                                        {0.2,-0.5,1.5,-7.5,0.13,0.1},                                        
                                      };

    private final double[] errs = {0.1,0.1,0.1,0.1,0.1,0.1};
    //2d
    private final String[] parNames2d = {"p0", "p1", "p2", "p3", "p4", "p5", "p6","p7", "p8", "p9", "p10", "p11"};
    private final double[] pars2d = {2.0,-0.09, -1.03, 0.23, 1.72, -2.94,1.04, -0.23, -5.72,3.71, -0.75, 1.3};
    private final double[] pars2derrs = {1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, 1,1, 1, 0.1};
    private final WireIneffAnal _wia;
    public void intializeFits() {
        //init gui
        for(int l = 0; l<6; l++) {
            _wia.getCalib().addEntry(0,l+1,0);
            for(int p = 0; p<parNames2d.length; p++) {
                _wia.getCalib().setDoubleValue(0.0, parNames2d[p], 0, l+1, 0);
            }
        }
        _wia.getCalib().fireTableDataChanged();
        
        //2d init
        for(int li = 0; li<2; li++) {
            fitParams2d.put(new Coordinate(li), new MinuitPar());
            for(int p = 0; p<parNames2d.length; p++) {
                fitParams2d.get(new Coordinate(li)).add(parNames2d[p], pars2d[p], pars2derrs[p]);
                _wia.getCalib().setDoubleValue(fitParams2d.get(new Coordinate(li)).value(p), parNames2d[p], 0, li+3, 0);
            }
        }
        _wia.getCalib().fireTableDataChanged();
        
        //1d init
        for(int l = 0; l<2; l++) {
            fitParams.put(new Coordinate(l), new MinuitPar());
            fitParams.put(new Coordinate(l+4), new MinuitPar());
            for(int p = 0; p<pars1d[0].length; p++) {
                fitParams.get(new Coordinate(l)).add(parNames[p], pars1d[l][p], errs[p]); //region 1
                fitParams.get(new Coordinate(l+4)).add(parNames[p], pars1d[l+2][p], errs[p]); //refion 3
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(p), parNames[p], 0, l+1, 0);
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l+4)).value(p), parNames[p], 0, l+5, 0);
            }
        }
        _wia.getCalib().fireTableDataChanged();
    }
  

    public void updateTable() {
        
        //2d fit
        for(int li = 0; li<2; li++) {
            for(int p = 0; p<parNames2d.length; p++) {
                _wia.getCalib().setDoubleValue(fitParams2d.get(new Coordinate(li)).value(p), parNames2d[p], 0, li+3, 0);
            }
        }
        //1d fit
        for(int l = 0; l<2; l++) {
            for(int p = 0; p<parNames.length; p++) {
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l)).value(p), parNames[p], 0, l+1, 0);
                _wia.getCalib().setDoubleValue(fitParams.get(new Coordinate(l+4)).value(p), parNames[p], 0, l+5, 0);
            }
            
        }
        
        _wia.getCalib().fireTableDataChanged(); 
    }
    
    public FitManager(WireIneffAnal wia) {
       _wia = wia;
    }
    
    double[][] R2Ubounds = new double[][]{ 
                                       {1.5, 0.0, -0.5, 0.75, 2.5, -0.5, 3.0, 0.0, 0.0, 5.0, 0.0, 2.0},
                                       {1.5, 0.0, -0.5, 0.75, 2.5, -0.5, 3.0, 0.0, 0.0, 5.0, 0.0, 2.0}  
                                        };   
    double[][] R2Lbounds = new double[][]{
                                       {0.5, -0.5, -1.5, 0.05, 1.0, -5.0, 0.5, -0.8, -6.0, 0.0, -1.2, 0.1},
                                       {0.5, -0.5, -1.5, 0.05, 1.0, -5.0, 0.5, -0.8, -6.0, 0.0, -1.2, 0.1}  
                                        };   
    
    public void run2DFit(int sl, Map<Coordinate, H2F>  g,
        Map<Coordinate, H2F>  gr) {
        Coordinate coord = new Coordinate(sl);
        FitFunction2D ff = new FitFunction2D(sl, g, gr);
        MnUserParameters params = fitParams2d.getOrDefault(coord, new MnUserParameters());
        System.out.println("2D INIT "+sl+" "+params.toString());
        MnScan scan = new MnScan(ff, params, 2);
        
        for(int b = 0; b<R2Ubounds[0].length; b++) {
            scan.setLimits(b, R2Lbounds[sl][b], R2Ubounds[sl][b]);
        }
        for (int iter = 0; iter < 100; iter++) {
            params = scan.minimize().userParameters();
        }
       //if(params!=null) System.err.printf("SCAN RESULTS SL %d, : %s%n", sl + 1,  params);
        MnMigrad migrad = new MnMigrad(ff, params, 2);
        //migrad.setCheckAnalyticalDerivatives(true);

        FunctionMinimum min = null;
        for (int iter = 0; iter < 10; iter++) {
            min = migrad.minimize();
            if(min!=null)
                params = min.userParameters();
        }
        for(int i = 0; i<12; i++) {
            for(int p =0; p< 12; p++) {
                if(p!=i) {
                    migrad.fix(p);
                    min = migrad.minimize();
                    migrad.release(p);
                    if(min!=null)
                        params = min.userParameters();
                }
                if(p==i) {
                    migrad.setLimits(i, 
                            min.userParameters().value(i) - min.userParameters().error(i), 
                            min.userParameters().value(i) + min.userParameters().error(i));
                }
            }
        }    
        for (int iter = 0; iter < 10; iter++) {
            min = migrad.minimize();
            if(min!=null)
                params = min.userParameters();
        }
        if(min!=null) System.err.printf("FIT RESULTS SL %d, : %s%n", sl + 1,  min);
        fitParams2d.put(coord, params);
            
    }
    double[][] Ubounds = new double[][]{
                                       {0.5, 5.0, 50.0, 0.0, 0.0, 5.0}, 
                                       {0.5, 5.0, 50.0, 0.0, 0.0, 5.0},
                                       {0.5, 1.0, 5.0, 0.0, 0.0, 0.5},
                                       {0.5, 1.0, 5.0, 0.0, 0.0, 0.5}  
                                        };   
    double[][] Lbounds = new double[][]{
                                       {0.1, 0.0, 0.0, -20.0, -100.0, 1.0}, 
                                       {0.1, 0.0, 0.0, -20.0, -100.0, 1.0},
                                       {0.02, -1.0, 0.0, -10.0, -0.5, 0.09},
                                       {0.02, -1.0, 0.0, -10.0, -0.5, 0.09}  
                                        };   
    // Runs the fitting procedure for a given superlayer.
    public void runFit(int sl, Map<Coordinate, GraphErrors>  prof) {
        Coordinate coord = new Coordinate(sl);
        FitFunction ff = new FitFunction(sl, prof);
        MnUserParameters params = fitParams.getOrDefault(coord, new MnUserParameters());
        
        MnScan scan = new MnScan(ff, params, 2);
        int slb=sl;
        if(sl>3) {
            slb=sl-2;
        }
        for(int b = 0; b<Ubounds[0].length; b++) {
            scan.setLimits(b, Lbounds[slb][b], Ubounds[slb][b]);
        }
        for (int iter = 0; iter < 100; iter++) {
            params = scan.minimize().userParameters();
        }
        //System.err.printf("SCAN RESULTS SL %d: %s%n", sl + 1,  params);
        MnMigrad migrad = new MnMigrad(ff, params, 1);
        migrad.setCheckAnalyticalDerivatives(true);
        
        FunctionMinimum min = null;
        for (int iter = 0; iter < 10; iter++) {
            min = migrad.minimize();
            if(min!=null)
                params = min.userParameters();
        }
        for(int i = 0; i<6; i++) {
            for(int p =0; p< 6; p++) {
                if(p!=i) {
                    migrad.fix(p);
                    min = migrad.minimize();
                    migrad.release(p);
                    if(min!=null)
                        params = min.userParameters();
                }
                if(p==i) {
                    migrad.setLimits(i, 
                            min.userParameters().value(i) - min.userParameters().error(i), 
                            min.userParameters().value(i) + min.userParameters().error(i));
                }
            }
        }    
        for (int iter = 0; iter < 10; iter++) {
            min = migrad.minimize();
            if(min!=null)
                params = min.userParameters();
        }
        System.err.printf("FIT RESULTS SL %d: %s%n", sl + 1,  params);
        fitParams.put(coord, params);
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
    public MnUserParameters getFitParams2d(int sl) {
        return fitParams2d.get(new Coordinate(sl));
    }
}
