/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.docasmear;

import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.FCNBase;
import org.jlab.groot.data.GraphErrors;

/**
 *
 * @author ziegler
 */
public class FitFunction implements FCNBase{

    private int i;
    private int j;
    public String fcn;
    private Map<Coordinate, GraphErrors> _timeResvstrkdocasProf;
    
    public FitFunction() {
        
    }
    public FitFunction(int i, int j, String fc, Map<Coordinate, GraphErrors> timeResvstrkdocasProf) {
        _timeResvstrkdocasProf = timeResvstrkdocasProf;
        this.i = i;
        this.j = j;
        this.fcn = fc;
    }
         
    public double eval(double x, double[] par) {
        double beta = DocaSmearAnal.Beta.get(new Coordinate(this.i,this.j)).getMean();
        beta=1;
        if(this.fcn.equalsIgnoreCase("fc1")) {
            return smearFcn(x, par, beta);
        } else {
            if(this.fcn.equalsIgnoreCase("fc3")) {
                return smearFcn3(x, par);
            } else {
                return smearFcn2(x, par, beta);
            }
        }
    }
    
    public double smearFcn (double x, double[] par, double beta) {
        
        double scale_factor = par[0];
        double a1           = par[1];
        double a2           = par[2];
        double a3           = par[3];
        double a4           = par[4];
        
        double doca_smearing  = scale_factor * ( ( Math.sqrt (x*x + a1 * beta*beta) - x ) 
                + a2 * Math.sqrt(x) + a3 * beta*beta  / (1 - x +a4)*(1 - x +a4) ) ;

        return doca_smearing;
    }
    public double smearFcn2 (double x, double[] par, double beta) {
        
        double scale_factor = par[0];
        double a1           = par[1];
        double a2           = par[2];
        double a3           = par[3];
        double a4           = par[4];
        
        double doca_smearing  = scale_factor * ( ( Math.sqrt (x*x + a1 * beta*beta) - x ) 
                + a2 * Math.sqrt(x) + a3 * beta*beta  / (1 - x +a4) ) ;

        return doca_smearing;
    }
    
    public double smearFcn3 (double x, double[] par) {
        
        double a0           = par[0];
        double a1           = par[1];
        double a2           = par[2];
        double a3           = par[3];
        double a4           = par[4];
        
        double doca_smearing  = a0 + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x ;

        return doca_smearing;
    }
    
    @Override
    public double valueOf(double[] par) {
        double chisq = 0;
        double delta = 0;
        if(_timeResvstrkdocasProf.get(new Coordinate(this.i, j)).getVectorX().size()>0){ 
            GraphErrors gr = _timeResvstrkdocasProf.get(new Coordinate(this.i, j));
            for (int ix =0; ix< gr.getDataSize(0); ix++) {
                double x = gr.getDataX(ix);
                double res = gr.getDataY(ix);
                double err = gr.getDataEY(ix);
                if(err>0 && x> 0.1) {
                    double smear = this.eval(x, par);
                    delta = (res - smear) / err; 
                    chisq += delta * delta;
                }
            }
        }
        return chisq;
    }
}
