/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

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
    private Map<Coordinate, GraphErrors> _ineffsvstrkdocasProf;
    
    public FitFunction() {
        
    }
    public FitFunction(int i, Map<Coordinate, GraphErrors> ineffsvstrkdocasProf) {
        _ineffsvstrkdocasProf = ineffsvstrkdocasProf;
        this.i = i;
    }

    /*     
    public double eval(double x, double[] par) {
        double value = par[0]*(par[1]/Math.pow(x*x + par[2], 2) + par[3]/Math.pow( (1-x) + par[4], 2));
        //gemc: double ddEff = dcc.iScale[SECI][SLI]*(dcc.P1[SECI][SLI]/pow(X*X + dcc.P2[SECI][SLI], 2) + dcc.P3[SECI][SLI]/pow( (1-X) + dcc.P4[SECI][SLI], 2));
	
        return value;
    }
    */
    public double eval(double x, double[] par) {
        // par[0] = global scale
        // par[1], par[2] = amplitude, offset near 0
        // par[3], par[4] = amplitude, offset near 1
        // par[5] = exponent near 0
        // par[6] = exponent near 1

        double term0 = par[1] / Math.pow(x*x + par[2], par[5]);
        double term1 = par[3] / Math.pow((1-x) + par[4], par[6]);
        return par[0] * (term0 + term1);
    }

    @Override
    public double valueOf(double[] par) {
        double chisq = 0;
        double delta = 0;
        if(_ineffsvstrkdocasProf.get(new Coordinate(this.i)).getVectorX().size()>0){ 
            GraphErrors gr = _ineffsvstrkdocasProf.get(new Coordinate(this.i));
            for (int ix =0; ix< gr.getDataSize(0); ix++) {
                double x = gr.getDataX(ix);
                double y = gr.getDataY(ix);
                double err = gr.getDataEY(ix);
                if(err>0) {
                    double f = this.eval(x, par);
                    delta = (y - f) / err; 
                    chisq += delta * delta;
                }
            }
        }
        return chisq;
    }
}
