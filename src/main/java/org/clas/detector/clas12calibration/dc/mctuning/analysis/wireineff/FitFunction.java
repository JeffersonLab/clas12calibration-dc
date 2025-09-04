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

public class FitFunction implements FCNBase {

    private int i;
    private Map<Coordinate, GraphErrors> _ineffsvstrkdocasProf;
    private GraphErrors _graph; // <-- new field for single graph

    // Default constructor
    public FitFunction() {
    }

    // Constructor for Map case
    public FitFunction(int i, Map<Coordinate, GraphErrors> ineffsvstrkdocasProf) {
        _ineffsvstrkdocasProf = ineffsvstrkdocasProf;
        this.i = i; 
    }

    // Constructor for single GraphErrors case
    public FitFunction(GraphErrors graph) {
        this._graph = graph;
    }

    public double eval(double x, double[] par) {
        double term0 = par[0] / Math.pow(x * x + par[1] * x + par[2], par[3]);
        double term1 = par[4]/ Math.pow((1 - x) * (1 - x) + par[5], 2);
        return term0 + term1;
    }

    
    /*
    double term0 = par[1] / Math.pow(x * x + par[2] * x + par[5], par[6]);
        double term1 = par[3] / Math.pow((1 - x) * (1 - x) + par[4], 2);
        return par[0] * (term0 + term1);
    */
    
    @Override
    public double valueOf(double[] par) {
        double chisq = 0;
        // Case 1: using Map<Coordinate, GraphErrors>
        if (_ineffsvstrkdocasProf != null) {
            GraphErrors gr = _ineffsvstrkdocasProf.get(new Coordinate(this.i));
            if (gr != null && gr.getVectorX().size() > 0) {
                chisq = computeChiSq(gr, par);
            }
        }

        // Case 2: using single GraphErrors
        else if (_graph != null && _graph.getVectorX().size() > 0) {
            chisq = computeChiSq(_graph, par); 
        }
        
        return chisq;
    }

    // Helper function to compute chi-square from a graph
    private double computeChiSq(GraphErrors gr, double[] par) {
        double chisq = 0;
        for (int ix = 0; ix < gr.getDataSize(0); ix++) {
            double x = gr.getDataX(ix);
            double y = gr.getDataY(ix);
            double err = gr.getDataEY(ix);
            if (err > 0) {
                double f = this.eval(x, par);
                double delta = (y - f) / err;
                chisq += delta * delta;
            }
        }
        return chisq;
    }
}


