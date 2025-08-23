/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.GraphErrors;

import java.util.HashMap;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
/**
 *
 * @author veronique
 */


public class HistogramManager {
    private final Map<Coordinate, H1F> histograms = new HashMap<>();
    private final Map<Coordinate, GraphErrors> profiles = new HashMap<>();

    // Initializes histograms and profiles for each superlayer.
    public void initialize(int superlayers, int nBins) {
        for (int sl = 0; sl < superlayers; sl++) {
            Coordinate coord = new Coordinate(sl);
            H1F hist = new H1F("ineff_vs_trkDoca_SL" + (sl + 1),
                               "Superlayer " + (sl + 1), nBins, 0.0, 1.0);
            hist.setTitleX("Track Doca (cm)");
            hist.setTitleY("Inefficiency (%)");
            getHistograms().put(coord, hist);
            getProfiles().put(coord, new GraphErrors());
        }
    }

    // Fills histograms based on the arrays of total and effective layer hits.
    public void fillFromArrays(int[][] totLayA, int[][] effLayA) {
        for (int sl = 0; sl < totLayA.length; sl++) {
            System.out.println("HISTOGRAM "+sl);
            for (int bin = 0; bin < totLayA[sl].length; bin++) {
                double efficiency = (totLayA[sl][bin] == 0) ? 0.0 :
                                    (effLayA[sl][bin] / (double) totLayA[sl][bin]) * 100.0;
                getHistograms().get(new Coordinate(sl)).setBinContent(bin, 100-efficiency);
                double err = 100*(effLayA[sl][bin]/(double)totLayA[sl][bin])*Math.sqrt(1/(double)effLayA[sl][bin] +1/(double)totLayA[sl][bin]);
                getHistograms().get(new Coordinate(sl)).setBinError(bin, err);
                System.out.println(bin+":"+getHistograms().get(new Coordinate(sl)).getDataY(bin)
                        +"+/-"+getHistograms().get(new Coordinate(sl)).getDataEY(bin));
                getProfiles().get(new Coordinate(sl)).addPoint(getHistograms().get(new Coordinate(sl)).getDataX(bin),
                        getHistograms().get(new Coordinate(sl)).getDataY(bin), 
                        0, getHistograms().get(new Coordinate(sl)).getDataEY(bin));
            }        
        }
    }

    // Accessor for histograms.
    public H1F getHistogram(int sl) {
        return getHistograms().get(new Coordinate(sl));
    }

    // Accessor for profiles.
    public GraphErrors getProfile(int sl) {
        return getProfiles().get(new Coordinate(sl));
    }

    /**
     * @return the histograms
     */
    public Map<Coordinate, H1F> getHistograms() {
        return histograms;
    }

    /**
     * @return the profiles
     */
    public Map<Coordinate, GraphErrors> getProfiles() {
        return profiles;
    }
}
