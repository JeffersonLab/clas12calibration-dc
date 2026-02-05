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
import static org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.SegmentAnalysis.cellSizes;
import static org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.WireIneffAnal.nBins;
import static org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.WireIneffAnal.Bsq;
import org.jlab.groot.data.H2F;
/**
 *
 * @author veronique
 */


public class HistogramManager {
    private final Map<Coordinate, H1F> histograms = new HashMap<>();
    private final Map<Coordinate, GraphErrors> profiles = new HashMap<>();
    private final Map<Coordinate, H1F> nhistograms = new HashMap<>();
    private final Map<Coordinate, GraphErrors> nprofiles = new HashMap<>();
    private final Map<Coordinate, H1F> bhistograms = new HashMap<>();
    private final Map<Coordinate, GraphErrors> bprofiles = new HashMap<>();
    private final Map<Coordinate, GraphErrors> bprofilesu = new HashMap<>();
    private final Map<Coordinate, GraphErrors> parsvvsB = new HashMap<>();
    private final Map<Coordinate, H2F> n2dhistograms = new HashMap<>();
    private final Map<Coordinate, H2F> n2dhistogramserrs = new HashMap<>();
    public static Map<Coordinate, H1F> Bhists = new HashMap<>();
    private final Map<Coordinate, H1F> n1dhistograms = new HashMap<>();
    // Initializes histograms and profiles for each superlayer.
    public void initialize(int superlayers, int nBins) {
        //2d histos:
        int NBbins = Bsq.length;
        double BSqbinWidth = Bsq[0];
        double Bmax = Bsq[NBbins-1]+BSqbinWidth;
        
        for (int s = 0; s < 6; s++) { //sector
            for(int sl =0; sl<6; sl++) { //superlayer
                    Coordinate coord = new Coordinate(s, sl);
                    H1F hsl = new H1F("s"+(s+1)+"sl"+(sl+1), "superlayer"+(sl+1), 6, 0.5, 6.5);
                    hsl.setOptStat(0);
                    hsl.setTitleX("Layer");
                    hsl.setTitleY("Efficiency [%]");
                    hsl.setLineColor(sl+2);
                    hsl.setLineWidth(3);
                    n1dhistograms.put(coord, hsl);
                }
            }
        
        for (int s = 0; s < 2; s++) { //regions 2, rsl 1 and 2
            for(int b =0; b<NBbins; b++) {
                for(int k =0; k< nBins; k++) {
                    Coordinate coord = new Coordinate(s, k, b);
                    H1F nhist = new H1F("br2", "B", 10, Bsq[b]-BSqbinWidth, Bsq[b]+BSqbinWidth);
                    Bhists.put(coord, nhist);
                }
            }
            H2F nhist = new H2F("r2sl1", "Inefficiency as a function of trkDoca and B", nBins, 0.0, 1.0, NBbins, 0, Bmax);
            Coordinate coord = new Coordinate(s);
            getN2dhistograms().put(coord, nhist);
            getN2dhistograms().get(coord).setTitleY("B^2 (T^2)");
            getN2dhistograms().get(coord).setTitleX("Norm. trkDoca");
            H2F nhiste = new H2F("r2sl1", "Ineff. Errs as a function of trkDoca and B", nBins, 0.0, 1.0, NBbins, 0, Bmax);
            getN2dhistogramserrs().put(coord, nhiste);
        }
        //
        for (int s = 0; s < 2; s++) {
            for (int p = 0; p < 6; p++) {
                Coordinate coord = new Coordinate(s, p);
                getParsvvsB().put(coord, new GraphErrors());
                getParsvvsB().get(coord).setMarkerColor(9);
                getParsvvsB().get(coord).setLineColor(9);
                getParsvvsB().get(coord).setTitle("Parameter "+(p+1));
                getParsvvsB().get(coord).setTitleX("B^2 (T^2)");
                getParsvvsB().get(coord).setTitleY("parameter value");
            }
        }
        for (int sl = 0; sl < superlayers; sl++) {
            //double cellSize = 2.0*Constants.getInstance().wpdist[sl];
            
            Coordinate coord = new Coordinate(sl);
            
            getProfiles().put(coord, new GraphErrors());
            getProfiles().get(coord).setMarkerColor(sl+1);
            getProfiles().get(coord).setLineColor(sl+1);
            getProfiles().get(coord).setTitleX("trkDOCA (cm)");
            getProfiles().get(coord).setTitleY("Inefficiency (%)");
            
            H1F nhist = new H1F("ineff_vs_normtrkDoca_SL" + (sl + 1),
                               "Superlayer " + (sl + 1), nBins, 0.0, 1.0);
            getNormalizedHistograms().put(coord, nhist);
            getNormalizedProfiles().put(coord, new GraphErrors());
            getNormalizedProfiles().get(coord).setMarkerColor(sl+1);
            getNormalizedProfiles().get(coord).setLineColor(sl+1);
            getNormalizedProfiles().get(coord).setTitleX("Normalized trkDOCA");
            getNormalizedProfiles().get(coord).setTitleY("Inefficiency (%)");
        }
        for(int rsl=0; rsl<2; rsl++) {
            for(int bbin = 0; bbin<Bsq.length; bbin++) {
                H1F nhist = new H1F("ineff_vs_normtrkDoca_Bbin" + (bbin + 1),
                                   "R2-Superlayer " + (rsl + 1), nBins, 0.0, 1.0);
                getBHistograms().put(new Coordinate(rsl,bbin), nhist);
                getBProfiles().put(new Coordinate(rsl, bbin), new GraphErrors());
                getBProfiles().get(new Coordinate(rsl, bbin)).setMarkerColor(bbin+1);
                getBProfiles().get(new Coordinate(rsl, bbin)).setMarkerStyle(2+rsl);
                getBProfiles().get(new Coordinate(rsl, bbin)).setLineColor(bbin+1);
                getBProfiles().get(new Coordinate(rsl, bbin)).setTitleX("Normalized trkDOCA");
                getBProfiles().get(new Coordinate(rsl, bbin)).setTitleY("Inefficiency (%)");
                getBProfilesU().put(new Coordinate(rsl, bbin), new GraphErrors());
                getBProfilesU().get(new Coordinate(rsl, bbin)).setMarkerStyle(2+rsl);
                getBProfilesU().get(new Coordinate(rsl, bbin)).setMarkerColor(bbin+1);
                getBProfilesU().get(new Coordinate(rsl, bbin)).setLineColor(bbin+1);
                getBProfilesU().get(new Coordinate(rsl, bbin)).setTitleX("trkDOCA (cm)");
                getBProfilesU().get(new Coordinate(rsl, bbin)).setTitleY("Inefficiency (%)");
            }
        }
    }

    // Fills histograms based on the arrays of total and effective layer hits.
    public void fillFromArrays(int[][][] totLayA, int[][][] effLayA,
            int[][][] totLayB, int[][][] effLayB) {
        for (int s = 0; s < 6; s++) {
            for (int sl = 0; sl < 6; sl++) {
                for (int l = 0; l < 6; l++) {
                    if( totLayB[s][sl][l] >0 ) {
                        double beff = (effLayB[s][sl][l] / (double) totLayB[s][sl][l]) * 100.0;
                        double berr = 100*(effLayB[s][sl][l]/(double)totLayB[s][sl][l])*
                              Math.sqrt(1/(double)effLayB[s][sl][l] +1/(double)totLayB[s][sl][l]);
                        getN1dhistograms().get(new Coordinate(s, sl)).setBinContent(l, beff);
                        getN1dhistograms().get(new Coordinate(s, sl)).setBinError(l, berr);
                        
                     
                    }
                }
            }
             
         }
        for (int sl = 0; sl < 6; sl++) {
            for (int bin = 0; bin < nBins; bin++) {
                double efficiency=0;
                double err=0;
                double totBNumerator = 0;
                double totBDenominator = 0;
                for(int bbin = 0; bbin<Bsq.length; bbin++) {
                    if(totLayA[sl][bin][bbin]!=0) {
                        totBNumerator+=effLayA[sl][bin][bbin];
                        totBDenominator+=totLayA[sl][bin][bbin];
                        double beff = (effLayA[sl][bin][bbin] / (double) totLayA[sl][bin][bbin]) * 100.0;
                        double berr = 100*(effLayA[sl][bin][bbin]/(double)totLayA[sl][bin][bbin])*
                                 Math.sqrt(1/(double)effLayA[sl][bin][bbin] +1/(double)totLayA[sl][bin][bbin]);
                        if(beff!=0 && (sl==2 || sl==3)) {
                            getN2dhistograms().get(new Coordinate(sl-2)).setBinContent(bin, bbin, 100-beff);
                            getN2dhistogramserrs().get(new Coordinate(sl-2)).setBinContent(bin, bbin, berr);
                            
                            getBHistograms().get(new Coordinate(sl-2, bbin)).setBinContent(bin, 100-beff);
                            getBHistograms().get(new Coordinate(sl-2, bbin)).setBinError(bin, berr);
                            getBProfiles().get(new Coordinate(sl-2, bbin)).
                                    addPoint(getBHistograms().get(new Coordinate(sl-2, bbin)).getDataX(bin),
                                    getBHistograms().get(new Coordinate(sl-2, bbin)).getDataY(bin), 
                                    0, getBHistograms().get(new Coordinate(sl-2, bbin)).getDataEY(bin));
                            getBProfilesU().get(new Coordinate(sl-2, bbin)).
                                    addPoint(getBHistograms().get(new Coordinate(sl-2, bbin)).getDataX(bin)*cellSizes[sl],
                                    getBHistograms().get(new Coordinate(sl-2, bbin)).getDataY(bin), 
                                    0, getBHistograms().get(new Coordinate(sl-2, bbin)).getDataEY(bin));
                        }
                    }
                }
                if(totBDenominator!=0) {
                    efficiency=100.0*totBNumerator/totBDenominator;
                    err = 100*(totBNumerator/totBDenominator)*
                                Math.sqrt(1/totBNumerator +1/totBDenominator);
                }
                    //efficiency = (totLayA[sl][bin][bbin] ==0) ? 0.0 :
                    //                    (effLayA[sl][bin][bbin] / (double) totLayA[sl][bin][bbin]) * 100.0;
                    //double err = 100*(effLayA[sl][bin][bbin]/(double)totLayA[sl][bin][bbin])*
                    //            Math.sqrt(1/(double)effLayA[sl][bin][bbin] +1/(double)totLayA[sl][bin][bbin]);
                if(efficiency!=0.0) {
                    getNormalizedHistograms().get(new Coordinate(sl)).setBinContent(bin, 100-efficiency);

                    getNormalizedHistograms().get(new Coordinate(sl)).setBinError(bin, err);
                    getNormalizedProfiles().get(new Coordinate(sl)).addPoint(getNormalizedHistograms().get(new Coordinate(sl)).getDataX(bin),
                            getNormalizedHistograms().get(new Coordinate(sl)).getDataY(bin), 
                            0, getNormalizedHistograms().get(new Coordinate(sl)).getDataEY(bin));

                    getProfiles().get(new Coordinate(sl)).addPoint(getNormalizedHistograms().get(new Coordinate(sl)).getDataX(bin)*cellSizes[sl],
                            getNormalizedHistograms().get(new Coordinate(sl)).getDataY(bin), 
                            0, getNormalizedHistograms().get(new Coordinate(sl)).getDataEY(bin));

                }
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
    
    // Accessor for histograms.
    public H1F getNormalizedHistogram(int sl) {
        return getNormalizedHistograms().get(new Coordinate(sl));
    }

    // Accessor for profiles.
    public GraphErrors getNormalizedProfile(int sl) {
        return getNormalizedProfiles().get(new Coordinate(sl));
    }

    /**
     * @return the histograms
     */
    public Map<Coordinate, H1F> getNormalizedHistograms() {
        return nhistograms;
    }

    /**
     * @return the profiles
     */
    public Map<Coordinate, GraphErrors> getNormalizedProfiles() {
        return nprofiles;
    }
    
    // Accessor for histograms.
    public H1F getBHistogram(int s, int b) {
        return getBHistograms().get(new Coordinate(s, b));
    }

    // Accessor for profiles.
    public GraphErrors getBProfile(int s, int b) {
        return getBProfiles().get(new Coordinate(s, b));
    }
    

    /**
     * @return the histograms
     */
    public Map<Coordinate, H1F> getBHistograms() {
        return bhistograms;
    }

    /**
     * @return the profiles
     */
    public Map<Coordinate, GraphErrors> getBProfiles() {
        return bprofiles;
    }
    
    // Accessor for profiles.
    public GraphErrors getBProfileU(int s, int b) {
        return getBProfilesU().get(new Coordinate(s, b));
    }
    /**
     * @return the profiles
     */
    public Map<Coordinate, GraphErrors> getBProfilesU() {
        return bprofilesu;
    }

    /**
     * @return the parsvvsB
     */
    public Map<Coordinate, GraphErrors> getParsvvsB() {
        return parsvvsB;
    }

    /**
     * @return the n2dhistograms
     */
    public Map<Coordinate, H2F> getN2dhistograms() {
        return n2dhistograms;
    }

    /**
     * @return the n2dhistogramserrs
     */
    public Map<Coordinate, H2F> getN2dhistogramserrs() {
        return n2dhistogramserrs;
    }
    
    /**
     * @return the n1dhistograms
     */
    public Map<Coordinate, H1F> getN1dhistograms() {
        return n1dhistograms;
    }

}
