/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import java.util.HashMap;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.AnalysisMonitor;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

/**
 *
 * @author veronique
 **/
public class WireIneffAnal extends AnalysisMonitor {
    
    //public static double[] Bsq = new double[] {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25}; 
    public static double[] Bsq = new double[] {0.5,1.5,2.5,3.5,4.5};
    //public static double[] Bsq = new double[] {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25};
    private final HistogramManager histogramMgr = new HistogramManager();
    private final FitManager fitMgr ;
    private final EventProcessor eventProc = new EventProcessor();
    private OutputManager outputMgr;
    private int[][][] totLayA;   // Total layer hits
    private int[][][] effLayA;   // Effective layer hits
    public static String variation = "default";
    
    public WireIneffAnal(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames("Inefficiency vs TrackDoca");
        this.init(false, "p0:p1:p2:p3:p4:p5:p6:p7:p8:p9:p10:p11"); 
        fitMgr = new FitManager(this);
    }

    public static final int nBins = 50;//100;
    @Override
    public void createHistos() {
        histogramMgr.initialize(6, nBins); // 6 superlayers, nBins bins
        totLayA = new int[6][nBins][Bsq.length];  // Initialize total layer hits array
        effLayA = new int[6][nBins][Bsq.length];  // Initialize effective layer hits array
        DataGroup tr = new DataGroup(6,2);
        for(int i =0; i<6; i++) {
            tr.addDataSet(histogramMgr.getNormalizedHistogram(i), 0);
            if(i==2 || i==3) {
                for(int bbin = 0; bbin<Bsq.length; bbin++) {
                    tr.addDataSet(histogramMgr.getBHistogram(i-2,bbin), 0);
                }
            }
            this.getDataGroup().add(tr, 0, i+1, 0); 
        }
        
        this.getDataGroup().add(tr, 0,0,0);
        this.getCalib().fireTableDataChanged();
    }
    
    @Override
    public void processEvent(DataEvent ev) {
        eventProc.processEvent(ev, totLayA, effLayA,nBins,this.fitMgr);
    }
    final Map<Coordinate, FitLine2D> fitLines2d = new HashMap<>();
    final Map<Coordinate, FitLine2D> fitLines2dcolor = new HashMap<>();
    @Override
    public void analysis() {
        histogramMgr.fillFromArrays(totLayA, effLayA);
        
        // Fit and profile each superlayer
        for (int sl = 0; sl < 2; sl++) {
            fitMgr.runFit(sl, histogramMgr.getNormalizedProfiles()); // Fit the profile data
            fitMgr.runFit(sl+4, histogramMgr.getNormalizedProfiles()); // Fit the profile data
            //outputMgr.writeFitParams(sl, fitMgr.getFitParams(sl)); // Write fit parameters to output
            //plot(sl);
//            if(sl==2 || sl==3) {
//                for(int bbin = 0; bbin<Bsq.length; bbin++) {
//                    H1F hb = histogramMgr.getBHistogram(sl-2, bbin);
//                    if(hb.getIntegral()>0) { System.out.println("OK to fit....");
//                        GraphErrors bprofile = histogramMgr.getBProfile(sl-2, bbin);
//                        MnUserParameters mnp = fitMgr.r2Fit(sl,bprofile);
//                        
//                        for(int p =0; p<6; p++) { 
//                             histogramMgr.getParsvvsB().get(new Coordinate(sl-2,p)).addPoint(Bsq[bbin], mnp.value(p+1), 0, mnp.error(p+1));
//                        }
//                    }
//                }
//            }
        }
        for(int sli = 0; sli<2; sli++) {
            fitMgr.run2DFit(sli, histogramMgr.getN2dhistograms(), 
                histogramMgr.getN2dhistogramserrs());
            MnUserParameters params = fitMgr.getFitParams2d(sli);
            for(int j = 0; j<Bsq.length; j++) {
                fitLines2d.put(new Coordinate(sli, j), 
                        new FitLine2D("f2d"+sli+j, sli, j, params));
                fitLines2d.get(new Coordinate(sli, j)).setLineWidth(4);
                if(j==0) fitLines2d.get(new Coordinate(sli, j)).setLineColor(0);
                fitLines2dcolor.put(new Coordinate(sli, j), 
                        new FitLine2D("f2d"+sli+j, sli, j, params));
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineStyle(3);
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineWidth(4);
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineColor(j+1);
                
            }
            
        } 
        
        fitMgr.updateTable();
        TCanvas can1 = new TCanvas("SL 3 B-field dependence of parameters", 800, 800);
        TCanvas can2 = new TCanvas("SL 4 B-field dependence of parameters", 800, 800);
        
        can1.cd(0);
        can1.draw(histogramMgr.getN2dhistograms().get(new Coordinate(0)));
        
   
        can2.cd(0);
        can2.draw(histogramMgr.getN2dhistograms().get(new Coordinate(1)));
        
      /*  
        for(int p =0; p<6; p++) {
            TCanvas can01 = new TCanvas("SL 3 B-field dependence of parameters", 800, 800);
            can01.cd(0);
            can01.draw(histogramMgr.getParsvvsB().get(new Coordinate(0,p)));
        }
         for(int p =0; p<6; p++) {
            TCanvas can02 = new TCanvas("SL 4 B-field dependence of parameters", 800, 800);
            can02.cd(0);
            can02.draw(histogramMgr.getParsvvsB().get(new Coordinate(1,p)));
        }  
        */
        fitMgr.updateTable();
    }

    public void plot(int sl) {
        // Plot the profiles and the fits for each superlayer
        GraphErrors nprofile = histogramMgr.getNormalizedProfile(sl);
        GraphErrors profile = histogramMgr.getProfile(sl);
        
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").clear();
        if(sl==2 || sl==3) {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").divide(2, 2);
        } else {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").divide(2, 1);
        }
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(0);
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(nprofile);  // Draw the profile
        if(sl<2 || sl>3) {
            FitLine fitLine = fitMgr.getFitLine(sl);
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLine, "same");  // Draw the fit line
        }
        
        if(sl==2 || sl==3) {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(2);
            for(int bbin = 0; bbin<Bsq.length; bbin++) {
                H1F hb = histogramMgr.getBHistogram(sl-2, bbin);
                if(hb.getIntegral()>0) {
                    GraphErrors bprofile = histogramMgr.getBProfile(sl-2, bbin);
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(bprofile, "same");
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLines2d.get(new Coordinate(sl-2, bbin)), "same");
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLines2dcolor.get(new Coordinate(sl-2, bbin)), "same");
                }
            }
        }
            
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(1);
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(profile);  // Draw the profile
        if(sl==2 || sl==3) {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(3);
            for(int bbin = 0; bbin<Bsq.length; bbin++) {
                H1F hb = histogramMgr.getBHistogram(sl-2, bbin);
                if(hb.getIntegral()>0) {
                    GraphErrors bprofile = histogramMgr.getBProfileU(sl-2, bbin);
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(bprofile, "same");
                }
            }
        }
    }

    public void finish() {
        // Finalize and close output
        if (outputMgr != null) {
            outputMgr.close();
        }
    }

    // Set the output manager
    public void setOutputManager(OutputManager outputMgr) {
        this.outputMgr = outputMgr;
    }

    // Renaming of output file with run number and timestamp
    public void renameOutput(int runNumber, int iteration) {
        if (outputMgr != null) {
            outputMgr.renameOutput(runNumber, iteration);
        }
    }
    
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_component = (String) cc.getValueAt(row, 2);
        System.out.println(str_sector + " " + str_layer + " " + str_component);
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);
       int component = Integer.parseInt(str_component);

       if(group.hasItem(sector,layer,component)==true){
           this.plot(layer-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       
   
    }
}
